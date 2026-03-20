"""
SOC converter: convert projwfc SOC output to (l, ml) basis fatbands.
Extracted verbatim from qep.py (convert_soc_proj_to_ml).
"""
import pathlib
import re
from collections import defaultdict
from typing import Dict, List, Tuple, IO, Generator, Union

_CG_CACHE: Dict[Tuple[int, float, float, int], float] = {}

def _cg_prob(l: int, j: float, mj: float, ml: int) -> float:
    from sympy.physics.wigner import clebsch_gordan as _CG
    from sympy import Rational as _R

    prob = 0.0
    for m_s in (-0.5, +0.5):
        if abs(ml - (mj - m_s)) < 1e-8:
            prob += float(_CG(l, _R(1, 2), j, ml, m_s, mj) ** 2)
    return prob

def _cg_cached(l: int, j: float, mj: float, ml: int) -> float:
    key = (l, j, mj, ml)
    if key not in _CG_CACHE:
        _CG_CACHE[key] = _cg_prob(l, j, mj, ml)
    return _CG_CACHE[key]

_RX_SOC = re.compile(
    r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]+)\s*\)\s*,\s*"
    r"wfc\s*(\d+)\s*\(l=(\d+)\s+j=([0-9./]+)\s+m_j=\s*([0-9.\-+/]+)\)"
)
_RX_COL = re.compile(
    r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]+)\s*\)\s*,\s*"
    r"wfc\s*(\d+)\s*\(l=(\d+)\s*m=\s*(\d+)\)"
)

IdxInfo   = Dict[str, Union[int, float, str]]
GroupName = Tuple[int, str, int, int]

def parse_frac(s):
    if '/' in s:
        num, den = s.split('/')
        return float(num)/float(den)
    return float(s)

def _parse_state_table(text: str) -> Tuple[Dict[int, IdxInfo], bool]:
    idx2info: Dict[int, IdxInfo] = {}
    is_soc = False

    for line in text.splitlines():
        m_soc = _RX_SOC.search(line)
        if m_soc:
            is_soc = True
            gidx, atom, elem, wfc, l, j_str, mj_str = m_soc.groups()

            j  = parse_frac(j_str)
            mj = parse_frac(mj_str)
            gidx, atom, wfc, l = map(int, (gidx, atom, wfc, l))

            idx2info[gidx] = dict(atom=atom, elem=elem.strip(), wfc=wfc, l=l, j=j, mj=mj)
            continue

        m_col = _RX_COL.search(line)
        if m_col:
            gidx, atom, elem, wfc, l, m = m_col.groups()
            gidx, atom, wfc, l, m = map(int, (gidx, atom, wfc, l, m))
            idx2info[gidx] = dict(atom=atom, elem=elem.strip(), wfc=wfc, l=l, m=m)

    return idx2info, is_soc

_K_RX      = re.compile(r"^\s*k\s*=\s*[0-9.eE+\-]+\s+[0-9.eE+\-]+\s+[0-9.eE+\-]+")
_ENERGY_RX = re.compile(r"====\s*e\(\s*\d+\)\s*=\s*([-\d.Ee+]+)\s*eV")
_COEFF_RX  = re.compile(r"([0-9.]+)\*\[\#\s*([0-9]+)\]")

def _stream_states(fh: IO[str]) -> Generator[Tuple[int, float, Dict[int, float]], None, None]:
    ik = 0
    collecting = False
    current_E  = None
    current_w  = defaultdict(float)

    for line in fh:
        if _K_RX.match(line):
            if collecting and current_E is not None:
                yield ik, current_E, current_w
                collecting, current_E = False, None
                current_w = defaultdict(float)
            ik += 1
            continue

        mE = _ENERGY_RX.match(line)
        if mE:
            if collecting and current_E is not None:
                yield ik, current_E, current_w
            current_E  = float(mE.group(1))
            current_w  = defaultdict(float)
            collecting = True
            continue

        if collecting:
            for amp, idx in _COEFF_RX.findall(line):
                current_w[int(idx)] += float(amp) ** 2

    if collecting and current_E is not None:
        yield ik, current_E, current_w

_L2SYM = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}

def _orb_sym(l: int) -> str:
    return _L2SYM.get(l, f"l{l}")

def _make_filename(outdir: pathlib.Path, grp: GroupName) -> pathlib.Path:
    atom, elem, wfc, l = grp
    return outdir / f"fatbands.pdos_atm#{atom}({elem})_wfc#{wfc}({_orb_sym(l)})"

def convert_soc_proj_to_ml(proj_out: Union[str, pathlib.Path],
                           outdir: Union[str, pathlib.Path] = "MLM_pdos",
                           *,
                           overwrite: bool = True,
                           verbose: bool = True) -> None:
    """
    Convert Quantum-ESPRESSO proj.out (with SOC) to (l, ml) basis fatbands.
    """
    proj_out = pathlib.Path(proj_out)
    outdir   = pathlib.Path(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    text = proj_out.read_text(errors="ignore")
    idx2info, is_soc = _parse_state_table(text)

    if not idx2info:
        raise ValueError("Could not parse 'state #' table. Is this a valid projwfc output?")

    orbital_groups: List[GroupName] = []
    seen = set()
    for info in idx2info.values():
        key = (info["atom"], info["elem"], info["wfc"], info["l"])
        if key not in seen:
            seen.add(key)
            orbital_groups.append(key)

    open_files: Dict[GroupName, IO[str]] = {}
    for grp in orbital_groups:
        path = _make_filename(outdir, grp)
        if path.exists() and not overwrite:
            raise FileExistsError(path)
        fh = path.open("w")
        _, _, _, l = grp
        m_cols = "".join(f"   pdos(E)_{i+1}" for i in range(2 * l + 1))
        header = f"# ik    E (eV)   ldos(E){m_cols}\n"
        fh.write(header)
        open_files[grp] = fh

    with proj_out.open("r", errors="ignore") as fh:
        for ik, E, weights in _stream_states(fh):
            ldos_map = defaultdict(float)
            ml_map   = defaultdict(lambda: defaultdict(float))

            for gidx, w in weights.items():
                info = idx2info[gidx]
                grp  = (info["atom"], info["elem"], info["wfc"], info["l"])
                l    = info["l"]

                ldos_map[grp] += w

                if is_soc:
                    j  = info["j"]
                    mj = info["mj"]
                    for ml in range(-l, l + 1):
                        cg2 = _cg_cached(l, j, mj, ml)
                        if cg2 > 1e-8:
                            ml_map[grp][ml] += w * cg2
                else:
                    m = info.get("m", 1)
                    ml_map[grp][m - 1] += w

            for grp in orbital_groups:
                fh_out = open_files[grp]
                l = grp[3]
                ldos = ldos_map[grp]

                row_str = f"{ik:4d}  {E:10.5f}  {ldos:.6e}"
                for ml in range(-l, l + 1):
                    val = ml_map[grp][ml]
                    row_str += f"  {val:.6e}"
                
                fh_out.write(row_str + "\n")

            if verbose and ik % 10 == 0:
                print(f"Processed k-point {ik}", end="\r")

    for fh in open_files.values():
        fh.close()

    if verbose:
        print(f"\nDone - wrote {len(open_files)} (l, m_l)-projected files to '{outdir}'.")
