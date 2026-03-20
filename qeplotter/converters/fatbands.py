"""
Converter: standardize projwfc output for fatband plotting.
Extracted verbatim from qep.py (convert_consistent).
"""
import pathlib
import re
from collections import defaultdict
from typing import Dict, List, Tuple, IO, Generator, Union

_STATE_RX = re.compile(
    r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]+)\s*\)"
    r"\s*,\s*wfc\s*(\d+)\s*\(l=(\d+)\s+j=([0-9.]+)\s+m_j=\s*([-\d.]+)\)"
)

IdxInfo   = Dict[str, Union[int, float, str]]
GroupKey  = Tuple[int, str, int, int, float]

def _parse_state_table(text: str) -> Tuple[Dict[int, IdxInfo], Dict[GroupKey, List[int]]]:
    spin_orbit_coupled = "j=" in text

    if spin_orbit_coupled:
        state_rx = re.compile(
            r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]+)\s*\)"
            r"\s*,\s*wfc\s*(\d+)\s*\(l=(\d+)\s+j=([0-9.]+)\s+m_j=\s*([-\d.]+)\)"
        )
    else:
        state_rx = re.compile(
            r"state #\s*(\d+):\s*atom\s*(\d+)\s*\(\s*([A-Za-z]+)\s*\)\s*,\s*"
            r"wfc\s*(\d+)\s*\(l=(\d+)\s*m=\s*\d+\)"
        )

    idx2info: Dict[int, IdxInfo] = {}
    group2idx: Dict[Tuple, List[Tuple[float, int]]] = defaultdict(list)

    for line in text.splitlines():
        m = state_rx.search(line)
        if not m:
            continue

        if spin_orbit_coupled:
            gidx, atom, elem, wfc, l, j, mj = m.groups()
            gidx, atom, wfc, l = map(int, (gidx, atom, wfc, l))
            j, mj = float(j), float(mj)
            info = dict(atom=atom, elem=elem.strip(), wfc=wfc, l=l, j=j, mj=mj)
            key = (atom, elem.strip(), wfc, l, j)
            group2idx[key].append((mj, gidx))
        else:
            gidx, atom, elem, wfc, l = m.groups()
            gidx, atom, wfc, l = map(int, (gidx, atom, wfc, l))
            info = dict(atom=atom, elem=elem.strip(), wfc=wfc, l=l)
            key = (atom, elem.strip(), wfc, l)
            group2idx[key].append((0.0, gidx))

        idx2info[int(gidx)] = info

    group_sorted: Dict[Tuple, List[int]] = {}
    for key, lst in group2idx.items():
        lst.sort(key=lambda t: t[0])
        group_sorted[key] = [g for _, g in lst]

    return idx2info, group_sorted

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


def _make_filename(outdir: pathlib.Path, key: Tuple) -> pathlib.Path:
    atom, elem, wfc, l = key[:4]
    j = key[4] if len(key) > 4 else None

    if j is not None:
        jstr = f"{j:.1f}".rstrip("0").rstrip(".")
        filename = f"fatbands.pdos_atm#{atom}({elem})_wfc#{wfc}({_orb_sym(l)}_j{jstr})"
    else:
        filename = f"fatbands.pdos_atm#{atom}({elem})_wfc#{wfc}({_orb_sym(l)})"

    return outdir / filename


def convert_consistent(proj_out: Union[str, pathlib.Path],
                       outdir: Union[str, pathlib.Path] = "BMS_pdos",
                       *,
                       overwrite: bool = True,
                       verbose: bool = True) -> None:
    """
    Standardizes 'projwfc' output files for plotting.
    """
    proj_out = pathlib.Path(proj_out)
    outdir   = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    text = proj_out.read_text(errors="ignore")
    idx2info, group2idx = _parse_state_table(text)
    all_groups = list(group2idx)

    open_files: Dict[GroupKey, Tuple[IO[str], List[int]]] = {}
    for key in all_groups:
        path = _make_filename(outdir, key)
        if path.exists() and not overwrite:
            raise FileExistsError(path)
        fh = path.open("w")
        idx_list = group2idx[key]
        header = "# ik    E (eV)   ldos(E)" + "".join(
            f"   pdos(E)_{i+1}" for i in range(len(idx_list))
        ) + "\n"
        fh.write(header)
        open_files[key] = (fh, idx_list)

    with proj_out.open("r", errors="ignore") as fh:
        for ik, E, weights in _stream_states(fh):
            ldos = {k: 0.0 for k in all_groups}
            cols = {k: [0.0] * len(group2idx[k]) for k in all_groups}

            for gidx, w in weights.items():
                info = idx2info[gidx]
                if "j" in info:
                    key = (info["atom"], info["elem"], info["wfc"], info["l"], info["j"])
                else:
                    key = (info["atom"], info["elem"], info["wfc"], info["l"])

                col  = group2idx[key].index(gidx)
                ldos[key]         += w
                cols[key][col]     = w

            for key in all_groups:
                fh_out, _ = open_files[key]
                row = f"{ik:4d}  {E:10.5f}  {ldos[key]:.6e}" + "".join(
                      f"  {c:.6e}" for c in cols[key]) + "\n"
                fh_out.write(row)

            if verbose and ik % 10 == 0:
                print(f"Processed k-point {ik}", end="\r")

    for fh, _ in open_files.values():
        fh.close()
    if verbose:
        print(f"\nDone - wrote {len(open_files)} projector files to '{outdir}'.")
