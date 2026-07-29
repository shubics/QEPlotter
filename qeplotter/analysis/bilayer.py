"""
Bilayer stacking analysis tool.
Extracted verbatim from qep.py.
"""
import ast
import operator
import re
from collections import Counter
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Union
import numpy as np
from ase.formula import Formula

from qeplotter.core.utils import (
    BOHR_TO_ANGSTROM, PLANAR_TOL, _A1, _A2, _A3,
    cart_from_frac, strip_number,
)


def gather_blocks(text: str) -> List[Tuple[str, List[str]]]:
    blocks = []
    current: List[str] = []
    label = None
    auto = 1
    for line in text.splitlines():
        m = re.match(r'^\s*>>>\s*(\S+)', line)
        if m:
            if current and any(s.strip() for s in current):
                blocks.append((label if label else str(auto), current))
                auto += 1
            label = m.group(1)
            current = []
        else:
            current.append(line)
    if current and any(s.strip() for s in current):
        blocks.append((label if label else str(auto), current))
    return blocks

# ---------- QE PARSER ---------- #
_QE_BINARY_OPERATORS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}
_QE_UNARY_OPERATORS = {ast.UAdd: operator.pos, ast.USub: operator.neg}


def _parse_qe_number(value: str) -> float:
    """Parse a numeric QE expression without using ``eval``.

    QE accepts compact arithmetic such as ``1/3`` and ``2*3^(-1/2)`` in
    coordinates. Only numeric constants and the documented arithmetic
    operators are accepted here.
    """
    normalized = re.sub(r"(?<=\d)[dD](?=[+-]?\d)", "e", value)
    expression = ast.parse(normalized.replace("^", "**"), mode="eval")

    def evaluate(node):
        if isinstance(node, ast.Expression):
            return evaluate(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.BinOp) and type(node.op) in _QE_BINARY_OPERATORS:
            return _QE_BINARY_OPERATORS[type(node.op)](
                evaluate(node.left), evaluate(node.right)
            )
        if isinstance(node, ast.UnaryOp) and type(node.op) in _QE_UNARY_OPERATORS:
            return _QE_UNARY_OPERATORS[type(node.op)](evaluate(node.operand))
        raise ValueError(f"Unsupported QE numeric expression: {value}")

    return float(evaluate(expression))


def _card_option(line: str, card: str):
    """Return a QE card option written as ``{x}``, ``(x)`` or plain ``x``."""
    match = re.search(
        rf"\b{re.escape(card)}\b\s*(?:\{{\s*([A-Za-z_]+)\s*\}}"
        rf"|\(\s*([A-Za-z_]+)\s*\)|([A-Za-z_]+))?",
        line,
        re.I,
    )
    if not match:
        return None
    return next((item.lower() for item in match.groups() if item), None)


def _namelist_value(lines: Sequence[str], name: str):
    pattern = re.compile(
        rf"(?<![A-Za-z0-9_]){re.escape(name)}\s*=\s*([0-9.dDeE+\-/^*()]+)",
        re.I,
    )
    for line in lines:
        match = pattern.search(line.split("!")[0])
        if match:
            return _parse_qe_number(match.group(1))
    return None


def parse_qe_block(lines: Sequence[str]) -> Tuple[np.ndarray, List[str], np.ndarray]:
    ibrav = None
    celldm = {}
    cell_parameters = None
    cell_units = None
    atpos_units = 'alat'
    species, coords = [], []

    for ln in lines:
        m = re.search(r'\bibrav\s*=\s*(-?\d+)', ln, re.I)
        if m: ibrav = int(m.group(1))
        for k in range(1,7):
            mk = re.search(
                rf'celldm\(\s*{k}\s*\)\s*=\s*([0-9.dDeE+\-/^*()]+)',
                ln, re.I)
            if mk:
                celldm[k] = _parse_qe_number(mk.group(1))

    conventional = {
        name: _namelist_value(lines, name)
        for name in ("A", "B", "C", "cosAB", "cosAC", "cosBC")
    }
    has_conventional = conventional["A"] is not None
    if has_conventional and 1 in celldm:
        raise ValueError("Specify either celldm or A/B/C parameters, not both.")
    alat_angstrom = (
        celldm[1] * BOHR_TO_ANGSTROM if 1 in celldm
        else conventional["A"]
    )
    nat_value = _namelist_value(lines, "nat")
    nat = int(nat_value) if nat_value is not None else None

    # CELL_PARAMETERS
    for i, ln in enumerate(lines):
        if "CELL_PARAMETERS" in ln.upper():
            cell_units = _card_option(ln, "CELL_PARAMETERS")
            if i + 3 >= len(lines):
                raise ValueError("CELL_PARAMETERS must contain three vectors.")
            mat = [
                [_parse_qe_number(value) for value in lines[j].split()[:3]]
                for j in range(i + 1, i + 4)
            ]
            cell_parameters = np.array(mat)
            break

    # ATOMIC_POSITIONS
    at_start = None
    for i, ln in enumerate(lines):
        if "ATOMIC_POSITIONS" in ln.upper():
            atpos_units = _card_option(ln, "ATOMIC_POSITIONS") or "alat"
            at_start = i+1
            break

    if at_start is not None:
        for ln in lines[at_start:]:
            t = ln.strip()
            if not t: continue
            upper = t.upper()
            if (upper.startswith('K_POINTS') or upper.startswith('CELL_PARAMETERS') or
                upper.startswith('ATOMIC_SPECIES') or t.startswith('/') or t.startswith('&')):
                break
            if t.startswith(("!", "#")):
                continue
            toks = t.split()
            if len(toks) < 4: continue
            sp, x, y, z = toks[:4]
            species.append(sp)
            coords.append([
                _parse_qe_number(x), _parse_qe_number(y), _parse_qe_number(z)
            ])
            if nat is not None and len(species) >= nat:
                break
    coords = np.array(coords)

    if len(species) == 0:
        return None, [], np.zeros((0,3))

    # cell
    if cell_parameters is not None:
        cell = cell_parameters.copy()
        if cell_units == 'alat' or (cell_units is None and alat_angstrom is not None):
            if alat_angstrom is None:
                raise ValueError("CELL_PARAMETERS alat requires celldm(1) or A.")
            cell *= alat_angstrom
        elif cell_units == 'bohr':
            cell *= BOHR_TO_ANGSTROM
        elif cell_units == 'angstrom':
            pass
        elif cell_units is None:
            # QE's deprecated no-option fallback is bohr when no alat is given.
            cell *= BOHR_TO_ANGSTROM
        else:
            raise ValueError(f"Unsupported CELL_PARAMETERS unit: {cell_units}")
    else:
        if ibrav is None:
            ibrav = _guess_ibrav_from_celldm(celldm)
            if ibrav is None:
                raise ValueError("Neither CELL_PARAMETERS nor a usable ibrav was found.")
        if has_conventional:
            cell = ibrav2cell_from_conventional(ibrav, conventional)
        else:
            cell = ibrav2cell(ibrav, celldm)

    # frac
    if atpos_units in ('crystal', 'crystal_sg'):
        frac = coords
    elif atpos_units == 'alat':
        if alat_angstrom is None:
            # For CELL_PARAMETERS without an explicit alat, QE defines alat from
            # the length of the first cell vector.
            alat_angstrom = float(np.linalg.norm(cell[0]))
        cart = coords * alat_angstrom
        frac = cart @ np.linalg.inv(cell)
    elif atpos_units == 'bohr':
        cart = coords * BOHR_TO_ANGSTROM
        frac = cart @ np.linalg.inv(cell)
    elif atpos_units == 'angstrom':
        cart = coords
        frac = cart @ np.linalg.inv(cell)
    else:
        raise ValueError(f"Unsupported ATOMIC_POSITIONS unit: {atpos_units}")

    return cell, species, frac

def _guess_ibrav_from_celldm(cd: dict) -> Union[int,None]:
    if 1 in cd and 3 in cd and len(cd)==2:
        return 4  # hex
    if 1 in cd and 2 in cd and 3 in cd and len(cd)==3:
        return 8  # simple orthorhombic
    if 1 in cd and len(cd)==1:
        return 1  # cubic
    return None

def _positive_sqrt(value, label):
    if value < -1e-12:
        raise ValueError(f"Invalid lattice parameters: {label}² is negative.")
    return np.sqrt(max(0.0, value))


def _ibrav_cell(ibrav, a, b, c, cos_ab=0.0, cos_ac=0.0, cos_bc=0.0):
    """Build QE primitive vectors in Å following INPUT_PW."""
    if ibrav == 0:
        raise ValueError("ibrav=0 requires CELL_PARAMETERS.")
    if ibrav == 1:
        return np.diag([a, a, a])
    if ibrav == 2:
        return np.array([[-a/2, 0, a/2], [0, a/2, a/2], [-a/2, a/2, 0]])
    if ibrav == 3:
        return np.array([[a/2, a/2, a/2], [-a/2, a/2, a/2], [-a/2, -a/2, a/2]])
    if ibrav == -3:
        return np.array([[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]])
    if ibrav == 4:
        return np.array([[a, 0, 0], [-a/2, np.sqrt(3)*a/2, 0], [0, 0, c]])
    if ibrav in (5, -5):
        tx = _positive_sqrt((1 - cos_ab) / 2, "tx")
        ty = _positive_sqrt((1 - cos_ab) / 6, "ty")
        tz = _positive_sqrt((1 + 2*cos_ab) / 3, "tz")
        if ibrav == 5:
            return a * np.array([[tx, -ty, tz], [0, 2*ty, tz], [-tx, -ty, tz]])
        u = tz - 2*np.sqrt(2)*ty
        v = tz + np.sqrt(2)*ty
        return a / np.sqrt(3) * np.array([[u, v, v], [v, u, v], [v, v, u]])
    if ibrav == 6:
        return np.diag([a, a, c])
    if ibrav == 7:
        return np.array([[a/2, -a/2, c/2], [a/2, a/2, c/2], [-a/2, -a/2, c/2]])
    if ibrav == 8:
        return np.diag([a, b, c])
    if ibrav == 9:
        return np.array([[a/2, b/2, 0], [-a/2, b/2, 0], [0, 0, c]])
    if ibrav == -9:
        return np.array([[a/2, -b/2, 0], [a/2, b/2, 0], [0, 0, c]])
    if ibrav == 91:
        return np.array([[a, 0, 0], [0, b/2, -c/2], [0, b/2, c/2]])
    if ibrav == 10:
        return np.array([[a/2, 0, c/2], [a/2, b/2, 0], [0, b/2, c/2]])
    if ibrav == 11:
        return np.array([[a/2, b/2, c/2], [-a/2, b/2, c/2], [-a/2, -b/2, c/2]])
    if ibrav == 12:
        sin_ab = _positive_sqrt(1 - cos_ab**2, "sin(gamma)")
        return np.array([[a, 0, 0], [b*cos_ab, b*sin_ab, 0], [0, 0, c]])
    if ibrav == -12:
        sin_ac = _positive_sqrt(1 - cos_ac**2, "sin(beta)")
        return np.array([[a, 0, 0], [0, b, 0], [c*cos_ac, 0, c*sin_ac]])
    if ibrav == 13:
        sin_ab = _positive_sqrt(1 - cos_ab**2, "sin(gamma)")
        return np.array([[a/2, 0, -c/2], [b*cos_ab, b*sin_ab, 0], [a/2, 0, c/2]])
    if ibrav == -13:
        sin_ac = _positive_sqrt(1 - cos_ac**2, "sin(beta)")
        return np.array([[a/2, b/2, 0], [-a/2, b/2, 0], [c*cos_ac, 0, c*sin_ac]])
    if ibrav == 14:
        sin_ab = _positive_sqrt(1 - cos_ab**2, "sin(gamma)")
        z_squared = (
            1 + 2*cos_bc*cos_ac*cos_ab
            - cos_bc**2 - cos_ac**2 - cos_ab**2
        )
        return np.array([
            [a, 0, 0],
            [b*cos_ab, b*sin_ab, 0],
            [c*cos_ac, c*(cos_bc-cos_ac*cos_ab)/sin_ab,
             c*_positive_sqrt(z_squared, "triclinic z")/sin_ab],
        ])
    raise ValueError(f"Unsupported QE ibrav={ibrav}.")


def ibrav2cell(ibrav: int, cd: dict) -> np.ndarray:
    if 1 not in cd:
        raise ValueError("celldm(1) is required when CELL_PARAMETERS is absent.")
    a = cd[1] * BOHR_TO_ANGSTROM
    b = a * cd.get(2, 1.0)
    c = a * cd.get(3, 1.0)
    return _ibrav_cell(
        ibrav, a, b, c,
        cos_ab=cd.get(6, cd.get(4, 0.0) if ibrav in (5, -5, 12, 13) else 0.0),
        cos_ac=cd.get(5, 0.0),
        cos_bc=cd.get(4, 0.0) if ibrav == 14 else 0.0,
    )


def ibrav2cell_from_conventional(ibrav: int, values: dict) -> np.ndarray:
    a = values["A"]
    if a is None:
        raise ValueError("A is required when conventional lattice parameters are used.")
    b = values["B"] if values["B"] is not None else a
    c = values["C"] if values["C"] is not None else a
    return _ibrav_cell(
        ibrav, a, b, c,
        cos_ab=values["cosAB"] or 0.0,
        cos_ac=values["cosAC"] or 0.0,
        cos_bc=values["cosBC"] or 0.0,
    )


def compute_all_z_distances(cart: np.ndarray):
    for i,j in combinations(range(len(cart)),2):
        yield i,j,abs(cart[j,2]-cart[i,2])

def compute_all_distances(cart: np.ndarray):
    for i,j in combinations(range(len(cart)),2):
        yield i,j,float(np.linalg.norm(cart[j]-cart[i]))

def split_layers(frac: np.ndarray) -> Tuple[List[int], List[int]]:
    if len(frac)==0: return [],[]
    z = frac[:,2]
    order = np.argsort(z)
    gaps  = np.diff(z[order])
    if len(gaps)==0: return list(range(len(frac))), []
    k   = np.argmax(gaps)+1
    thr = (z[order][k-1]+z[order][k])/2
    lower = [i for i,v in enumerate(z) if v<=thr]
    upper = [i for i,v in enumerate(z) if v> thr]
    return lower, upper


def detect_bilayer(cell: np.ndarray, species: List[str], frac: np.ndarray) -> dict:
    """Conservatively decide whether a periodic structure contains two layers.

    The largest cyclic z-gap is treated as vacuum.  A second gap can separate
    two slabs only when the separating gap is physically meaningful and larger
    than either slab's internal corrugation.  Composition is intentionally not
    required to match, so heterobilayers are supported.
    """
    result = {"is_bilayer": False, "lower": [], "upper": [],
              "spacing": None, "reason": ""}
    if len(frac) < 4 or len(species) != len(frac):
        result["reason"] = "At least four labelled atoms are required."
        return result

    # Use fractional c as the layer axis, unwrap immediately after the largest
    # cyclic gap, then measure distances along the actual c vector.
    z = np.mod(np.asarray(frac, dtype=float)[:, 2], 1.0)
    order = np.argsort(z)
    sorted_z = z[order]
    cyclic_gaps = np.append(np.diff(sorted_z), 1.0 + sorted_z[0] - sorted_z[-1])
    vacuum_i = int(np.argmax(cyclic_gaps))
    start = (vacuum_i + 1) % len(z)
    unwrapped_order = np.roll(order, -start)
    unwrapped = np.mod(z[unwrapped_order] - z[unwrapped_order[0]], 1.0)
    internal_gaps = np.diff(unwrapped)
    if not len(internal_gaps):
        result["reason"] = "No separable layers were found."
        return result

    split = int(np.argmax(internal_gaps)) + 1
    first = unwrapped_order[:split].tolist()
    second = unwrapped_order[split:].tolist()
    vacuum_gap = float(cyclic_gaps[vacuum_i])
    layer_gap = float(internal_gaps[split - 1])
    cell = np.asarray(cell, dtype=float)
    normal = np.cross(cell[0], cell[1])
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm < 1e-8:
        result["reason"] = "The in-plane cell vectors do not define a layer plane."
        return result
    normal /= normal_norm
    cell_height = abs(float(np.dot(cell[2], normal)))
    first_internal = np.diff(unwrapped[:split])
    second_internal = np.diff(unwrapped[split:])
    internal_max = max(float(first_internal.max()) if len(first_internal) else 0.0,
                       float(second_internal.max()) if len(second_internal) else 0.0)

    if (vacuum_gap * cell_height < 2.0 or vacuum_gap < 1.2 * layer_gap or
            layer_gap * cell_height < 1.5 or
            (internal_max and layer_gap < 1.35 * internal_max)):
        result["reason"] = "No clear vacuum-separated two-layer slab was detected."
        return result

    cart = np.asarray(frac) @ cell
    projected = cart @ normal
    first_mean = float(projected[first].mean())
    second_mean = float(projected[second].mean())
    lower, upper = (first, second) if first_mean <= second_mean else (second, first)
    result.update(is_bilayer=True, lower=lower, upper=upper,
                  spacing=abs(second_mean - first_mean), reason="")
    return result

def custom_labeling(species: List[str]) -> List[str]:
    cnt, labels = {}, []
    for s in species:
        cnt[s]=cnt.get(s,0)+1
        labels.append(f"{s}{cnt[s]}")
    return labels

_STACKING_DESCRIPTIONS = {
    "AA": "core over core and interface over interface",
    "AB": "upper interface site over lower core; the opposite sites are hollow",
    "BA": "upper core over lower interface site; the opposite sites are hollow",
    "AA′": "both core/interface pairs are cross-aligned",
    "AB′": "core over core; interface sites are staggered",
    "A′B": "interface over interface; core sites are staggered",
}


def _composition_key(elements: Sequence[str]) -> Tuple[Tuple[str, int], ...]:
    counts = Counter(elements)
    divisor = 0
    for count in counts.values():
        divisor = gcd(divisor, count)
    divisor = divisor or 1
    return tuple(sorted((element, count // divisor)
                        for element, count in counts.items()))


def _formula(elements: Sequence[str]) -> str:
    if not elements:
        return "—"
    return Formula.from_list(list(elements)).reduce()[0].format("metal")


def _layer_planes(
    cell: np.ndarray,
    frac: np.ndarray,
    indices: Sequence[int],
    plane_tolerance: float = 0.45,
) -> List[List[int]]:
    """Group one slab into atomic planes along the real layer normal."""
    cell = np.asarray(cell, dtype=float)
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    heights = np.asarray(frac, dtype=float) @ cell @ normal
    ordered = sorted(indices, key=lambda index: heights[index])
    planes: List[List[int]] = []
    for index in ordered:
        if not planes:
            planes.append([index])
            continue
        previous_height = float(np.mean(heights[planes[-1]]))
        if abs(float(heights[index]) - previous_height) <= plane_tolerance:
            planes[-1].append(index)
        else:
            planes.append([index])
    return planes


def _inplane_distances(
    cell: np.ndarray,
    frac: np.ndarray,
    first: int,
    others: Sequence[int],
) -> np.ndarray:
    delta = np.asarray(frac, dtype=float)[list(others), :2] - frac[first, :2]
    delta -= np.round(delta)
    vectors = delta @ np.asarray(cell, dtype=float)[:2]
    return np.linalg.norm(vectors, axis=1)


def _alignment_fraction(
    cell: np.ndarray,
    frac: np.ndarray,
    first: Sequence[int],
    second: Sequence[int],
    tolerance: float,
) -> float:
    if not first or not second:
        return 0.0
    aligned = sum(
        float(np.min(_inplane_distances(cell, frac, index, second))) <= tolerance
        for index in first
    )
    return aligned / len(first)


def _sets_aligned(
    cell: np.ndarray,
    frac: np.ndarray,
    first: Sequence[int],
    second: Sequence[int],
    tolerance: float,
) -> bool:
    # Checking both directions prevents one accidental coincidence in a large
    # or incommensurate supercell from being called a global high-symmetry site.
    return (
        _alignment_fraction(cell, frac, first, second, tolerance) >= 0.8
        and _alignment_fraction(cell, frac, second, first, tolerance) >= 0.8
    )


def _tmd_roles(
    cell: np.ndarray,
    frac: np.ndarray,
    planes: Sequence[Sequence[int]],
    is_lower: bool,
    tolerance: float,
):
    """Return core and interface planes for a trigonal-prismatic MX2-like layer."""
    if len(planes) != 3:
        return None
    bottom, core, top = [list(plane) for plane in planes]
    if not (len(bottom) == len(core) == len(top)):
        return None
    # In a trigonal-prismatic MX2/Janus monolayer the two surface planes have
    # the same in-plane projection. Octahedral or more complex slabs are left
    # unlabelled instead of being forced into the AA/AB vocabulary.
    if not _sets_aligned(cell, frac, bottom, top, tolerance):
        return None
    return {
        "core": core,
        "interface": top if is_lower else bottom,
        "bottom": bottom,
        "top": top,
    }


def _representative_shift(
    cell: np.ndarray,
    frac: np.ndarray,
    lower_core: Sequence[int],
    upper_core: Sequence[int],
) -> Tuple[float, float]:
    shifts = []
    for upper_index in upper_core:
        delta = frac[upper_index, :2] - frac[list(lower_core), :2]
        delta -= np.round(delta)
        distances = np.linalg.norm(delta @ np.asarray(cell)[:2], axis=1)
        shifts.append(delta[int(np.argmin(distances))])
    if not shifts:
        return (float("nan"), float("nan"))
    shift = np.median(np.asarray(shifts), axis=0)
    return (float(shift[0]), float(shift[1]))


def _face_label(elements: Sequence[str], indices: Sequence[int]) -> str:
    return _formula([elements[index] for index in indices])


def analyse_stacking(
    cell: np.ndarray,
    species: List[str],
    frac: np.ndarray,
    planar_tolerance: float = PLANAR_TOL,
) -> Dict[str, object]:
    """Describe a bilayer without forcing unsupported structures into AA/AB.

    The six-site convention is ordered from lower to upper layer:
    R-type ``AA/AB/BA`` and H-type ``AA′/AB′/A′B``.  It is applied only to
    commensurate trigonal-prismatic MX2-like layers.  Heterobilayer and Janus
    identity, layer order, and the two interface terminations remain explicit.
    """
    cell = np.asarray(cell, dtype=float)
    frac = np.asarray(frac, dtype=float)
    elements = [strip_number(label) for label in species]
    detected = detect_bilayer(cell, species, frac)
    result: Dict[str, object] = {
        "is_bilayer": bool(detected["is_bilayer"]),
        "label": "Not applicable",
        "family": "—",
        "confidence": "none",
        "spacing": detected.get("spacing"),
        "bilayer_type": "—",
        "lower_formula": "—",
        "upper_formula": "—",
        "interface": "—",
        "shift": None,
        "description": detected.get("reason", ""),
        "reason": detected.get("reason", ""),
    }
    if not detected["is_bilayer"]:
        return result

    lower, upper = detected["lower"], detected["upper"]
    lower_elements = [elements[index] for index in lower]
    upper_elements = [elements[index] for index in upper]
    lower_planes = _layer_planes(cell, frac, lower)
    upper_planes = _layer_planes(cell, frac, upper)
    lower_bottom, lower_top = lower_planes[0], lower_planes[-1]
    upper_bottom, upper_top = upper_planes[0], upper_planes[-1]
    lower_janus = set(elements[index] for index in lower_bottom) != set(
        elements[index] for index in lower_top)
    upper_janus = set(elements[index] for index in upper_bottom) != set(
        elements[index] for index in upper_top)
    same_composition = _composition_key(lower_elements) == _composition_key(
        upper_elements)
    if lower_janus or upper_janus:
        bilayer_type = "Janus homobilayer" if same_composition else "Janus heterobilayer"
    else:
        bilayer_type = "Homobilayer" if same_composition else "Heterobilayer"

    lower_formula = _formula(lower_elements)
    upper_formula = _formula(upper_elements)
    interface = (
        f"{_face_label(elements, lower_top)} | "
        f"{_face_label(elements, upper_bottom)}"
    )
    result.update(
        bilayer_type=bilayer_type,
        lower_formula=lower_formula,
        upper_formula=upper_formula,
        interface=interface,
        description="Bilayer detected; canonical TMD registry is being checked.",
        reason="",
    )

    lower_roles = _tmd_roles(
        cell, frac, lower_planes, is_lower=True, tolerance=planar_tolerance)
    upper_roles = _tmd_roles(
        cell, frac, upper_planes, is_lower=False, tolerance=planar_tolerance)
    if lower_roles is None or upper_roles is None:
        result.update(
            label="General registry",
            family="Custom / unsupported lattice",
            confidence="safe fallback",
            description=(
                "The bilayer is valid, but it is not a commensurate "
                "trigonal-prismatic three-plane MX2-like pair. No AA/AB label "
                "was assigned."
            ),
            reason="Canonical six-site TMD convention is not applicable.",
        )
        return result

    lower_core = lower_roles["core"]
    lower_face = lower_roles["interface"]
    upper_core = upper_roles["core"]
    upper_face = upper_roles["interface"]
    shift = _representative_shift(cell, frac, lower_core, upper_core)
    result["shift"] = shift

    mm = _sets_aligned(
        cell, frac, upper_core, lower_core, planar_tolerance)
    xx = _sets_aligned(
        cell, frac, upper_face, lower_face, planar_tolerance)
    # Naming convention is ordered lower -> upper. This makes AB and BA
    # deterministic even when exchanging the two layers is not a symmetry.
    upper_core_lower_face = _sets_aligned(
        cell, frac, upper_core, lower_face, planar_tolerance)
    upper_face_lower_core = _sets_aligned(
        cell, frac, upper_face, lower_core, planar_tolerance)

    flags = (mm, xx, upper_core_lower_face, upper_face_lower_core)
    patterns = {
        (True, True, False, False): ("AA", "R-type (parallel)"),
        (False, False, False, True): ("AB", "R-type (parallel)"),
        (False, False, True, False): ("BA", "R-type (parallel)"),
        (False, False, True, True): ("AA′", "H-type (antiparallel)"),
        (True, False, False, False): ("AB′", "H-type (antiparallel)"),
        (False, True, False, False): ("A′B", "H-type (antiparallel)"),
    }
    match = patterns.get(flags)
    if match is None:
        result.update(
            label="General registry",
            family="Twisted / translated / incommensurate",
            confidence="safe fallback",
            description=(
                "No global six-site high-symmetry registry matched. The "
                "fractional core shift is reported without inventing an AA label."
            ),
            reason="No canonical high-symmetry coincidence pattern matched.",
        )
        return result

    label, family = match
    result.update(
        label=label,
        family=family,
        confidence="high",
        description=_STACKING_DESCRIPTIONS[label],
        reason="",
    )
    return result


def classify_stacking(cell: np.ndarray, species: List[str], frac: np.ndarray) -> str:
    """Backward-compatible short stacking label."""
    return str(analyse_stacking(cell, species, frac)["label"])


def analyse_file(path: Union[str, Path]):
    """
    Analyzes a Quantum ESPRESSO input/output file structure.
    Determines lattice parameters, stacking sequence (for 2D), and interlayer distances.
    """
    text = Path(path).read_text()
    blocks = gather_blocks(text)
    for tag, blk in blocks:
        try:
            cell, species, frac = parse_qe_block(blk)
            if len(species)==0:
                print(f"[{tag}] EMPTY BLOCK, SKIPPED.")
                continue
        except Exception as e:
            print(f"[{tag}] ERROR: {e}")
            continue

        cart = cart_from_frac(cell, frac)
        a, c = np.linalg.norm(cell[0]), np.linalg.norm(cell[2])
        stacking = classify_stacking(cell, species, frac)
        labels   = custom_labeling(species)

        print("="*40)
        print(f"[{tag}]  a={a:.3f} Å  c={c:.3f} Å  → stacking: {stacking}")
        print("-"*40)
        print("All ΔZ (vertical) distances (Å):")
        for i,j,dz in compute_all_z_distances(cart):
            print(f"  {labels[i]}-{labels[j]}: {dz:.3f}")
        print("\nAll 3D distances (Å):")
        for i,j,d in compute_all_distances(cart):
            print(f"  {labels[i]}-{labels[j]}: {d:.3f}")
        print()
