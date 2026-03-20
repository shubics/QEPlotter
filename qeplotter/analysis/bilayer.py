"""
Bilayer stacking analysis tool.
Extracted verbatim from qep.py.
"""
import re
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import List, Sequence, Tuple, Union
import numpy as np

from qeplotter.core.utils import BOHR_TO_ANGSTROM, PLANAR_TOL, SHIFT_TOL, _A1, _A2, _A3, cart_from_frac


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
def parse_qe_block(lines: Sequence[str]) -> Tuple[np.ndarray, List[str], np.ndarray]:
    ibrav = None
    celldm = {}
    cell_parameters = None
    cell_units  = 'alat'
    atpos_units = 'alat'
    species, coords = [], []

    for ln in lines:
        m = re.search(r'\bibrav\s*=\s*(-?\d+)', ln, re.I)
        if m: ibrav = int(m.group(1))
        for k in range(1,7):
            mk = re.search(rf'celldm\({k}\)\s*=\s*([0-9.eE+-]+)', ln, re.I)
            if mk: celldm[k] = float(mk.group(1))

    # CELL_PARAMETERS
    for i, ln in enumerate(lines):
        if "CELL_PARAMETERS" in ln.upper():
            m = re.search(r'cell_parameters\s*\{?(\w+)?\}?', ln, re.I)
            if m and m.group(1): cell_units = m.group(1).lower()
            mat = [list(map(float, lines[j].split()[:3])) for j in range(i+1, i+4)]
            cell_parameters = np.array(mat)
            break

    # ATOMIC_POSITIONS
    at_start = None
    for i, ln in enumerate(lines):
        if "ATOMIC_POSITIONS" in ln.upper():
            m = re.search(r'atomic_positions\s*\{?(\w+)?\}?', ln, re.I)
            if m and m.group(1): atpos_units = m.group(1).lower()
            at_start = i+1
            break

    if at_start is not None:
        for ln in lines[at_start:]:
            t = ln.strip()
            if not t: continue
            if (t.startswith('K_POINTS') or t.startswith('CELL_PARAMETERS') or
                t.startswith('ATOMIC_SPECIES') or t.startswith('/') or t.startswith('&')):
                break
            toks = t.split()
            if len(toks) < 4: continue
            sp, x, y, z = toks[:4]
            species.append(sp)
            coords.append([float(x), float(y), float(z)])
    coords = np.array(coords)

    if len(species) == 0:
        return None, [], np.zeros((0,3))

    # cell
    if cell_parameters is not None:
        cell = cell_parameters.copy()
        if cell_units == 'alat':
            if 1 not in celldm:
                raise ValueError("CELL_PARAMETERS {alat} var ama celldm(1) yok!")
            cell *= celldm[1] * BOHR_TO_ANGSTROM
        elif cell_units == 'bohr':
            cell *= BOHR_TO_ANGSTROM
    else:
        if ibrav is None:
            ibrav = _guess_ibrav_from_celldm(celldm)
            if ibrav is None:
                raise ValueError("Ne CELL_PARAMETERS var ne ibrav! (celldm kombinasyonu da yetersiz)")
        cell = ibrav2cell(ibrav, celldm)

    # frac
    if atpos_units in ('crystal','crystal_sg','alat'):
        frac = coords
    elif atpos_units == 'bohr':
        cart = coords * BOHR_TO_ANGSTROM
        frac = cart @ np.linalg.inv(cell)
    elif atpos_units == 'angstrom':
        cart = coords
        frac = cart @ np.linalg.inv(cell)
    elif atpos_units == 'cartesian':
        if cell_parameters is None and '1' in celldm and cell_units == 'alat':
            cart = coords * celldm.get(1,1.0) * BOHR_TO_ANGSTROM
        elif cell_parameters is None and cell_units == 'bohr':
            cart = coords * BOHR_TO_ANGSTROM
        else:
            cart = coords
        frac = cart @ np.linalg.inv(cell)
    else:
        raise ValueError(f"Desteklenmeyen ATOMIC_POSITIONS birimi: {atpos_units}")

    return cell, species, frac

def _guess_ibrav_from_celldm(cd: dict) -> Union[int,None]:
    if 1 in cd and 3 in cd and len(cd)==2:
        return 4  # hex
    if 1 in cd and 2 in cd and 3 in cd and len(cd)==3:
        return 8  # simple orthorhombic
    if 1 in cd and len(cd)==1:
        return 1  # cubic
    return None

def ibrav2cell(ibrav: int, cd: dict) -> np.ndarray:
    a = cd.get(1, 1.0)
    if   ibrav == 1:  cell = np.eye(3)*a
    elif ibrav == 2:  cell = a*np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
    elif ibrav == 3:  cell = a*np.array([[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5]])
    elif ibrav == 4:
        c = a*cd[3]
        cell = np.array([[0.5*a,-np.sqrt(3)/2*a,0],[0.5*a,np.sqrt(3)/2*a,0],[0,0,c]])
    elif ibrav == 5:
        alpha = np.arccos(cd[4])
        v = a*np.array([np.sin(alpha),0,np.cos(alpha)])
        cell = np.vstack([v, np.roll(v,1), np.roll(v,2)])
    elif ibrav == 6:
        cell = np.diag([a,a,a*cd[3]])
    elif ibrav == 7:
        c = a*cd[3]
        cell = np.array([[ a/2,-a/2, c/2],[ a/2, a/2, c/2],[-a/2,-a/2, c/2]])
    elif ibrav == 8:
        cell = np.diag([a,a*cd[2],a*cd[3]])
    elif ibrav == 9:
        b,c = a*cd[2], a*cd[3]
        cell = np.array([[a,0,0],[0,b,0],[a/2,b/2,c]])
    elif ibrav == 10:
        b,c = a*cd[2], a*cd[3]
        cell = np.array([[0,b/2,c/2],[a/2,0,c/2],[a/2,b/2,0]])
    elif ibrav == 11:
        b,c = a*cd[2], a*cd[3]
        cell = np.array([[a/2,b/2,c/2],[-a/2,b/2,c/2],[0,-b/2,c/2]])
    elif ibrav == 12:
        b,c = a*cd[2], a*cd[3]; beta = cd[4]
        cell = np.array([[a,0,0],[0,b,0],[c*np.cos(beta),0,c*np.sin(beta)]])
    elif ibrav == 13:
        b,c = a*cd[2], a*cd[3]; beta = cd[4]
        cell = np.array([[a/2,-b/2,0],[a/2,b/2,0],[c*np.cos(beta),0,c*np.sin(beta)]])
    elif ibrav == 14:
        b,c = a*cd[2], a*cd[3]; cosb,cosa,cosg = cd[4], cd[5], cd[6]
        sing = np.sqrt(1-cosg**2)
        cell = np.array([[a,0,0],
                         [b*cosg, b*sing, 0],
                         [c*cosa, c*(cosb-cosa*cosg)/sing,
                          c*np.sqrt(1-cosa**2-((cosb-cosa*cosg)/sing)**2)]])
    else:
        raise ValueError(f"ibrav={ibrav} tanımlı değil.")
    return cell*BOHR_TO_ANGSTROM


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

def custom_labeling(species: List[str]) -> List[str]:
    cnt, labels = {}, []
    for s in species:
        cnt[s]=cnt.get(s,0)+1
        labels.append(f"{s}{cnt[s]}")
    return labels

def classify_stacking(cell: np.ndarray, species: List[str], frac: np.ndarray) -> str:
    if len(frac)==0: return "NO_ATOMS"
    metal = min(Counter(species), key=species.count)
    lower, upper = split_layers(frac)
    cart = cart_from_frac(cell, frac)
    pairs=set()
    for i in upper:
        for j in lower:
            d = cart[i]-cart[j]; d[2]=0
            if np.linalg.norm(d)<PLANAR_TOL:
                pairs.add((species[i]==metal, species[j]==metal))
    mm = any(pi and pj for pi,pj in pairs)
    xx = any((not pi) and (not pj) for pi,pj in pairs)
    mx = any(pi and (not pj) for pi,pj in pairs)
    xm = any((not pi) and pj for pi,pj in pairs)
    if mx and xm: return 'AA′'
    if xx and not(mm or mx or xm): return 'A′B'
    if xm and not(mm or mx or xx): return 'AB'
    if mm and not(xx or mx or xm): return 'AB′'
    disp=[]
    for i in upper:
        if species[i]!=metal: continue
        dv = frac[lower,:2]-frac[i,:2]; dv -= np.round(dv)
        disp.append(dv[np.argmin(np.linalg.norm(dv,axis=1))])
    if not disp: return 'AA'
    Δ = np.mod(np.mean(disp,axis=0),1)
    CANON={"AA":np.array([0,0]),"AB":np.array([1/3,2/3]),"AB′":np.array([2/3,1/3])}
    name,dist=min(((k,np.linalg.norm(Δ-v)) for k,v in CANON.items()), key=lambda x:x[1])
    return name if dist<SHIFT_TOL else "AA"


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
