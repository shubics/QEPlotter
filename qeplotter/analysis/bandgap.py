"""
Band gap detection tool.
Extracted verbatim from qep.py.
"""
from pathlib import Path
import numpy as np
import re


def parse_kpoints_crystal_b(kpt_file):
    lines = Path(kpt_file).read_text().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith('K_POINTS'):
            break
    lines = lines[i + 2:]
    weights, labels = [], []
    pat = re.compile(r"!(\S+)")
    for ln in lines:
        ln = ln.strip()
        if not ln:
            break
        parts = ln.split()
        weights.append(int(parts[3]))
        m = pat.search(ln)
        labels.append(m.group(1) if m else f"pt{len(labels)}")
    edges = [0]
    for w in weights:
        edges.append(edges[-1] + w)
    return weights, labels, edges


def parse_bandgnu_blocks(band_file):
    """
    Parses a QE 'bands.dat.gnu' file into K-distance and Energy arrays.
    """
    lines = Path(band_file).read_text().splitlines()
    blocks, curr = [], []
    for ln in lines:
        if not ln.strip():
            if curr:
                blocks.append(curr)
                curr = []
        else:
            curr.append(ln)
    if curr:
        blocks.append(curr)
    nk, nb = len(blocks[0]), len(blocks)
    kdist = np.zeros(nk)
    E = np.zeros((nk, nb))
    for b, blk in enumerate(blocks):
        if len(blk) != nk:
            raise ValueError(f"Band {b} has {len(blk)} pts, expected {nk}")
        for i, ln in enumerate(blk):
            p = ln.split()
            if b == 0:
                kdist[i] = float(p[0])
            E[i, b] = float(p[1])
    return kdist, E


def segment_for_index(idx, edges, labels):
    for i in range(len(edges)):
        if idx == edges[i]:
            return labels[i] if i < len(labels) else f"pt{i}"
        if i < len(edges) - 1 and edges[i] < idx < edges[i + 1]:
            return f"{labels[i]}–{labels[i + 1]}"
    return "?"


def detect_band_gap(band_file, kpt_file, fermi_level=None):
    """
    Analyzes the band structure to detect VBM, CBM, and Band Gap.
    """
    kdist, E = parse_bandgnu_blocks(band_file)
    weights, labels, edges = parse_kpoints_crystal_b(kpt_file)

    if fermi_level is None:
        vbm = np.max(E)
        cbm_candidates = E[E > vbm + 1e-6]
        if cbm_candidates.size == 0:
            return "ERROR: No CBM found above VBM."
        cbm = np.min(cbm_candidates)
        gap = cbm - vbm
        fermi_level = 0.5 * (vbm + cbm)
    else:
        E_rel = E - fermi_level
        vbm = np.max(E[E_rel <= 1e-7])
        cbm = np.min(E[E_rel > 0])
        gap = cbm - vbm

    vbm_idx = np.argwhere(np.isclose(E, vbm)).tolist()[0]
    cbm_idx = np.argwhere(np.isclose(E, cbm)).tolist()[0]
    kv, bv = vbm_idx
    kc, bc = cbm_idx

    is_direct = (kv == kc)

    seg_v = segment_for_index(kv, edges, labels)
    seg_c = segment_for_index(kc, edges, labels)

    result = f"{Path(band_file).stem}: Gap={gap:.3f} eV at Fermi={fermi_level:.3f} eV ({'direct' if is_direct else 'indirect'})\n"
    result += f"  VBM: E={vbm:.3f} eV (kpt {kv}, band {bv}, {seg_v})\n"
    result += f"  CBM: E={cbm:.3f} eV (kpt {kc}, band {bc}, {seg_c})\n"
    if not is_direct:
        result += f"  Indirect: VBM at {seg_v}, CBM at {seg_c}\n"

    print(result)
