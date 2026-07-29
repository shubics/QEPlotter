"""
Band gap detection tool.
Extracted verbatim from qep.py.
"""
from pathlib import Path
import numpy as np
import re


def parse_kpoints_crystal_b(kpt_file):
    lines = Path(kpt_file).read_text().splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith('K_POINTS'):
            start = i + 2
            break
    if start is None:
        raise ValueError("K_POINTS header not found")
    try:
        declared = int(lines[start - 1].split()[0])
    except (ValueError, IndexError):
        declared = None
    entries = []
    pat = re.compile(r"!(\S+)")
    for ln in lines[start:]:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split()
        if len(parts) < 4:
            continue
        m = pat.search(ln)
        entries.append((int(parts[3]), m.group(1) if m else f"pt{len(entries)}"))
        if declared is not None and len(entries) >= declared:
            break
    if len(entries) < 2:
        raise ValueError("crystal_b path needs at least two vertices")
    weights = [weight for weight, _label in entries[:-1]]
    labels = [label for _weight, label in entries]
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
    if labels and edges and idx == edges[-1] - 1:
        return labels[-1]
    for i in range(len(edges)):
        if idx == edges[i]:
            return labels[i] if i < len(labels) else f"pt{i}"
        if i < len(edges) - 1 and edges[i] < idx < edges[i + 1]:
            return f"{labels[i]}–{labels[i + 1]}"
    return "?"


def _bands_cross_fermi(energies, fermi_level, band_axis, tolerance=1e-7):
    """Return True when at least one sampled band spans the Fermi level."""
    energies = np.asarray(energies, dtype=float)
    minimum = np.min(energies, axis=band_axis)
    maximum = np.max(energies, axis=band_axis)
    return bool(np.any(
        (minimum <= fermi_level + tolerance)
        & (maximum >= fermi_level - tolerance)
        & ((minimum < fermi_level - tolerance)
           | (maximum > fermi_level + tolerance))
    ))


def _edge_indices(energies, value, atol=1e-6):
    return np.argwhere(np.isclose(energies, value, atol=atol))


def _is_direct_gap(vbm_indices, cbm_indices, k_column):
    vbm_k = set(int(row[k_column]) for row in vbm_indices)
    cbm_k = set(int(row[k_column]) for row in cbm_indices)
    shared = sorted(vbm_k & cbm_k)
    return bool(shared), shared


def detect_band_gap(band_file, kpt_file, fermi_level=None):
    """
    Analyzes the band structure to detect VBM, CBM, and Band Gap.
    """
    kdist, E = parse_bandgnu_blocks(band_file)
    weights, labels, edges = parse_kpoints_crystal_b(kpt_file)

    if fermi_level is None:
        message = (
            "ERROR: A Fermi level (or occupation data) is required to "
            "distinguish valence and conduction bands."
        )
        print(message)
        return {"metallic": None, "error": message}

    fermi_level = float(fermi_level)
    if _bands_cross_fermi(E, fermi_level, band_axis=0):
        message = (
            f"{Path(band_file).stem}: Metallic — at least one band crosses "
            f"the Fermi level ({fermi_level:.3f} eV)."
        )
        print(message)
        return {"metallic": True, "fermi_level": fermi_level}

    E_rel = E - fermi_level
    below = E[E_rel <= 1e-7]
    above = E[E_rel > 1e-7]
    if below.size == 0 or above.size == 0:
        message = "ERROR: Could not find states on both sides of the Fermi level."
        print(message)
        return {"metallic": None, "error": message}
    vbm = np.max(below)
    cbm = np.min(above)
    gap = cbm - vbm

    vbm_indices = _edge_indices(E, vbm)
    cbm_indices = _edge_indices(E, cbm)
    kv, bv = vbm_indices[0]
    kc, bc = cbm_indices[0]

    is_direct, shared_k = _is_direct_gap(vbm_indices, cbm_indices, k_column=0)
    if shared_k:
        kv = kc = shared_k[0]

    seg_v = segment_for_index(kv, edges, labels)
    seg_c = segment_for_index(kc, edges, labels)

    result = f"{Path(band_file).stem}: Gap={gap:.3f} eV at Fermi={fermi_level:.3f} eV ({'direct' if is_direct else 'indirect'})\n"
    result += f"  VBM: E={vbm:.3f} eV (kpt {kv}, band {bv}, {seg_v})\n"
    result += f"  CBM: E={cbm:.3f} eV (kpt {kc}, band {bc}, {seg_c})\n"
    if not is_direct:
        result += f"  Indirect: VBM at {seg_v}, CBM at {seg_c}\n"

    print(result)
    return {
        "metallic": False,
        "gap": float(gap),
        "fermi_level": fermi_level,
        "vbm_e": float(vbm),
        "cbm_e": float(cbm),
        "vbm_k": int(kv),
        "cbm_k": int(kc),
        "is_direct": is_direct,
    }


# ==============================
# BAND GAP ANNOTATION (for plotting)
# Extracted verbatim from qep.py (_parse_scf_gap, _find_band_gap, _annotate_band_gap).
# ==============================

def _parse_scf_gap(scf_file):
    """
    Parse HOMO/LUMO from QE scf.out file.

    Returns
    -------
    tuple (homo, lumo) or None
    """
    if scf_file is None:
        return None
    try:
        with open(scf_file, 'r', errors='ignore') as f:
            content = f.read()

        # Metal: "the Fermi energy is X eV"
        m = re.search(r"the Fermi energy is\s+([-+]?\d*\.?\d+)\s+eV", content)
        if m:
            ef = float(m.group(1))
            return (ef, ef)  # no gap for metals, but return for reference

        # Insulator: "highest occupied, lowest unoccupied level (ev): HOMO LUMO"
        m2 = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)", content)
        if m2:
            return (float(m2.group(1)), float(m2.group(2)))

        # "highest occupied level (ev): X"
        m3 = re.search(r"highest occupied level\s*\(ev\):\s+([-+]?\d*\.?\d+)", content)
        if m3:
            return (float(m3.group(1)), None)

        return None
    except Exception:
        return None


def _find_band_gap(x_dist, band_energies, fermi_level=None, shift_fermi=False, scf_file=None):
    """
    Detect VBM and CBM from band data.

    If scf_file is provided, uses HOMO/LUMO from SCF output for more accurate detection.
    Otherwise falls back to Fermi-level based detection.

    Parameters
    ----------
    x_dist : np.ndarray, shape (N_k,)
    band_energies : np.ndarray, shape (N_bands, N_k)
        Already Fermi-shifted if shift_fermi was applied.
    fermi_level : float or None
    shift_fermi : bool
    scf_file : str or None
        Path to scf.out file for accurate HOMO/LUMO reference.

    Returns
    -------
    dict or None
    """
    scf_data = _parse_scf_gap(scf_file) if scf_file else None

    if scf_data is not None:
        homo_abs, lumo_abs = scf_data
        if lumo_abs is None or homo_abs == lumo_abs:
            # Metal or no LUMO — no gap
            print(f"SCF: metallic or no gap (HOMO={homo_abs})")
            return None

        # Apply Fermi shift to HOMO/LUMO if band_energies is shifted
        if shift_fermi and fermi_level is not None:
            homo = homo_abs - fermi_level
            lumo = lumo_abs - fermi_level
        else:
            homo = homo_abs
            lumo = lumo_abs

        gap = lumo - homo

        # Find k-point in band data closest to HOMO/LUMO energies
        # VBM: max energy <= homo (with tolerance)
        below_mask = band_energies <= homo + 0.05
        above_mask = band_energies >= lumo - 0.05

        if not np.any(below_mask) or not np.any(above_mask):
            print(f"SCF gap={gap:.3f} eV but could not locate VBM/CBM on band path")
            return None

        # Find the closest band energy to HOMO
        vbm_e = np.max(band_energies[below_mask])
        cbm_e = np.min(band_energies[above_mask])

        vbm_idx = _edge_indices(band_energies, vbm_e, atol=1e-4)
        cbm_idx = _edge_indices(band_energies, cbm_e, atol=1e-4)

        if len(vbm_idx) == 0 or len(cbm_idx) == 0:
            return None

        vbm_band, vbm_k = vbm_idx[0]
        cbm_band, cbm_k = cbm_idx[0]
        is_direct, shared_k = _is_direct_gap(
            vbm_idx, cbm_idx, k_column=1
        )
        if shared_k:
            vbm_k = cbm_k = shared_k[0]

        # Use the band structure VBM/CBM for the gap value (on the k-path)
        actual_gap = cbm_e - vbm_e

        print(f"SCF reference: HOMO={homo_abs:.4f}, LUMO={lumo_abs:.4f}, SCF gap={gap:.3f} eV")

        return {
            'gap': actual_gap,
            'vbm_e': vbm_e,
            'cbm_e': cbm_e,
            'vbm_x': x_dist[vbm_k],
            'cbm_x': x_dist[cbm_k],
            'vbm_k': int(vbm_k),
            'cbm_k': int(cbm_k),
            'is_direct': is_direct,
        }

    # Fallback: Fermi-level based detection
    if fermi_level is None:
        return None

    # Reference energy: 0 if shifted, fermi_level if not
    e_ref = 0.0 if shift_fermi else fermi_level

    if _bands_cross_fermi(
        band_energies, e_ref, band_axis=1
    ):
        print("Metallic band crossing detected on the sampled path.")
        return None

    # All energies below or at Fermi
    below_mask = band_energies <= e_ref + 1e-7
    above_mask = band_energies > e_ref

    if not np.any(below_mask) or not np.any(above_mask):
        return None  # metallic or no states

    vbm_e = np.max(band_energies[below_mask])
    cbm_e = np.min(band_energies[above_mask])
    gap = cbm_e - vbm_e

    if gap <= 0:
        return None  # metallic

    # Find k-point positions of VBM and CBM
    vbm_idx = _edge_indices(band_energies, vbm_e)
    cbm_idx = _edge_indices(band_energies, cbm_e)

    if len(vbm_idx) == 0 or len(cbm_idx) == 0:
        return None

    vbm_band, vbm_k = vbm_idx[0]
    cbm_band, cbm_k = cbm_idx[0]
    is_direct, shared_k = _is_direct_gap(vbm_idx, cbm_idx, k_column=1)
    if shared_k:
        vbm_k = cbm_k = shared_k[0]

    return {
        'gap': gap,
        'vbm_e': vbm_e,
        'cbm_e': cbm_e,
        'vbm_x': x_dist[vbm_k],
        'cbm_x': x_dist[cbm_k],
        'vbm_k': int(vbm_k),
        'cbm_k': int(cbm_k),
        'is_direct': is_direct,
    }


def _annotate_band_gap(ax, gap_info):
    """
    Draw a band gap arrow annotation on the plot axis.

    Parameters
    ----------
    ax : matplotlib Axes
    gap_info : dict from _find_band_gap
    """
    if gap_info is None:
        return

    vbm_e = gap_info['vbm_e']
    cbm_e = gap_info['cbm_e']
    vbm_x = gap_info['vbm_x']
    cbm_x = gap_info['cbm_x']
    gap = gap_info['gap']
    is_direct = gap_info['is_direct']

    # Mark VBM and CBM points
    ax.plot(vbm_x, vbm_e, 'o', color='blue', markersize=6, zorder=10)
    ax.plot(cbm_x, cbm_e, 'o', color='blue', markersize=6, zorder=10)

    if is_direct:
        # Direct gap: vertical double-headed arrow
        ax.annotate(
            '', xy=(vbm_x, cbm_e), xytext=(vbm_x, vbm_e),
            arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5),
            zorder=10
        )
        # Label next to arrow
        mid_e = 0.5 * (vbm_e + cbm_e)
        ax.text(
            vbm_x + 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0]),
            mid_e,
            f'$E_g$ = {gap:.3f} eV',
            fontsize=9, fontweight='bold', color='blue',
            ha='left', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='blue'),
            zorder=11
        )
    else:
        # Indirect gap: angled arrow from VBM to CBM
        ax.annotate(
            '', xy=(cbm_x, cbm_e), xytext=(vbm_x, vbm_e),
            arrowprops=dict(arrowstyle='->', color='blue', lw=1.5, ls='--'),
            zorder=10
        )
        # Label at midpoint
        mid_x = 0.5 * (vbm_x + cbm_x)
        mid_e = 0.5 * (vbm_e + cbm_e)
        ax.text(
            mid_x, mid_e,
            f'$E_g$ = {gap:.3f} eV\n(indirect)',
            fontsize=9, fontweight='bold', color='blue',
            ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='blue'),
            zorder=11
        )

    gap_type = "direct" if is_direct else "indirect"
    print(f"Band gap: {gap:.3f} eV ({gap_type}), VBM at x={vbm_x:.3f} E={vbm_e:.3f}, CBM at x={cbm_x:.3f} E={cbm_e:.3f}")
