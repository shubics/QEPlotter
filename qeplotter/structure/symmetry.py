"""
Space-group / symmetry detection via spglib.

Lightweight (runs in milliseconds on a single structure) so it stays on the
server without burdening it.
"""
import numpy as np
import spglib


# International space-group number ranges -> crystal system.
_CRYSTAL_SYSTEMS = [
    (1, 2, "triclinic"),
    (3, 15, "monoclinic"),
    (16, 74, "orthorhombic"),
    (75, 142, "tetragonal"),
    (143, 167, "trigonal"),
    (168, 194, "hexagonal"),
    (195, 230, "cubic"),
]


def crystal_system(number):
    """Map an international space-group number (1-230) to its crystal system."""
    for lo, hi, name in _CRYSTAL_SYSTEMS:
        if lo <= number <= hi:
            return name
    return "unknown"


def _ds_get(dataset, key, default=None):
    """spglib changed dataset from dict to an object across versions; support both."""
    if dataset is None:
        return default
    if isinstance(dataset, dict):
        return dataset.get(key, default)
    return getattr(dataset, key, default)


def get_spacegroup(atoms, symprec=1e-3, angle_tolerance=-1.0):
    """
    Detect space-group / symmetry information for an ASE ``Atoms``.

    Parameters
    ----------
    atoms : ase.Atoms
    symprec : float
        Distance tolerance (Å) for symmetry finding. Larger = more forgiving.
    angle_tolerance : float
        Angle tolerance in degrees (spglib convention; <0 means automatic).

    Returns
    -------
    dict
        Keys: ``number``, ``international``, ``hall``, ``hall_number``,
        ``pointgroup``, ``crystal_system``, ``n_symmetry_ops``, ``symprec``.
        On failure returns ``{'error': <msg>, 'symprec': symprec}``.
    """
    cell = (
        np.asarray(atoms.cell[:]),
        np.asarray(atoms.get_scaled_positions()),
        np.asarray(atoms.get_atomic_numbers()),
    )
    try:
        ds = spglib.get_symmetry_dataset(
            cell, symprec=symprec, angle_tolerance=angle_tolerance
        )
    except Exception as e:
        return {"error": str(e), "symprec": symprec}

    if ds is None:
        return {"error": "spglib could not determine symmetry", "symprec": symprec}

    number = int(_ds_get(ds, "number", 0))
    rotations = _ds_get(ds, "rotations", None)
    n_ops = int(len(rotations)) if rotations is not None else None

    return {
        "number": number,
        "international": _ds_get(ds, "international", "?"),
        "hall": _ds_get(ds, "hall", "?"),
        "hall_number": _ds_get(ds, "hall_number", None),
        "pointgroup": _ds_get(ds, "pointgroup", "?"),
        "crystal_system": crystal_system(number),
        "n_symmetry_ops": n_ops,
        "symprec": symprec,
    }
