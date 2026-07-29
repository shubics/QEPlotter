"""
Multi-format structure I/O.

Reads CIF, POSCAR/CONTCAR, Quantum ESPRESSO (``.in`` / ``.out``) and any other
format ASE understands, returning a single common object: an ASE ``Atoms``.

For QE inputs that use ``ibrav``/``celldm`` (which ASE's reader does not always
handle) we fall back to the project's own parser in
:mod:`qeplotter.analysis.bilayer` (``parse_qe_block``), so nothing is duplicated.
"""
import os
import re

import numpy as np
from ase import Atoms
from ase.io import read as ase_read

from qeplotter.core.utils import strip_number

# Human-facing mapping of file extension -> ASE format hint.
SUPPORTED_FORMATS = {
    ".cif": "cif",
    ".poscar": "vasp",
    ".vasp": "vasp",
    ".contcar": "vasp",
    ".xsf": "xsf",
    ".xyz": "extxyz",
    ".pwi": "espresso-in",
    ".in": "espresso-in",
    ".pwo": "espresso-out",
    ".out": "espresso-out",
}

_QE_EXTS = {".in", ".pwi", ".out", ".pwo"}


def _atoms_from_qe_block(path):
    """Fallback QE parser (handles ibrav/celldm) -> ASE Atoms, or None."""
    from qeplotter.analysis.bilayer import gather_blocks, parse_qe_block

    text = open(path, "r", errors="ignore").read()
    blocks = gather_blocks(text)
    # gather_blocks always returns at least one block (the whole file when there
    # are no ">>>" separators), so this also covers plain single-structure files.
    for _tag, blk in blocks:
        try:
            cell, species, frac = parse_qe_block(blk)
        except Exception:
            continue
        if species is None or len(species) == 0:
            continue
        symbols = [strip_number(s) for s in species]
        atoms = Atoms(symbols=symbols, scaled_positions=np.asarray(frac),
                      cell=np.asarray(cell), pbc=True)
        # keep the original QE labels (e.g. "Mo1", "S2") for layer/stacking tools
        atoms.set_array("qe_labels", np.array(species, dtype=object))
        return atoms
    return None


def read_structure(path, fmt=None):
    """
    Read a crystal structure file into an ASE ``Atoms`` object.

    Parameters
    ----------
    path : str
        Path to the structure file (CIF / POSCAR / QE / XSF / XYZ ...).
    fmt : str, optional
        Explicit ASE format. If ``None`` it is guessed from the extension and,
        failing that, ASE's own auto-detection is used.

    Returns
    -------
    ase.Atoms

    Raises
    ------
    ValueError
        If the file cannot be parsed by any available reader.
    """
    ext = os.path.splitext(path)[1].lower()
    guess = fmt or SUPPORTED_FORMATS.get(ext)

    errors = []

    # 1) QE files: try our ibrav-aware parser first (more robust than ASE here).
    if ext in _QE_EXTS:
        try:
            atoms = _atoms_from_qe_block(path)
            if atoms is not None and len(atoms) > 0:
                return atoms
        except Exception as e:  # pragma: no cover - defensive
            errors.append(f"QE parser: {e}")

    # 2) ASE with the guessed format.
    if guess is not None:
        try:
            return ase_read(path, format=guess)
        except Exception as e:
            errors.append(f"ASE[{guess}]: {e}")

    # 3) ASE auto-detection.
    try:
        return ase_read(path)
    except Exception as e:
        errors.append(f"ASE[auto]: {e}")

    raise ValueError(
        f"Could not parse structure file '{os.path.basename(path)}'. "
        f"Tried: {'; '.join(errors)}"
    )


def structure_summary(atoms):
    """
    Return a dict of basic structural information.

    Keys: ``formula``, ``natoms``, ``a``, ``b``, ``c`` (Å),
    ``alpha``, ``beta``, ``gamma`` (deg), ``volume`` (Å³), ``pbc``.
    """
    a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
    return {
        "formula": atoms.get_chemical_formula(empirical=False),
        "reduced_formula": atoms.get_chemical_formula(mode="metal", empirical=True),
        "natoms": len(atoms),
        "a": float(a),
        "b": float(b),
        "c": float(c),
        "alpha": float(alpha),
        "beta": float(beta),
        "gamma": float(gamma),
        "volume": float(atoms.get_volume()),
        "pbc": tuple(bool(p) for p in atoms.pbc),
    }
