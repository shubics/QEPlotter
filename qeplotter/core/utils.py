"""
Core utilities: constants, helper functions shared across modules.
Extracted verbatim from qep.py.
"""
import re
import numpy as np
from ase.data import atomic_numbers

# ---- Constants (from bilayer section of qep.py) ----
BOHR_TO_ANGSTROM = 0.529177210903
PLANAR_TOL = 0.25
SHIFT_TOL  = 0.12
_A1 = np.array([0.5, -np.sqrt(3)/2, 0.0])
_A2 = np.array([0.5,  np.sqrt(3)/2, 0.0])
_A3 = np.array([0.0,   0.0,       1.0])


def strip_number(atom_label):
    """Return the chemical symbol at the start of a QE species label.

    QE permits labels such as ``Fe1``, ``Fe_up`` and ``C-h``. Chemical symbols
    are case-insensitive in QE, so normalize them against ASE's element table.
    """
    prefix = re.match(r"[A-Za-z]{1,2}", str(atom_label))
    if not prefix:
        return str(atom_label)
    candidate = prefix.group(0).capitalize()
    if candidate in atomic_numbers:
        return candidate
    fallback = candidate[0]
    return fallback if fallback in atomic_numbers else str(atom_label)


def cart_from_frac(cell: np.ndarray, frac: np.ndarray) -> np.ndarray:
    return frac @ cell
