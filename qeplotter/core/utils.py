"""
Core utilities: constants, helper functions shared across modules.
Extracted verbatim from qep.py.
"""
import re
import numpy as np

# ---- Constants (from bilayer section of qep.py) ----
BOHR_TO_ANGSTROM = 0.529177
PLANAR_TOL = 0.25
SHIFT_TOL  = 0.12
_A1 = np.array([0.5, -np.sqrt(3)/2, 0.0])
_A2 = np.array([0.5,  np.sqrt(3)/2, 0.0])
_A3 = np.array([0.0,   0.0,       1.0])


def strip_number(atom_label):
    return re.sub(r"\d+$", "", atom_label)


def cart_from_frac(cell: np.ndarray, frac: np.ndarray) -> np.ndarray:
    return frac @ cell
