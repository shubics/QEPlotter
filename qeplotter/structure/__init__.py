"""
Generic crystal-structure subpackage.

Format-agnostic structure handling (CIF / POSCAR / Quantum ESPRESSO and more),
built on ASE + spglib. The heavy 3D rendering is delegated to 3Dmol.js, which
runs entirely client-side in the browser (see :mod:`qeplotter.structure.viz`),
so the server only does lightweight parsing/symmetry work.

Public API
----------
read_structure       : read a structure file into an ASE ``Atoms`` object.
structure_summary    : formula, lattice parameters, volume.
get_spacegroup       : spglib space-group / symmetry information.
find_bonds           : PBC-aware bond table (lengths).
find_angles          : PBC-aware bond-angle table.
build_3dmol_html     : client-side 3Dmol.js viewer HTML for an ``Atoms`` object.
"""
from qeplotter.structure.io import read_structure, structure_summary, SUPPORTED_FORMATS
from qeplotter.structure.symmetry import get_spacegroup
from qeplotter.structure.bonds import (
    find_bonds, find_angles, analyse_bonds, analyse_angles, atom_labels,
)
from qeplotter.structure.viz import build_3dmol_html

__all__ = [
    "read_structure",
    "structure_summary",
    "SUPPORTED_FORMATS",
    "get_spacegroup",
    "find_bonds",
    "find_angles",
    "analyse_bonds",
    "analyse_angles",
    "atom_labels",
    "build_3dmol_html",
]
