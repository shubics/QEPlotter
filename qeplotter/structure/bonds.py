"""
PBC-aware bond-length and bond-angle analysis.

Uses ASE neighbour lists with covalent-radii based cutoffs, so periodic images
are handled correctly (a bond can cross a cell boundary). Returns tidy pandas
DataFrames ready to drop into a Streamlit table.
"""
import itertools

import numpy as np
import pandas as pd
from ase.neighborlist import NeighborList, natural_cutoffs


def _label(atoms, i):
    """Per-atom label like 'Mo1' (1-based, per element).

    Prefer an explicit QE label only when it carries an index/suffix
    (e.g. 'Mo1', 'Fe_up'); a bare element label ('Mo') is not unique, so we
    fall back to a generated 1-based per-element index.
    """
    sym = atoms[i].symbol
    if "qe_labels" in atoms.arrays:
        lbl = str(atoms.get_array("qe_labels")[i])
        if lbl and lbl != sym:
            return lbl
    same = [k for k, a in enumerate(atoms) if a.symbol == sym]
    return f"{sym}{same.index(i) + 1}"


def atom_labels(atoms):
    """Return the same stable labels used by bond and angle tables."""
    return [_label(atoms, index) for index in range(len(atoms))]


def _neighbor_list(atoms, tol):
    cutoffs = natural_cutoffs(atoms, mult=tol)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
    nl.update(atoms)
    return nl


def _image_label(label, offset):
    """Make periodic copies explicit instead of showing duplicate-looking rows."""
    offset = tuple(int(value) for value in offset)
    if not any(offset):
        return label
    image = ",".join(f"{value:+d}" for value in offset)
    return f"{label} [{image}]"


def analyse_bonds(atoms, tol=1.2, decimals=3, max_rows=2000,
                  include_periodic=True):
    """Return the public bond table together with aligned 3D geometry records."""
    nl = _neighbor_list(atoms, tol)
    pos = atoms.get_positions()
    cell = np.asarray(atoms.get_cell())
    seen, records = set(), []

    for i in range(len(atoms)):
        idx, offsets = nl.get_neighbors(i)
        for j, raw_offset in zip(idx, offsets):
            offset = tuple(int(value) for value in raw_offset)
            if not include_periodic and any(offset):
                continue
            if i < j:
                key = (i, int(j), offset)
            elif j < i:
                key = (int(j), i, tuple(-value for value in offset))
            else:
                inverse = tuple(-value for value in offset)
                key = (i, i, min(offset, inverse))
            if key in seen:
                continue
            seen.add(key)

            start = np.asarray(pos[i], dtype=float)
            end = np.asarray(pos[j] + np.asarray(offset) @ cell, dtype=float)
            distance = float(np.linalg.norm(end - start))
            element_i, element_j = atoms[i].symbol, atoms[j].symbol
            records.append({
                "atom_1": _label(atoms, i),
                "atom_2": _image_label(_label(atoms, j), offset),
                "elements": "-".join(sorted((element_i, element_j))),
                "length (Å)": round(distance, decimals),
                "indices": [int(i), int(j)],
                "point_labels": [_label(atoms, i),
                                 _image_label(_label(atoms, j), offset)],
                "offset": list(offset),
                "points": [start.tolist(), end.tolist()],
                "point_elements": [element_i, element_j],
            })

    records.sort(key=lambda record: (record["length (Å)"], record["atom_1"],
                                     record["atom_2"]))
    records = records[:max_rows]
    columns = ["atom_1", "atom_2", "elements", "length (Å)"]
    frame = pd.DataFrame([{column: record[column] for column in columns}
                          for record in records], columns=columns)
    return frame, records


def find_bonds(atoms, tol=1.2, decimals=3, max_rows=2000,
               include_periodic=True):
    """
    Detect bonds via covalent-radii cutoffs (scaled by ``tol``), PBC-aware.

    Parameters
    ----------
    atoms : ase.Atoms
    tol : float
        Multiplier on summed covalent radii. ~1.1-1.3 is typical.
    decimals : int
        Rounding for the length column.
    max_rows : int
        Safety cap on the number of returned bonds.

    Returns
    -------
    pandas.DataFrame
        Columns: ``atom_1``, ``atom_2``, ``elements``, ``length (Å)``.
        Sorted by length. Each undirected bond appears once.
    """
    return analyse_bonds(atoms, tol=tol, decimals=decimals,
                         max_rows=max_rows,
                         include_periodic=include_periodic)[0]


def analyse_angles(atoms, tol=1.2, decimals=2, max_rows=4000,
                   include_periodic=True):
    """Return the public angle table together with aligned PBC geometry."""
    nl = _neighbor_list(atoms, tol)
    pos = atoms.get_positions()
    cell = np.asarray(atoms.get_cell())
    seen, records = set(), []

    for i in range(len(atoms)):
        idx, offsets = nl.get_neighbors(i)
        neighbours = [(int(j), tuple(int(value) for value in offset),
                       np.asarray(pos[j] + offset @ cell - pos[i], dtype=float))
                      for j, offset in zip(idx, offsets)
                      if include_periodic or not np.any(offset)]
        for (j, offset_j, vector_j), (k, offset_k, vector_k) in itertools.combinations(neighbours, 2):
            endpoint_j = (j, offset_j)
            endpoint_k = (k, offset_k)
            key = (i, min(endpoint_j, endpoint_k), max(endpoint_j, endpoint_k))
            if endpoint_j == endpoint_k or key in seen:
                continue
            seen.add(key)
            norm_j, norm_k = np.linalg.norm(vector_j), np.linalg.norm(vector_k)
            if norm_j < 1e-6 or norm_k < 1e-6:
                continue
            cosine = float(np.clip(np.dot(vector_j, vector_k) /
                                   (norm_j * norm_k), -1.0, 1.0))
            angle = float(np.degrees(np.arccos(cosine)))

            endpoints = sorted([
                (j, offset_j, vector_j, _label(atoms, j), atoms[j].symbol),
                (k, offset_k, vector_k, _label(atoms, k), atoms[k].symbol),
            ], key=lambda item: (item[3], item[1]))
            j, offset_j, vector_j, label_j, symbol_j = endpoints[0]
            k, offset_k, vector_k, label_k, symbol_k = endpoints[1]
            vertex_position = np.asarray(pos[i], dtype=float)
            records.append({
                "vertex": _label(atoms, i),
                "neighbor_1": _image_label(label_j, offset_j),
                "neighbor_2": _image_label(label_k, offset_k),
                "elements": "-".join((symbol_j, atoms[i].symbol, symbol_k)),
                "angle (°)": round(angle, decimals),
                "indices": [int(j), int(i), int(k)],
                "point_labels": [_image_label(label_j, offset_j),
                                 _label(atoms, i),
                                 _image_label(label_k, offset_k)],
                "offsets": [list(offset_j), [0, 0, 0], list(offset_k)],
                "points": [(vertex_position + vector_j).tolist(),
                           vertex_position.tolist(),
                           (vertex_position + vector_k).tolist()],
                "point_elements": [symbol_j, atoms[i].symbol, symbol_k],
            })

    records.sort(key=lambda record: (record["vertex"], record["angle (°)"],
                                     record["neighbor_1"], record["neighbor_2"]))
    records = records[:max_rows]
    columns = ["vertex", "neighbor_1", "neighbor_2", "elements", "angle (°)"]
    frame = pd.DataFrame([{column: record[column] for column in columns}
                          for record in records], columns=columns)
    return frame, records


def find_angles(atoms, tol=1.2, decimals=2, max_rows=4000,
                include_periodic=True):
    """
    Detect bond angles (j-i-k, vertex at i) for all bonded neighbour pairs.

    PBC-aware: each neighbour is taken in the correct periodic image.

    Returns
    -------
    pandas.DataFrame
        Columns: ``vertex``, ``neighbor_1``, ``neighbor_2``, ``elements``,
        ``angle (°)``. Sorted by vertex then angle.
    """
    return analyse_angles(atoms, tol=tol, decimals=decimals,
                          max_rows=max_rows,
                          include_periodic=include_periodic)[0]
