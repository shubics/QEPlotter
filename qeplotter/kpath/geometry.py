"""Reciprocal-space geometry only; contains no band-path policy."""
from itertools import product

import numpy as np
from scipy.spatial import Voronoi


def reciprocal_rows(cell):
    """Return reciprocal row vectors including the 2π factor (Å⁻¹)."""
    return 2.0 * np.pi * np.linalg.inv(np.asarray(cell, dtype=float)).T


def first_brillouin_zone(reciprocal):
    """Construct the reciprocal Wigner--Seitz cell by Voronoi tessellation."""
    shifts = list(product(range(-2, 3), repeat=3))
    points = np.asarray(shifts, dtype=float) @ reciprocal
    origin_index = shifts.index((0, 0, 0))
    voronoi = Voronoi(points)
    faces_raw = []
    for pair, ridge in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        if origin_index not in pair or -1 in ridge or len(ridge) < 3:
            continue
        faces_raw.append([voronoi.vertices[index] for index in ridge])
    if not faces_raw:
        raise ValueError("Could not construct a finite first Brillouin zone")

    vertices, faces, lookup = [], [], {}
    for face in faces_raw:
        indices = []
        for vertex in face:
            key = tuple(np.round(vertex, 9))
            if key not in lookup:
                lookup[key] = len(vertices)
                vertices.append(np.asarray(vertex, dtype=float))
            indices.append(lookup[key])
        faces.append(indices)
    edges = set()
    for face in faces:
        for first, second in zip(face, face[1:] + face[:1]):
            edges.add(tuple(sorted((first, second))))
    return {"vertices": np.asarray(vertices), "faces": faces,
            "edges": sorted(edges)}
