"""Native conventional band-path orchestration.

The geometry engine is separate.  Paths follow Setyawan--Curtarolo/AFLOW and
are transformed explicitly from the AFLOW canonical primitive reciprocal basis
to the uploaded cell reciprocal basis.  See ``recipes_sc.py`` for formulae.
"""
import re

import numpy as np
import spglib
from ase.cell import Cell
from ase.lattice import identify_lattice

from qeplotter.kpath.geometry import first_brillouin_zone, reciprocal_rows
from qeplotter.kpath.recipes_sc import G as GAMMA, get_sc_recipe
from qeplotter.structure.symmetry import crystal_system


_PEARSON = {
    "CUB":"cP", "FCC":"cF", "BCC":"cI", "TET":"tP", "BCT":"tI",
    "ORC":"oP", "ORCF":"oF", "ORCI":"oI", "ORCC":"oC",
    "HEX":"hP", "RHL":"hR", "MCL":"mP", "MCLC":"mC", "TRI":"aP",
}


def _get(dataset, key):
    return dataset[key] if isinstance(dataset, dict) else getattr(dataset, key)


def _interpolate_path(points_cart, path, inverse_input, reference_distance):
    explicit = []
    for segment_index, (start, end) in enumerate(path):
        first, last = points_cart[start], points_cart[end]
        length = float(np.linalg.norm(last - first))
        count = max(2, int(np.ceil(length / max(reference_distance, 1e-4))) + 1)
        for index, fraction in enumerate(np.linspace(0.0, 1.0, count)):
            cart = first + fraction * (last - first)
            explicit.append({"segment": segment_index, "start": start, "end": end,
                             "label": start if index == 0 else end if index == count - 1 else "",
                             "frac": cart @ inverse_input, "cart": cart})
    return explicit


def primary_path(path):
    if not path:
        return []
    primary = [path[0]]
    for segment in path[1:]:
        if segment[0] != primary[-1][1]:
            break
        primary.append(segment)
    return primary


def parse_path_expression(expression, available_labels):
    aliases = {str(label).upper(): label for label in available_labels}
    aliases.update({"Γ": GAMMA, "G": GAMMA, "GAMMA": GAMMA})
    segments = []
    for raw_branch in expression.split("|"):
        tokens = [token.strip() for token in re.split(r"\s*(?:-|–|—|→)\s*", raw_branch)
                  if token.strip()]
        if len(tokens) < 2:
            raise ValueError(f"Each branch needs at least two points: '{raw_branch.strip()}'")
        labels = []
        for token in tokens:
            key = token.upper()
            if key not in aliases:
                choices = ", ".join("Γ" if label == GAMMA else label
                                    for label in available_labels)
                raise ValueError(f"Unknown point '{token}'. Available points: {choices}")
            labels.append(aliases[key])
        segments.extend(zip(labels, labels[1:]))
    if not segments:
        raise ValueError("The custom path is empty")
    return list(segments)


def with_path(result, path, reference_distance=None):
    available = result["point_coords_cart"]
    clean_path = [(str(start), str(end)) for start, end in path]
    for start, end in clean_path:
        if start not in available or end not in available:
            raise ValueError(f"Path segment {start}-{end} uses an unavailable point")
        if start == end:
            raise ValueError(f"Path segment {start}-{end} has zero length")
    updated = dict(result)
    updated["path"] = clean_path
    spacing = result["reference_distance"] if reference_distance is None else float(reference_distance)
    inverse_input = np.linalg.inv(result["reciprocal_input"])
    updated["explicit"] = _interpolate_path(available, clean_path, inverse_input, spacing)
    updated["reference_distance"] = spacing
    return updated


def recommend_kpath(atoms, symprec=1e-3, reference_distance=0.05):
    """Return a conventional SC path, explicit basis metadata and BZ geometry."""
    input_cell = np.asarray(atoms.cell[:], dtype=float)
    structure = (input_cell, np.asarray(atoms.get_scaled_positions(), dtype=float),
                 np.asarray(atoms.get_atomic_numbers(), dtype=int))
    dataset = spglib.get_symmetry_dataset(structure, symprec=symprec)
    if dataset is None:
        raise ValueError("spglib could not determine the crystal symmetry")

    primitive = spglib.standardize_cell(
        structure, to_primitive=True, no_idealize=True, symprec=symprec)
    if primitive is None:
        raise ValueError("spglib could not extract a primitive cell")
    primitive_cell = np.asarray(primitive[0], dtype=float)
    metric_eps = max(2e-4, 2.0 * symprec / max(np.linalg.norm(primitive_cell, axis=1)))
    try:
        lattice, canonical_to_primitive = identify_lattice(
            Cell(primitive_cell), eps=metric_eps)
    except Exception as error:
        raise ValueError(f"AFLOW primitive-cell canonicalisation failed: {error}") from error

    recipe = get_sc_recipe(lattice)
    canonical_points = recipe["points"]
    # ASE's identify_lattice guarantees: primitive ≈ op.T @ canonical_cell,
    # and reciprocal fractional row coordinates transform as k_prim=k_AFLOW@op.
    points_primitive = {label: coords @ canonical_to_primitive
                        for label, coords in canonical_points.items()}
    reciprocal_primitive = reciprocal_rows(primitive_cell)
    reciprocal_input = reciprocal_rows(input_cell)
    points_cart = {label: coords @ reciprocal_primitive
                   for label, coords in points_primitive.items()}
    inverse_input = np.linalg.inv(reciprocal_input)
    point_coords = {label: cart @ inverse_input for label, cart in points_cart.items()}

    symmetry = spglib.get_symmetry(structure, symprec=symprec)
    has_inversion = bool(symmetry is not None and any(
        np.array_equal(rotation, -np.eye(3, dtype=int))
        for rotation in symmetry["rotations"]))
    number = int(_get(dataset, "number"))
    bravais = _PEARSON[lattice.name]
    extended = bravais + recipe["variant"][-1] if recipe["variant"][-1].isdigit() else bravais + "1"
    explicit = _interpolate_path(points_cart, recipe["path"], inverse_input,
                                 reference_distance)

    mapped_canonical = canonical_to_primitive.T @ np.asarray(lattice.tocell())
    # A global rigid rotation is irrelevant for fractional coordinates; compare
    # basis metrics rather than Cartesian components.
    canonical_metric = mapped_canonical @ mapped_canonical.T
    primitive_metric = primitive_cell @ primitive_cell.T
    mapping_error = float(np.max(np.abs(canonical_metric - primitive_metric)))
    if mapping_error > max(5e-5, 4 * symprec):
        raise ValueError(f"Primitive-basis metric residual is too large ({mapping_error:.3e} Å²)")

    return {
        "spacegroup_number": number,
        "spacegroup_international": str(_get(dataset, "international")),
        "crystal_system": crystal_system(number), "bravais_lattice": bravais,
        "bravais_lattice_extended": extended, "recipe_variant": recipe["variant"],
        "convention": "Setyawan–Curtarolo (AFLOW, 2010)",
        "method": "conventional-recipe", "warning": None,
        "recipe_parameters": recipe["parameters"],
        "point_coords_aflow": canonical_points,
        "point_coords_primitive": points_primitive,
        "point_coords": point_coords, "point_coords_cart": points_cart,
        "path": recipe["path"], "explicit": explicit,
        "reference_distance": float(reference_distance),
        "reciprocal_input": reciprocal_input,
        "reciprocal_primitive": reciprocal_primitive,
        "primitive_lattice": primitive_cell,
        "canonical_primitive_lattice": np.asarray(lattice.tocell()),
        "canonical_to_primitive": np.asarray(canonical_to_primitive),
        "basis_mapping_error": mapping_error,
        "coordinate_basis": "fractional reciprocal coordinates of uploaded input cell",
        "reciprocal_convention": "row vectors, B=2π(A⁻¹)ᵀ, k_cart=k_frac·B",
        "bz": first_brillouin_zone(reciprocal_primitive),
        "has_inversion": has_inversion,
    }


def format_qe_kpoints(result, explicit=False):
    if result.get("method") != "conventional-recipe":
        raise ValueError("QE band-path export requires a conventional recipe")
    if explicit:
        rows = ["K_POINTS crystal", str(len(result["explicit"]))]
        for item in result["explicit"]:
            x, y, z = item["frac"]
            suffix = f" ! {item['label']}" if item["label"] else ""
            rows.append(f" {x: .10f} {y: .10f} {z: .10f} 1.0{suffix}")
        return "\n".join(rows) + "\n"
    counts = {}
    for item in result["explicit"]:
        counts[item["segment"]] = counts.get(item["segment"], 0) + 1
    vertices = []
    for index, (start, end) in enumerate(result["path"]):
        vertices.extend([(start, result["point_coords"][start], counts[index]),
                         (end, result["point_coords"][end], 1)])
    rows = ["K_POINTS crystal_b", str(len(vertices))]
    for label, coords, count in vertices:
        x, y, z = coords
        rows.append(f" {x: .10f} {y: .10f} {z: .10f} {count:d} ! {label}")
    return "\n".join(rows) + "\n"
