"""Native symmetry reduction of uniform reciprocal-space meshes.

The crystal symmetry search is delegated to spglib.  Grid generation, mesh
stabiliser detection, orbit construction, representative selection, weights,
and Quantum ESPRESSO export are implemented here.

All grid equivalence operations use integer addresses.  Floating-point
coordinates are created only for presentation and export, so boundary points
and half-grid shifts are never matched with a numerical tolerance.
"""

from itertools import product
from math import lcm, prod

import numpy as np
import spglib

from qeplotter.kpath.geometry import first_brillouin_zone, reciprocal_rows
from qeplotter.structure.symmetry import crystal_system


DEFAULT_MAX_GRID_POINTS = 500_000


def _dataset_value(dataset, key):
    return dataset[key] if isinstance(dataset, dict) else getattr(dataset, key)


def _validate_triplet(values, name, allowed=None):
    try:
        raw = tuple(values)
    except TypeError as error:
        raise ValueError(f"{name} must contain exactly three integers") from error
    if len(raw) != 3:
        raise ValueError(f"{name} must contain exactly three integers")
    parsed = []
    for value in raw:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} values must be integers, not booleans")
        numeric = float(value)
        if not np.isfinite(numeric):
            raise ValueError(f"{name} values must be finite integers")
        integer = int(numeric)
        if numeric != integer:
            raise ValueError(f"{name} values must be integers")
        parsed.append(integer)
    parsed = tuple(parsed)
    if allowed is not None and any(value not in allowed for value in parsed):
        choices = " or ".join(str(value) for value in sorted(allowed))
        raise ValueError(f"{name} values must be {choices}")
    return parsed


def _validate_inputs(mesh, shift, symprec, max_grid_points):
    mesh = _validate_triplet(mesh, "k-point mesh")
    shift = _validate_triplet(shift, "k-point shift", allowed={0, 1})
    if any(value < 1 for value in mesh):
        raise ValueError("k-point mesh values must be positive")
    symprec = float(symprec)
    if not np.isfinite(symprec) or symprec <= 0:
        raise ValueError("symmetry tolerance must be a positive finite number")
    total = int(prod(mesh))
    if total > int(max_grid_points):
        raise ValueError(
            f"The requested grid contains {total:,} points; the safe limit is "
            f"{int(max_grid_points):,}. Reduce the mesh dimensions."
        )
    return mesh, shift, symprec, total


def _reciprocal_rotation(direct_rotation):
    """Return the integer reciprocal action R=(W^-1)^T."""
    direct = np.asarray(direct_rotation, dtype=int)
    reciprocal = np.rint(np.linalg.inv(direct).T).astype(int)
    if not np.array_equal(direct.T @ reciprocal, np.eye(3, dtype=int)):
        raise ValueError("spglib returned a non-unimodular symmetry rotation")
    return reciprocal


def _integer_grid_frame(mesh, shift):
    """Put all three grid axes on one exact cyclic integer denominator."""
    denominator = lcm(*(2 * value for value in mesh))
    steps = np.asarray([denominator // value for value in mesh], dtype=np.int64)
    offsets = np.asarray(
        [flag * denominator // (2 * value)
         for value, flag in zip(mesh, shift)],
        dtype=np.int64,
    )
    return int(denominator), steps, offsets


def _preserves_grid(rotation, denominator, steps, offsets):
    """Test whether a reciprocal rotation stabilises this shifted mesh."""
    rotation = np.asarray(rotation, dtype=np.int64)
    moved_origin = rotation @ offsets - offsets
    if any(int(moved_origin[axis]) % int(steps[axis])
           for axis in range(3)):
        return False
    for target in range(3):
        for source in range(3):
            delta = int(rotation[target, source]) * int(steps[source])
            if delta % int(steps[target]):
                return False
    # The congruence checks above imply preservation on the cyclic cell.  Keep
    # denominator in the signature to document the exact ambient group.
    return denominator > 0


def _unique_reciprocal_rotations(direct_rotations):
    rotations = {}
    for direct in direct_rotations:
        reciprocal = _reciprocal_rotation(direct)
        rotations[tuple(int(value) for value in reciprocal.flat)] = reciprocal
    return list(rotations.values())


def _grid_indices(mesh):
    axes = np.meshgrid(*(np.arange(value, dtype=np.int64)
                         for value in mesh), indexing="ij")
    return np.stack([axis.ravel() for axis in axes], axis=1)


def _transform_grid(indices, rotation, mesh, denominator, steps, offsets):
    scaled = offsets + indices * steps
    transformed = np.mod(scaled @ np.asarray(rotation, dtype=np.int64).T,
                         denominator)
    delta = np.mod(transformed - offsets, denominator)
    remainder = np.mod(delta, steps)
    if np.any(remainder):
        raise RuntimeError("An accepted symmetry operation left the selected grid")
    return np.mod(delta // steps, np.asarray(mesh, dtype=np.int64))


def _center_fractional(indices, mesh, shift):
    values = (
        np.asarray(indices, dtype=float)
        + 0.5 * np.asarray(shift, dtype=float)
    ) / np.asarray(mesh, dtype=float)
    # A deterministic half-open cell makes +0.5 and -0.5 duplicates impossible.
    return np.mod(values + 0.5, 1.0) - 0.5


def _spglib_cell(atoms):
    lattice = np.asarray(atoms.cell[:], dtype=float)
    if lattice.shape != (3, 3) or not np.all(np.isfinite(lattice)):
        raise ValueError("The structure must contain a finite 3×3 unit cell")
    if abs(float(np.linalg.det(lattice))) <= 1e-12:
        raise ValueError("A non-zero periodic unit cell is required")
    positions = np.mod(
        np.asarray(atoms.get_scaled_positions(wrap=False), dtype=float), 1.0)
    numbers = np.asarray(atoms.get_atomic_numbers(), dtype=int)
    if len(positions) == 0:
        raise ValueError("The structure does not contain any atoms")
    return lattice, positions, numbers


def _validate_result(result):
    total = result["total_grid_points"]
    multiplicities = np.asarray(
        [point["multiplicity"] for point in result["points"]], dtype=int)
    if int(multiplicities.sum()) != total:
        raise RuntimeError("Internal error: k-point multiplicities do not cover the full grid")
    weights = np.asarray(
        [point["normalized_weight"] for point in result["points"]], dtype=float)
    if not np.isclose(float(weights.sum()), 1.0, atol=5e-14):
        raise RuntimeError("Internal error: normalized k-point weights do not sum to one")
    mapping = np.asarray(result["ir_index_by_grid_index"], dtype=int)
    if mapping.shape != (total,) or mapping.min() != 0:
        raise RuntimeError("Internal error: invalid full-grid to IBZ mapping")
    if mapping.max(initial=-1) != result["irreducible_count"] - 1:
        raise RuntimeError("Internal error: incomplete irreducible index range")


def irreducible_kmesh(
    atoms,
    mesh,
    shift=(0, 0, 0),
    symprec=1e-5,
    time_reversal=True,
    max_grid_points=DEFAULT_MAX_GRID_POINTS,
):
    """Reduce a uniform k-grid using crystal symmetry operations.

    spglib is used only to identify the structure's spatial symmetry.  This
    function performs the reciprocal rotation conversion, shifted-grid
    compatibility test, orbit reduction, weights, and representative mapping.

    Parameters
    ----------
    atoms
        ASE ``Atoms`` with a finite unit cell.
    mesh
        Three positive mesh dimensions ``(nk1, nk2, nk3)``.
    shift
        QE half-step flags; every value is 0 or 1.
    symprec
        spglib structure-symmetry distance tolerance in Å.
    time_reversal
        Add ``k -> -k`` equivalence.  Use only when the calculation preserves
        time-reversal symmetry.
    max_grid_points
        Guard against accidental grids that are too large for an interactive
        application.
    """
    mesh, shift, symprec, total = _validate_inputs(
        mesh, shift, symprec, max_grid_points)
    lattice, positions, numbers = _spglib_cell(atoms)
    structure = (lattice, positions, numbers)

    dataset = spglib.get_symmetry_dataset(structure, symprec=symprec)
    symmetry = spglib.get_symmetry(structure, symprec=symprec)
    if dataset is None or symmetry is None:
        raise ValueError("spglib could not determine the crystal symmetry")

    spatial_rotations = _unique_reciprocal_rotations(symmetry["rotations"])
    denominator, steps, offsets = _integer_grid_frame(mesh, shift)
    compatible_spatial = [
        rotation for rotation in spatial_rotations
        if _preserves_grid(rotation, denominator, steps, offsets)
    ]
    if not compatible_spatial:
        compatible_spatial = [np.eye(3, dtype=int)]

    transforms = {
        tuple(int(value) for value in rotation.flat): rotation
        for rotation in compatible_spatial
    }
    if bool(time_reversal):
        for rotation in compatible_spatial:
            reversed_rotation = -rotation
            if _preserves_grid(
                    reversed_rotation, denominator, steps, offsets):
                transforms[tuple(int(value)
                                 for value in reversed_rotation.flat)] = reversed_rotation
    transforms = list(transforms.values())

    indices = _grid_indices(mesh)
    representative_by_flat = np.arange(total, dtype=np.int64)
    for rotation in transforms:
        transformed = _transform_grid(
            indices, rotation, mesh, denominator, steps, offsets)
        transformed_flat = np.ravel_multi_index(
            tuple(transformed[:, axis] for axis in range(3)), mesh)
        representative_by_flat = np.minimum(
            representative_by_flat, transformed_flat)

    representatives, ir_index, multiplicities = np.unique(
        representative_by_flat, return_inverse=True, return_counts=True)
    representative_indices = indices[representatives]
    representative_fractional = _center_fractional(
        representative_indices, mesh, shift)

    points = []
    for serial, (flat_index, address, fractional, multiplicity) in enumerate(
            zip(representatives, representative_indices,
                representative_fractional, multiplicities), start=1):
        points.append({
            "index": serial,
            "grid_index": int(flat_index),
            "address": tuple(int(value) for value in address),
            "frac": np.asarray(fractional, dtype=float),
            "multiplicity": int(multiplicity),
            # Compatibility alias for early QEPlotter 2.0 callers.
            "weight": int(multiplicity),
            "normalized_weight": float(multiplicity) / total,
        })

    reciprocal = reciprocal_rows(lattice)
    number = int(_dataset_value(dataset, "number"))
    result = {
        "mesh": mesh,
        "shift": shift,
        "symprec": symprec,
        "time_reversal": bool(time_reversal),
        "spacegroup_number": number,
        "spacegroup_international": str(
            _dataset_value(dataset, "international")),
        "crystal_system": crystal_system(number),
        "total_grid_points": total,
        "irreducible_count": len(points),
        "points": points,
        "grid_indices": indices,
        "representative_by_grid_index": representative_by_flat,
        "ir_index_by_grid_index": ir_index,
        "reciprocal_input": reciprocal,
        "bz": first_brillouin_zone(reciprocal),
        "coordinate_basis":
            "fractional coordinates of the uploaded cell reciprocal basis",
        "reciprocal_convention":
            "row vectors, B=2π(A⁻¹)ᵀ, k_cart=k_frac·B",
        "engine": "QEPlotter native integer-orbit reducer",
        "symmetry_source": "spglib spatial operations",
        "detected_spacegroup_operations": int(len(symmetry["rotations"])),
        "detected_unique_point_rotations": len(spatial_rotations),
        "compatible_spatial_rotations": len(compatible_spatial),
        "equivalence_transforms": len(transforms),
        "full_crystal_symmetry_preserved":
            len(compatible_spatial) == len(spatial_rotations),
        "dropped_spatial_rotations":
            len(spatial_rotations) - len(compatible_spatial),
    }
    _validate_result(result)
    return result


def orbit_members(result, irreducible_index):
    """Return all full-grid members represented by one IBZ point."""
    index = int(irreducible_index)
    if not 0 <= index < result["irreducible_count"]:
        raise IndexError("irreducible k-point index is out of range")
    mask = np.asarray(result["ir_index_by_grid_index"]) == index
    addresses = np.asarray(result["grid_indices"])[mask]
    fractional = _center_fractional(
        addresses, result["mesh"], result["shift"])
    return [
        {
            "grid_index": int(flat),
            "address": tuple(int(value) for value in address),
            "frac": frac,
        }
        for flat, address, frac in zip(
            np.flatnonzero(mask), addresses, fractional)
    ]


def full_grid_points(result):
    """Return the complete uniform grid with its irreducible mapping."""
    fractional = _center_fractional(
        result["grid_indices"], result["mesh"], result["shift"])
    return [
        {
            "grid_index": index,
            "address": tuple(int(value) for value in address),
            "frac": frac,
            "irreducible_index": int(result["ir_index_by_grid_index"][index]) + 1,
        }
        for index, (address, frac) in enumerate(
            zip(result["grid_indices"], fractional))
    ]


def format_qe_automatic(result):
    """Return the equivalent QE ``K_POINTS automatic`` card."""
    mesh = " ".join(str(value) for value in result["mesh"])
    shift = " ".join(str(value) for value in result["shift"])
    rows = []
    if not result["time_reversal"]:
        rows.extend([
            "! QEPlotter: k <-> -k equivalence was disabled.",
            "! Set noinv=.true. in &SYSTEM to reproduce this choice in pw.x.",
        ])
    rows.extend(["K_POINTS automatic", f" {mesh} {shift}"])
    return "\n".join(rows) + "\n"


def format_qe_ir_kpoints(result, weight_mode="normalized"):
    """Return explicit weighted IBZ points as a QE ``K_POINTS crystal`` card."""
    if weight_mode not in {"normalized", "multiplicity"}:
        raise ValueError(
            "weight_mode must be 'normalized' or 'multiplicity'")
    rows = []
    if not result["time_reversal"]:
        rows.extend([
            "! QEPlotter: k <-> -k equivalence was disabled.",
            "! Set noinv=.true. in &SYSTEM to reproduce this choice in pw.x.",
        ])
    rows.extend(["K_POINTS crystal", str(result["irreducible_count"])])
    for point in result["points"]:
        x, y, z = point["frac"]
        weight = (
            point["normalized_weight"]
            if weight_mode == "normalized"
            else point["multiplicity"]
        )
        rows.append(
            f" {x: .10f} {y: .10f} {z: .10f} {weight: .12g}"
            f" ! multiplicity={point['multiplicity']}"
        )
    return "\n".join(rows) + "\n"
