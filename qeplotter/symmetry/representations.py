"""Γ-point orbital and displacement representations.

The implementation deliberately separates three ideas that are often mixed in
hand calculations:

* spglib finds the primitive crystal and its symmetry operations;
* a finite-group engine obtains the irreducible characters directly from the
  point group (there is no hard-coded character table or Seek-path call);
* the selected atomic/orbital basis is transformed and projected into symmetry
  adapted linear combinations (SALCs).

Only structure-derived Γ representations are handled here.  Classifying actual
electronic eigenstates still requires wavefunctions from the electronic-
structure calculation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import spglib
from ase import Atoms
from ase.data import chemical_symbols


class SymmetryAnalysisError(ValueError):
    """Raised when a structure cannot support a reliable symmetry analysis."""


_BASIS_COMPONENTS = {
    "s": ("s",),
    "p": ("pₓ", "pᵧ", "p_z"),
    "d": ("d_z²", "d_x²-y²", "d_xy", "d_xz", "d_yz"),
    "displacement": ("uₓ", "uᵧ", "u_z"),
}

_BASIS_TITLES = {
    "s": "s orbitals",
    "p": "p orbitals (pₓ, pᵧ, p_z)",
    "d": "d orbitals (five real components)",
    "displacement": "atomic displacements (uₓ, uᵧ, u_z)",
}


def _dataset_value(dataset, name, default=None):
    if dataset is None:
        return default
    if isinstance(dataset, dict):
        return dataset.get(name, default)
    return getattr(dataset, name, default)


def _clean_number(value: complex, tolerance: float = 1e-7):
    value = complex(value)
    real = 0.0 if abs(value.real) < tolerance else value.real
    imag = 0.0 if abs(value.imag) < tolerance else value.imag
    if abs(real - round(real)) < tolerance:
        real = float(round(real))
    if abs(imag - round(imag)) < tolerance:
        imag = float(round(imag))
    return complex(real, imag)


def format_character(value: complex) -> str:
    """Return a compact, human-readable character or SALC coefficient."""
    value = _clean_number(value)
    if abs(value.imag) < 1e-8:
        number = value.real
        if abs(number - round(number)) < 1e-8:
            return str(int(round(number)))
        return f"{number:.4f}".rstrip("0").rstrip(".")
    if abs(value.real) < 1e-8:
        return f"{value.imag:.4f}i".replace("1.0000i", "i").replace("-1.0000i", "-i")
    sign = "+" if value.imag >= 0 else "-"
    return f"{value.real:.4f}{sign}{abs(value.imag):.4f}i"


def _orthogonal_cartesian_rotation(lattice: np.ndarray, rotation: np.ndarray):
    # spglib uses column fractional coordinates: x' = R x + t.
    cart = lattice.T @ rotation @ np.linalg.inv(lattice.T)
    u, _, vt = np.linalg.svd(cart)
    return u @ vt


def _matrix_order(rotation: np.ndarray, maximum: int = 24) -> int:
    current = np.eye(3, dtype=int)
    for order in range(1, maximum + 1):
        current = current @ rotation
        if np.array_equal(current, np.eye(3, dtype=int)):
            return order
    return 0


def _operation_kind(rotation: np.ndarray, cart_rotation: np.ndarray) -> str:
    if np.array_equal(rotation, np.eye(3, dtype=int)):
        return "E"
    if np.array_equal(rotation, -np.eye(3, dtype=int)):
        return "i"
    determinant = int(round(np.linalg.det(rotation)))
    order = _matrix_order(rotation)
    if determinant > 0:
        return f"C{order}" if order else "proper rotation"
    if order == 2 and abs(np.trace(cart_rotation) - 1.0) < 1e-6:
        return "σ"
    return f"S{order}" if order else "improper rotation"


def _class_labels(
    classes: Sequence[Sequence[int]],
    rotations: Sequence[np.ndarray],
    cart_rotations: Sequence[np.ndarray],
) -> List[str]:
    raw = []
    for members in classes:
        kind = _operation_kind(rotations[members[0]], cart_rotations[members[0]])
        prefix = f"{len(members)} " if len(members) > 1 else ""
        raw.append(f"{prefix}{kind}")
    totals: Dict[str, int] = {}
    for label in raw:
        totals[label] = totals.get(label, 0) + 1
    seen: Dict[str, int] = {}
    labels = []
    for label in raw:
        seen[label] = seen.get(label, 0) + 1
        labels.append(
            f"{label} ({seen[label]})" if totals[label] > 1 else label
        )
    return labels


def _multiplication_table(rotations: Sequence[np.ndarray]):
    lookup = {tuple(np.asarray(rotation, dtype=int).ravel()): i
              for i, rotation in enumerate(rotations)}
    size = len(rotations)
    table = np.empty((size, size), dtype=int)
    for i, first in enumerate(rotations):
        for j, second in enumerate(rotations):
            key = tuple((first @ second).ravel())
            if key not in lookup:
                raise SymmetryAnalysisError(
                    "The detected point operations do not form a closed group."
                )
            table[i, j] = lookup[key]
    return table


def _conjugacy_classes(table: np.ndarray):
    size = len(table)
    identity = next(
        i for i in range(size)
        if np.array_equal(table[i], np.arange(size))
        and np.array_equal(table[:, i], np.arange(size))
    )
    inverse = np.empty(size, dtype=int)
    for element in range(size):
        candidates = np.where(
            (table[element] == identity) & (table[:, element] == identity)
        )[0]
        inverse[element] = int(candidates[0])

    unseen = set(range(size))
    classes = []
    while unseen:
        element = identity if identity in unseen else min(unseen)
        members = {
            int(table[table[h, element], inverse[h]]) for h in range(size)
        }
        classes.append(sorted(members))
        unseen.difference_update(members)
    return identity, inverse, classes


def _cluster_eigenvalues(values: np.ndarray, tolerance: float):
    groups: List[List[int]] = []
    for index, value in enumerate(values):
        if not groups or abs(value - values[groups[-1][-1]]) > tolerance:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def _complex_irreducible_characters(
    rotations: Sequence[np.ndarray],
    classes: Sequence[Sequence[int]],
    table: np.ndarray,
):
    """Calculate the complex irreducible character table from the regular rep."""
    size = len(rotations)
    class_sums = []
    for members in classes:
        matrix = np.zeros((size, size), dtype=complex)
        for element in members:
            for source in range(size):
                matrix[table[element, source], source] += 1.0
        class_sums.append(matrix)

    for attempt in range(12):
        rng = np.random.default_rng(872341 + attempt)
        central = np.zeros((size, size), dtype=complex)
        for class_sum in class_sums:
            coefficient = rng.normal() + 1j * rng.normal()
            central += (
                coefficient * class_sum
                + coefficient.conjugate() * class_sum.conjugate().T
            )
        central = (central + central.conjugate().T) / 2.0
        eigenvalues, eigenvectors = np.linalg.eigh(central)
        tolerance = max(1e-8, np.max(np.abs(eigenvalues)) * 1e-8)
        clusters = _cluster_eigenvalues(eigenvalues, tolerance)
        dimensions = [int(round(np.sqrt(len(indices)))) for indices in clusters]
        valid = (
            len(clusters) == len(classes)
            and all(dimension * dimension == len(indices)
                    for dimension, indices in zip(dimensions, clusters))
            and sum(dimension * dimension for dimension in dimensions) == size
        )
        if not valid:
            continue

        irreps = []
        for dimension, indices in zip(dimensions, clusters):
            vectors = eigenvectors[:, indices]
            characters = []
            for members, class_sum in zip(classes, class_sums):
                eigenvalue = (
                    np.trace(vectors.conjugate().T @ class_sum @ vectors)
                    / len(indices)
                )
                characters.append(
                    _clean_number(eigenvalue * dimension / len(members))
                )
            irreps.append({
                "dimension": dimension,
                "characters": np.asarray(characters, dtype=complex),
            })

        # The character rows must be orthonormal under the class inner product.
        gram = np.empty((len(irreps), len(irreps)), dtype=complex)
        weights = np.asarray([len(members) for members in classes])
        for i, left in enumerate(irreps):
            for j, right in enumerate(irreps):
                gram[i, j] = np.sum(
                    weights * left["characters"]
                    * right["characters"].conjugate()
                ) / size
        if np.allclose(gram, np.eye(len(irreps)), atol=2e-5):
            return irreps

    raise SymmetryAnalysisError(
        "Could not separate the numerical irreducible character table."
    )


def _character_sort_key(irrep):
    characters = irrep["characters"]
    trivial = int(not np.allclose(characters, 1.0, atol=1e-6))
    flattened = tuple(
        value
        for character in characters
        for value in (
            round(float(character.real), 7),
            round(float(character.imag), 7),
        )
    )
    return trivial, irrep["dimension"], flattened


def _combine_real_irreps(complex_irreps):
    """Combine complex-conjugate pairs into real crystallographic irreps."""
    used = set()
    real_irreps = []
    for index, irrep in enumerate(complex_irreps):
        if index in used:
            continue
        characters = irrep["characters"]
        if np.allclose(characters.imag, 0.0, atol=1e-7):
            used.add(index)
            real_irreps.append({
                "dimension": irrep["dimension"],
                "characters": characters.real.astype(complex),
                "members": (index,),
            })
            continue
        partner = next(
            (
                j for j, candidate in enumerate(complex_irreps)
                if j != index and j not in used
                and candidate["dimension"] == irrep["dimension"]
                and np.allclose(
                    candidate["characters"], characters.conjugate(), atol=1e-6
                )
            ),
            None,
        )
        if partner is None:
            raise SymmetryAnalysisError(
                "A complex irrep did not have its conjugate partner."
            )
        used.update((index, partner))
        real_irreps.append({
            "dimension": 2 * irrep["dimension"],
            "characters": (
                characters + complex_irreps[partner]["characters"]
            ),
            "members": (index, partner),
        })
    return sorted(real_irreps, key=_character_sort_key)


def _matches_characters(left, right, tolerance=1e-5):
    return np.allclose(
        np.asarray(left, dtype=complex),
        np.asarray(right, dtype=complex),
        atol=tolerance,
    )


def _assign_irrep_labels(
    pointgroup: str,
    irreps: List[dict],
    classes: Sequence[Sequence[int]],
    cart_rotations: Sequence[np.ndarray],
):
    """Use exact conventional labels where they can be assigned unambiguously."""
    polar = np.asarray([
        np.mean([np.trace(cart_rotations[index]) for index in members])
        for members in classes
    ])
    axial = np.asarray([
        np.mean([
            np.linalg.det(cart_rotations[index])
            * np.trace(cart_rotations[index])
            for index in members
        ])
        for members in classes
    ])
    inversion_class = next(
        (
            class_index for class_index, members in enumerate(classes)
            if np.allclose(cart_rotations[members[0]], -np.eye(3), atol=1e-6)
        ),
        None,
    )

    exact_cubic = pointgroup in {"-43m", "432", "m-3m"}
    counters: Dict[Tuple[str, str], int] = {}
    for serial, irrep in enumerate(irreps, start=1):
        dimension = irrep["dimension"]
        characters = irrep["characters"]
        parity = ""
        if inversion_class is not None:
            ratio = characters[inversion_class].real / dimension
            parity = "g" if ratio > 0 else "u"

        label = None
        source = "internal"
        if exact_cubic:
            if dimension == 1:
                if pointgroup == "-43m":
                    # In Td, A1 and A2 agree on every proper operation and
                    # differ only on S4/σd.
                    index = (
                        "1" if np.allclose(characters, 1.0, atol=1e-6)
                        else "2"
                    )
                else:
                    proper_values = [
                        characters[i] for i, members in enumerate(classes)
                        if np.linalg.det(cart_rotations[members[0]]) > 0
                    ]
                    index = (
                        "1" if np.allclose(proper_values, 1.0, atol=1e-6)
                        else "2"
                    )
                label = f"A{index}{parity}"
            elif dimension == 2:
                label = f"E{parity}"
            elif dimension == 3:
                # Td distinguishes polar T2 from axial T1. O/Oh use T1 for
                # the vector restriction to proper rotations.
                if pointgroup == "-43m":
                    index = "2" if _matches_characters(characters, polar) else "1"
                else:
                    proper = [
                        i for i, members in enumerate(classes)
                        if np.linalg.det(cart_rotations[members[0]]) > 0
                    ]
                    index = (
                        "1" if _matches_characters(
                            characters[proper], polar[proper]
                        ) else "2"
                    )
                label = f"T{index}{parity}"
            if label:
                source = "conventional"

        if label is None:
            family = {1: "A", 2: "E", 3: "T", 4: "G", 5: "H"}.get(
                dimension, f"D{dimension}"
            )
            key = (family, parity)
            counters[key] = counters.get(key, 0) + 1
            family_index = counters[key]
            label = f"{family}{family_index if counters[key] > 1 else ''}{parity}"

        irrep.update(
            label=label,
            gamma_label=f"Γ{serial}",
            label_source=source,
            polar_vector=_matches_characters(characters, polar),
            axial_vector=_matches_characters(characters, axial),
        )


def _d_orbital_matrix(rotation: np.ndarray):
    root2 = np.sqrt(2.0)
    root6 = np.sqrt(6.0)
    basis = np.asarray([
        np.diag([-1.0, -1.0, 2.0]) / root6,
        np.diag([1.0, -1.0, 0.0]) / root2,
        [[0, 1 / root2, 0], [1 / root2, 0, 0], [0, 0, 0]],
        [[0, 0, 1 / root2], [0, 0, 0], [1 / root2, 0, 0]],
        [[0, 0, 0], [0, 0, 1 / root2], [0, 1 / root2, 0]],
    ])
    matrix = np.empty((5, 5), dtype=float)
    for source, tensor in enumerate(basis):
        transformed = rotation @ tensor @ rotation.T
        for target, target_tensor in enumerate(basis):
            matrix[target, source] = np.sum(target_tensor * transformed)
    return matrix


def _orbital_matrix(basis: str, rotation: np.ndarray):
    if basis == "s":
        return np.ones((1, 1))
    if basis in {"p", "displacement"}:
        return rotation
    if basis == "d":
        return _d_orbital_matrix(rotation)
    raise SymmetryAnalysisError(f"Unsupported basis: {basis}")


def _primitive_structure(atoms: Atoms, symprec: float):
    if atoms.cell.volume <= 1e-8:
        raise SymmetryAnalysisError(
            "A non-zero periodic unit cell is required for crystal symmetry."
        )
    cell = (
        np.asarray(atoms.cell[:], dtype=float),
        np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float),
        np.asarray(atoms.get_atomic_numbers(), dtype=int),
    )
    primitive = spglib.standardize_cell(
        cell, to_primitive=True, no_idealize=False, symprec=symprec
    )
    if primitive is None:
        raise SymmetryAnalysisError(
            "spglib could not construct a primitive standardized cell."
        )
    lattice, positions, numbers = primitive
    result = Atoms(
        numbers=np.asarray(numbers, dtype=int),
        scaled_positions=np.asarray(positions, dtype=float),
        cell=np.asarray(lattice, dtype=float),
        pbc=True,
    )
    dataset = spglib.get_symmetry_dataset(
        (
            np.asarray(result.cell[:]),
            np.asarray(result.get_scaled_positions(wrap=True)),
            np.asarray(result.get_atomic_numbers()),
        ),
        symprec=symprec,
    )
    if dataset is None:
        raise SymmetryAnalysisError(
            "spglib could not determine symmetry for the primitive cell."
        )
    return result, dataset


def _unique_operations(dataset):
    rotations = np.asarray(_dataset_value(dataset, "rotations"), dtype=int)
    translations = np.asarray(
        _dataset_value(dataset, "translations"), dtype=float
    )
    unique_rotations = []
    unique_translations = []
    seen = set()
    for rotation, translation in zip(rotations, translations):
        key = tuple(rotation.ravel())
        if key in seen:
            continue
        seen.add(key)
        unique_rotations.append(rotation)
        unique_translations.append(translation)
    identity = next(
        i for i, rotation in enumerate(unique_rotations)
        if np.array_equal(rotation, np.eye(3, dtype=int))
    )
    order = [identity] + sorted(
        (i for i in range(len(unique_rotations)) if i != identity),
        key=lambda i: tuple(unique_rotations[i].ravel()),
    )
    return (
        np.asarray([unique_rotations[i] for i in order], dtype=int),
        np.asarray([unique_translations[i] for i in order], dtype=float),
    )


def _atom_permutations(
    atoms: Atoms,
    rotations: Sequence[np.ndarray],
    translations: Sequence[np.ndarray],
    tolerance: float,
):
    positions = np.asarray(atoms.get_scaled_positions(wrap=True))
    numbers = np.asarray(atoms.get_atomic_numbers())
    lattice = np.asarray(atoms.cell[:])
    permutations = []
    for rotation, translation in zip(rotations, translations):
        transformed = positions @ rotation.T + translation
        permutation = np.full(len(atoms), -1, dtype=int)
        for source, coordinate in enumerate(transformed):
            candidates = np.where(numbers == numbers[source])[0]
            difference = coordinate - positions[candidates]
            difference -= np.rint(difference)
            distances = np.linalg.norm(difference @ lattice, axis=1)
            best = int(np.argmin(distances))
            if distances[best] > max(5 * tolerance, 1e-5):
                raise SymmetryAnalysisError(
                    "A symmetry operation could not be mapped onto the atoms."
                )
            permutation[source] = int(candidates[best])
        if len(set(permutation.tolist())) != len(atoms):
            raise SymmetryAnalysisError(
                "A detected symmetry operation did not produce an atomic permutation."
            )
        permutations.append(permutation)
    return np.asarray(permutations, dtype=int)


def _representation_matrices(
    selected: Sequence[int],
    basis: str,
    permutations: np.ndarray,
    cart_rotations: Sequence[np.ndarray],
):
    selected = tuple(int(index) for index in selected)
    local_index = {atom_index: i for i, atom_index in enumerate(selected)}
    components = _BASIS_COMPONENTS[basis]
    block_size = len(components)
    dimension = len(selected) * block_size
    matrices = []
    for permutation, rotation in zip(permutations, cart_rotations):
        orbital = _orbital_matrix(basis, rotation)
        matrix = np.zeros((dimension, dimension), dtype=complex)
        for source_local, source_atom in enumerate(selected):
            target_atom = int(permutation[source_atom])
            if target_atom not in local_index:
                raise SymmetryAnalysisError(
                    "The chosen atoms are not closed under the point group. "
                    "Choose a complete symmetry orbit."
                )
            target_local = local_index[target_atom]
            row = slice(target_local * block_size, (target_local + 1) * block_size)
            column = slice(
                source_local * block_size, (source_local + 1) * block_size
            )
            matrix[row, column] = orbital
        matrices.append(matrix)
    return np.asarray(matrices)


def _basis_labels(atoms: Atoms, selected: Sequence[int], basis: str):
    counters: Dict[str, int] = {}
    atom_labels = {}
    for atom_index, symbol in enumerate(atoms.get_chemical_symbols()):
        counters[symbol] = counters.get(symbol, 0) + 1
        atom_labels[atom_index] = f"{symbol}{counters[symbol]}"
    return [
        f"{atom_labels[int(atom_index)]}:{component}"
        for atom_index in selected
        for component in _BASIS_COMPONENTS[basis]
    ]


def _salc_expressions(vectors: np.ndarray, labels: Sequence[str]):
    expressions = []
    for column in range(vectors.shape[1]):
        vector = vectors[:, column].copy()
        pivot = int(np.argmax(np.abs(vector)))
        if abs(vector[pivot]) > 1e-10:
            phase = vector[pivot] / abs(vector[pivot])
            vector /= phase
        terms = []
        coefficients = []
        for coefficient, label in zip(vector, labels):
            coefficient = _clean_number(coefficient)
            coefficients.append(coefficient)
            if abs(coefficient) < 1e-6:
                continue
            terms.append(f"{format_character(coefficient)}·{label}")
        expressions.append({
            "salc": f"SALC {column + 1}",
            "expression": " + ".join(terms).replace("+ -", "− "),
            "coefficients": coefficients,
        })
    return expressions


@dataclass
class _Orbit:
    id: int
    indices: Tuple[int, ...]
    element: str
    wyckoff: str
    site_symmetry: str

    @property
    def label(self):
        return (
            f"{self.element} · {len(self.indices)}{self.wyckoff} · "
            f"site {self.site_symmetry}"
        )


class GammaRepresentationAnalyzer:
    """Analyse structure-derived orbital/displacement representations at Γ."""

    supported_bases = tuple(_BASIS_COMPONENTS)
    basis_titles = dict(_BASIS_TITLES)

    def __init__(self, atoms: Atoms, symprec: float = 1e-3):
        self.symprec = float(symprec)
        self.atoms, self.dataset = _primitive_structure(atoms, self.symprec)
        self.rotations, self.translations = _unique_operations(self.dataset)
        self.cart_rotations = np.asarray([
            _orthogonal_cartesian_rotation(self.atoms.cell.array, rotation)
            for rotation in self.rotations
        ])
        self.table = _multiplication_table(self.rotations)
        _, _, self.classes = _conjugacy_classes(self.table)
        self.class_labels = _class_labels(
            self.classes, self.rotations, self.cart_rotations
        )
        self.complex_irreps = _complex_irreducible_characters(
            self.rotations, self.classes, self.table
        )
        self.real_irreps = _combine_real_irreps(self.complex_irreps)
        _assign_irrep_labels(
            self.pointgroup,
            self.real_irreps,
            self.classes,
            self.cart_rotations,
        )
        self.permutations = _atom_permutations(
            self.atoms,
            self.rotations,
            self.translations,
            self.symprec,
        )
        self.orbits = self._make_orbits()
        self._analysis_cache = {}

    @property
    def pointgroup(self):
        return str(_dataset_value(self.dataset, "pointgroup", "?"))

    @property
    def spacegroup(self):
        return str(_dataset_value(self.dataset, "international", "?"))

    @property
    def spacegroup_number(self):
        return int(_dataset_value(self.dataset, "number", 0))

    @property
    def operation_count(self):
        return len(self.rotations)

    @property
    def class_count(self):
        return len(self.classes)

    @property
    def has_internal_labels(self):
        return any(
            irrep["label_source"] != "conventional" for irrep in self.real_irreps
        )

    def _make_orbits(self):
        equivalent = np.asarray(
            _dataset_value(self.dataset, "equivalent_atoms"), dtype=int
        )
        wyckoffs = list(_dataset_value(self.dataset, "wyckoffs"))
        site_symbols = list(
            _dataset_value(
                self.dataset,
                "site_symmetry_symbols",
                ["?"] * len(self.atoms),
            )
        )
        result = []
        for orbit_id in sorted(set(equivalent.tolist())):
            indices = tuple(np.where(equivalent == orbit_id)[0].tolist())
            atomic_number = int(self.atoms.numbers[indices[0]])
            result.append(_Orbit(
                id=int(orbit_id),
                indices=indices,
                element=chemical_symbols[atomic_number],
                wyckoff=str(wyckoffs[indices[0]]),
                site_symmetry=str(site_symbols[indices[0]]).strip(),
            ))
        return result

    def orbit(self, orbit_id: int):
        for orbit in self.orbits:
            if orbit.id == int(orbit_id):
                return orbit
        raise SymmetryAnalysisError(f"Unknown symmetry orbit: {orbit_id}")

    def _class_characters(self, matrices: np.ndarray):
        traces = np.trace(matrices, axis1=1, axis2=2)
        return np.asarray([
            np.mean(traces[members]) for members in self.classes
        ])

    def _complex_multiplicities(self, class_characters):
        weights = np.asarray([len(members) for members in self.classes])
        result = []
        for irrep in self.complex_irreps:
            value = np.sum(
                weights * class_characters
                * irrep["characters"].conjugate()
            ) / self.operation_count
            if abs(value.imag) > 2e-5 or abs(value.real - round(value.real)) > 2e-5:
                raise SymmetryAnalysisError(
                    "The representation did not decompose into integer irreps. "
                    "Try a slightly different symmetry tolerance."
                )
            result.append(int(round(value.real)))
        return result

    def analyse(self, orbit_id: int, basis: str):
        if basis not in self.supported_bases:
            raise SymmetryAnalysisError(f"Unsupported basis: {basis}")
        cache_key = (int(orbit_id), basis)
        if cache_key in self._analysis_cache:
            return self._analysis_cache[cache_key]
        orbit = self.orbit(orbit_id)
        orbital_traces = np.asarray([
            np.trace(_orbital_matrix(basis, rotation))
            for rotation in self.cart_rotations
        ])
        fixed_counts = np.asarray([
            sum(
                int(permutation[atom] == atom) for atom in orbit.indices
            )
            for permutation in self.permutations
        ])
        operation_characters = fixed_counts * orbital_traces
        class_characters = np.asarray([
            np.mean(operation_characters[members]) for members in self.classes
        ])
        complex_multiplicities = self._complex_multiplicities(class_characters)
        labels = _basis_labels(self.atoms, orbit.indices, basis)

        decomposition = []
        for irrep in self.real_irreps:
            member_multiplicities = [
                complex_multiplicities[index] for index in irrep["members"]
            ]
            multiplicity = member_multiplicities[0]
            if any(value != multiplicity for value in member_multiplicities):
                raise SymmetryAnalysisError(
                    "A real representation contained unequal conjugate irreps."
                )
            if multiplicity == 0:
                continue
            rank = multiplicity * irrep["dimension"]
            decomposition.append({
                "gamma": irrep["gamma_label"],
                "label": irrep["label"],
                "display_label": (
                    irrep["label"]
                    if irrep["label_source"] == "conventional"
                    else f"{irrep['gamma_label']} ({irrep['label']}-like)"
                ),
                "label_source": irrep["label_source"],
                "dimension": irrep["dimension"],
                "multiplicity": multiplicity,
                "states": rank,
                "basis_hint": self._basis_hint(irrep),
                "members": irrep["members"],
            })

        represented_dimension = sum(
            row["states"] for row in decomposition
        )
        dimension = len(orbit.indices) * len(_BASIS_COMPONENTS[basis])
        if represented_dimension != dimension:
            raise SymmetryAnalysisError(
                "Irrep dimensions do not reproduce the selected basis dimension."
            )

        character_rows = []
        for class_index, (label, members) in enumerate(
            zip(self.class_labels, self.classes)
        ):
            representative = members[0]
            fixed = sum(
                int(self.permutations[representative, atom] == atom)
                for atom in orbit.indices
            )
            character_rows.append({
                "class": label,
                "operations": len(members),
                "fixed atoms": fixed,
                "χ(reducible)": format_character(
                    class_characters[class_index]
                ),
            })

        result = {
            "orbit": orbit,
            "basis": basis,
            "basis_title": self.basis_titles[basis],
            "dimension": dimension,
            "class_characters": class_characters,
            "complex_multiplicities": complex_multiplicities,
            "decomposition": decomposition,
            "character_rows": character_rows,
            "basis_labels": labels,
            "salcs": {},
        }
        self._analysis_cache[cache_key] = result
        return result

    def generate_salcs(self, result):
        """Generate projection-operator SALCs lazily for an analysis result."""
        if result["salcs"]:
            return result["salcs"]
        orbit = result["orbit"]
        basis = result["basis"]
        matrices = _representation_matrices(
            orbit.indices, basis, self.permutations, self.cart_rotations
        )
        salcs = {}
        for row in result["decomposition"]:
            irrep = next(
                candidate for candidate in self.real_irreps
                if candidate["gamma_label"] == row["gamma"]
            )
            projector = np.zeros_like(matrices[0])
            for member in irrep["members"]:
                complex_irrep = self.complex_irreps[member]
                projector += (
                    complex_irrep["dimension"] / self.operation_count
                    * sum(
                        complex_irrep["characters"][class_index].conjugate()
                        * matrices[operation]
                        for class_index, members in enumerate(self.classes)
                        for operation in members
                    )
                )
            projector = (projector + projector.conjugate().T) / 2.0
            values, vectors = np.linalg.eigh(projector)
            rank = row["states"]
            selected_vectors = vectors[:, np.argsort(values)[-rank:]]
            salcs[row["gamma"]] = _salc_expressions(
                selected_vectors, result["basis_labels"]
            )
        result["salcs"] = salcs
        return salcs

    @staticmethod
    def _basis_hint(irrep):
        hints = []
        if irrep.get("polar_vector"):
            hints.append("polar vector: x, y, z / p")
        if irrep.get("axial_vector"):
            hints.append("axial vector: Rₓ, Rᵧ, R_z")
        if np.allclose(irrep["characters"], 1.0, atol=1e-6):
            hints.append("fully symmetric / s-like")
        return "; ".join(hints) or "symmetry-adapted combination"

    @staticmethod
    def decomposition_text(result):
        terms = []
        for row in result["decomposition"]:
            prefix = f"{row['multiplicity']}×" if row["multiplicity"] > 1 else ""
            terms.append(f"{prefix}{row['label']}")
        return " ⊕ ".join(terms)

    def compatibility(
        self,
        source_result,
        bases: Iterable[str] = ("s", "p", "d"),
    ):
        source = {
            row["gamma"]: row["multiplicity"]
            for row in source_result["decomposition"]
        }
        rows = []
        for orbit in self.orbits:
            for basis in bases:
                candidate = self.analyse(orbit.id, basis)
                shared = []
                channels = 0
                for row in candidate["decomposition"]:
                    if row["gamma"] in source:
                        shared.append(row["display_label"])
                        channels += min(
                            source[row["gamma"]], row["multiplicity"]
                        ) * row["dimension"]
                rows.append({
                    "target orbit": orbit.label,
                    "target basis": self.basis_titles[basis],
                    "shared irreps": ", ".join(shared) if shared else "—",
                    "compatible channels": channels,
                    "symmetry allowed": "Yes" if shared else "No",
                })
        return rows

    def character_table_rows(self):
        rows = []
        for irrep in self.real_irreps:
            row = {
                "irrep": irrep["label"],
                "Γ label": irrep["gamma_label"],
                "dimension": irrep["dimension"],
            }
            row.update({
                label: format_character(character)
                for label, character in zip(
                    self.class_labels, irrep["characters"]
                )
            })
            rows.append(row)
        return rows
