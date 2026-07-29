import unittest

import numpy as np
import spglib
from ase import Atoms
from ase.build import bulk

from qeplotter.symmetry import GammaRepresentationAnalyzer
from qeplotter.symmetry.representations import (
    _combine_real_irreps,
    _complex_irreducible_characters,
    _conjugacy_classes,
    _multiplication_table,
)


def methane_crystal():
    """CH4 in a cubic P-43m cell: a compact Td regression structure."""
    offset = 0.12
    positions = np.asarray([
        [0.5, 0.5, 0.5],
        [0.5 + offset, 0.5 + offset, 0.5 + offset],
        [0.5 + offset, 0.5 - offset, 0.5 - offset],
        [0.5 - offset, 0.5 + offset, 0.5 - offset],
        [0.5 - offset, 0.5 - offset, 0.5 + offset],
    ])
    return Atoms(
        "CH4",
        scaled_positions=positions,
        cell=np.eye(3) * 8.0,
        pbc=True,
    )


def mos2_2h_monolayer():
    lattice_a = 3.18
    vacuum_c = 20.0
    sulfur_height = 0.078
    cell = np.asarray([
        [lattice_a, 0.0, 0.0],
        [-lattice_a / 2.0, np.sqrt(3.0) * lattice_a / 2.0, 0.0],
        [0.0, 0.0, vacuum_c],
    ])
    return Atoms(
        "MoS2",
        scaled_positions=[
            [0.0, 0.0, 0.5],
            [1 / 3, 2 / 3, 0.5 + sulfur_height],
            [1 / 3, 2 / 3, 0.5 - sulfur_height],
        ],
        cell=cell,
        pbc=True,
    )


class SymmetryRepresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analyzer = GammaRepresentationAnalyzer(methane_crystal())
        cls.carbon = next(
            orbit for orbit in cls.analyzer.orbits if orbit.element == "C"
        )
        cls.hydrogen = next(
            orbit for orbit in cls.analyzer.orbits if orbit.element == "H"
        )

    def test_complete_td_group_is_used(self):
        self.assertEqual(self.analyzer.pointgroup, "-43m")
        self.assertEqual(self.analyzer.operation_count, 24)
        self.assertEqual(sum(map(len, self.analyzer.classes)), 24)
        self.assertTrue(any("C3" in label for label in self.analyzer.class_labels))

    def test_four_ligand_s_orbitals_reduce_to_a1_plus_t2(self):
        result = self.analyzer.analyse(self.hydrogen.id, "s")
        components = {
            row["label"]: row["multiplicity"]
            for row in result["decomposition"]
        }
        self.assertEqual(components, {"A1": 1, "T2": 1})
        self.assertEqual(result["dimension"], 4)
        self.assertEqual(sum(row["states"] for row in result["decomposition"]), 4)

    def test_central_orbitals_match_the_ligand_channels(self):
        carbon_s = self.analyzer.analyse(self.carbon.id, "s")
        carbon_p = self.analyzer.analyse(self.carbon.id, "p")
        carbon_d = self.analyzer.analyse(self.carbon.id, "d")
        self.assertEqual(
            [row["label"] for row in carbon_s["decomposition"]], ["A1"]
        )
        self.assertEqual(
            [row["label"] for row in carbon_p["decomposition"]], ["T2"]
        )
        self.assertEqual(
            {row["label"] for row in carbon_d["decomposition"]}, {"E", "T2"}
        )

    def test_fully_symmetric_salc_is_normalised_and_equal_weight(self):
        result = self.analyzer.analyse(self.hydrogen.id, "s")
        self.analyzer.generate_salcs(result)
        a1 = next(
            row for row in result["decomposition"] if row["label"] == "A1"
        )
        coefficients = np.asarray(
            result["salcs"][a1["gamma"]][0]["coefficients"], dtype=complex
        )
        self.assertTrue(np.allclose(coefficients, 0.5))
        self.assertAlmostEqual(float(np.vdot(coefficients, coefficients).real), 1.0)

    def test_compatibility_reports_both_carbon_channels(self):
        ligand = self.analyzer.analyse(self.hydrogen.id, "s")
        rows = self.analyzer.compatibility(ligand, bases=("s", "p"))
        carbon_rows = {
            row["target basis"]: row for row in rows
            if row["target orbit"].startswith("C ")
        }
        self.assertIn("A1", carbon_rows["s orbitals"]["shared irreps"])
        self.assertIn("T2", carbon_rows["p orbitals (pₓ, pᵧ, p_z)"]["shared irreps"])

    def test_finite_group_engine_covers_all_32_crystallographic_point_groups(self):
        checked = set()
        for hall_number in range(1, 531):
            symmetry = spglib.get_symmetry_from_database(hall_number)
            rotations = []
            seen_rotations = set()
            for rotation in symmetry["rotations"]:
                key = tuple(rotation.ravel())
                if key not in seen_rotations:
                    seen_rotations.add(key)
                    rotations.append(np.asarray(rotation, dtype=int))
            pointgroup = spglib.get_pointgroup(rotations)[0].strip()
            if pointgroup in checked:
                continue
            identity = next(
                i for i, rotation in enumerate(rotations)
                if np.array_equal(rotation, np.eye(3, dtype=int))
            )
            order = [identity] + sorted(
                (i for i in range(len(rotations)) if i != identity),
                key=lambda i: tuple(rotations[i].ravel()),
            )
            rotations = [rotations[i] for i in order]
            table = _multiplication_table(rotations)
            _, _, classes = _conjugacy_classes(table)
            complex_irreps = _complex_irreducible_characters(
                rotations, classes, table
            )
            real_irreps = _combine_real_irreps(complex_irreps)
            self.assertEqual(
                sum(irrep["dimension"] ** 2 for irrep in complex_irreps),
                len(rotations),
            )
            self.assertTrue(real_irreps)
            checked.add(pointgroup)
        self.assertEqual(len(checked), 32)

    def test_real_materials_across_crystal_families(self):
        materials = {
            "diamond Si": (
                bulk("Si", "diamond", a=5.43), "Fd-3m", "m-3m"
            ),
            "rocksalt NaCl": (
                bulk("NaCl", "rocksalt", a=5.64), "Fm-3m", "m-3m"
            ),
            "hcp Mg": (
                bulk("Mg", "hcp", a=3.21, c=5.21), "P6_3/mmc", "6/mmm"
            ),
            "wurtzite ZnO": (
                bulk("ZnO", "wurtzite", a=3.25, c=5.21, u=0.382),
                "P6_3mc",
                "6mm",
            ),
            "2H monolayer MoS2": (
                mos2_2h_monolayer(), "P-6m2", "-6m2"
            ),
        }
        component_count = {"s": 1, "p": 3, "d": 5}
        for name, (atoms, expected_spacegroup, expected_pointgroup) in materials.items():
            with self.subTest(material=name):
                analyzer = GammaRepresentationAnalyzer(atoms)
                self.assertEqual(analyzer.spacegroup, expected_spacegroup)
                self.assertEqual(analyzer.pointgroup, expected_pointgroup)
                for orbit in analyzer.orbits:
                    for basis, components in component_count.items():
                        result = analyzer.analyse(orbit.id, basis)
                        expected_dimension = len(orbit.indices) * components
                        self.assertEqual(result["dimension"], expected_dimension)
                        self.assertEqual(
                            sum(row["states"] for row in result["decomposition"]),
                            expected_dimension,
                        )


if __name__ == "__main__":
    unittest.main()
