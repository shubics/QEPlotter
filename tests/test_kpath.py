import unittest
from types import SimpleNamespace

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.dft.kpoints import parse_path_string
from ase.lattice import MCLC, all_variants

from qeplotter.kpath import (
    build_bz_figure, format_qe_kpoints, parse_path_expression, primary_path,
    recommend_kpath, with_path,
)
from qeplotter.kpath.recipes_sc import get_sc_recipe


class NativeKPathTests(unittest.TestCase):
    def test_common_bravais_recipes(self):
        samples = [
            ("cP", bulk("Cu", "sc", a=3.6)),
            ("cF", bulk("Cu", "fcc", a=3.6)),
            ("cI", bulk("Fe", "bcc", a=2.86)),
            ("hP", bulk("Mg", "hcp", a=3.2, c=5.2)),
            ("tP", Atoms("Si", scaled_positions=[[0, 0, 0]],
                          cell=np.diag([3.0, 3.0, 5.0]), pbc=True)),
            ("oP", Atoms("Si", scaled_positions=[[0, 0, 0]],
                          cell=np.diag([3.0, 4.0, 5.0]), pbc=True)),
        ]
        for expected, atoms in samples:
            with self.subTest(bravais=expected):
                result = recommend_kpath(atoms, reference_distance=0.2)
                self.assertEqual(result["bravais_lattice"], expected)
                self.assertEqual(result["method"], "conventional-recipe")
                self.assertIn("GAMMA", result["point_coords"])
                self.assertGreaterEqual(len(result["bz"]["vertices"]), 8)

    def test_low_symmetry_uses_conventional_triclinic_recipe(self):
        atoms = Atoms("Si", scaled_positions=[[0, 0, 0]],
                      cell=[[3.1, 0.2, 0.1], [0.4, 4.2, 0.3],
                            [0.2, 0.5, 5.3]], pbc=True)
        result = recommend_kpath(atoms, reference_distance=0.2)
        self.assertEqual(result["method"], "conventional-recipe")
        self.assertTrue(result["recipe_variant"].startswith("TRI"))
        self.assertFalse(any(label.startswith("V")
                             for label in result["point_coords"]))

    def test_ti_variants_and_metric_parameters(self):
        cases = [(2.0, "BCT1", "tI1"), (5.0, "BCT2", "tI2")]
        for c, variant, extended in cases:
            atoms = Atoms("Si2", scaled_positions=[[0, 0, 0], [.5, .5, .5]],
                          cell=np.diag([3.0, 3.0, c]), pbc=True)
            result = recommend_kpath(atoms)
            self.assertEqual(result["recipe_variant"], variant)
            self.assertEqual(result["bravais_lattice_extended"], extended)
            self.assertIn("eta", result["recipe_parameters"])
            if variant == "BCT2":
                self.assertIn("zeta", result["recipe_parameters"])

    def test_all_sc_reference_variants(self):
        ignored_2d = {"OBL", "RECT", "CRECT", "HEX2D", "SQR", "LINE"}
        checked = set()
        for lattice in all_variants():
            if lattice.name in ignored_2d:
                continue
            recipe = get_sc_recipe(lattice)
            reference = lattice.get_special_points()
            for label, coords in reference.items():
                key = "GAMMA" if label == "G" else label
                np.testing.assert_allclose(recipe["points"][key], coords)
            expected_path = []
            for branch in parse_path_string(lattice.special_path):
                expected_path.extend(zip(branch, branch[1:]))
            expected_path = [("GAMMA" if a == "G" else a,
                              "GAMMA" if b == "G" else b)
                             for a, b in expected_path]
            self.assertEqual(recipe["path"], expected_path)
            checked.add(lattice.variant)
        # ASE's example generator lacks boundary variants MCLC2/MCLC4; exercise
        # those registry branches explicitly (their coordinate formulas share
        # MCLC1/MCLC3 but their recommended paths differ).
        base1 = MCLC(2, 3, 3, 50)
        base4 = MCLC(2, 2, 4, 60)
        for variant, base in (("MCLC2", base1), ("MCLC4", base4)):
            proxy = SimpleNamespace(name="MCLC", variant=variant,
                                    vars=base.vars)
            recipe = get_sc_recipe(proxy)
            self.assertEqual(recipe["variant"], variant)
            self.assertTrue(recipe["path"])
            checked.add(variant)
        self.assertGreaterEqual(len(checked), 25)

    def test_uploaded_reciprocal_basis_is_explicit_and_consistent(self):
        atoms = bulk("Cu", "fcc", a=3.6)
        result = recommend_kpath(atoms)
        for label, fractional in result["point_coords"].items():
            reconstructed = fractional @ result["reciprocal_input"]
            np.testing.assert_allclose(reconstructed,
                                       result["point_coords_cart"][label], atol=1e-10)
        self.assertIn("2π", result["reciprocal_convention"])
        self.assertLess(result["basis_mapping_error"], 1e-10)

    def test_engine_covers_all_reference_bravais_variants(self):
        ignored_2d = {"OBL", "RECT", "CRECT", "HEX2D", "SQR", "LINE"}
        covered_lattices = set()
        for lattice in all_variants():
            if lattice.name in ignored_2d:
                continue
            atoms = Atoms("H", scaled_positions=[[0, 0, 0]],
                          cell=lattice.tocell(), pbc=True)
            result = recommend_kpath(atoms, symprec=1e-5,
                                     reference_distance=0.2)
            self.assertEqual(result["method"], "conventional-recipe")
            self.assertFalse(any(label.startswith("V")
                                 for label in result["point_coords"]))
            covered_lattices.add(result["bravais_lattice"])
        self.assertEqual(len(covered_lattices), 14)

    def test_global_rotation_preserves_recipe_and_path_lengths(self):
        atoms = bulk("Cu", "fcc", a=3.6)
        original = recommend_kpath(atoms, reference_distance=0.2)
        rotated = atoms.copy()
        rotated.rotate(37.0, "z", rotate_cell=True)
        rotated.rotate(23.0, "x", rotate_cell=True)
        transformed = recommend_kpath(rotated, reference_distance=0.2)
        self.assertEqual(original["recipe_variant"], transformed["recipe_variant"])

        def lengths(result):
            return sorted(round(float(np.linalg.norm(
                result["point_coords_cart"][end] -
                result["point_coords_cart"][start])), 10)
                for start, end in result["path"])

        self.assertEqual(lengths(original), lengths(transformed))

    def test_qe_card_and_plotly_figure(self):
        result = recommend_kpath(bulk("Cu", "fcc", a=3.6),
                                 reference_distance=0.2)
        qe_card = format_qe_kpoints(result)
        lines = qe_card.splitlines()
        self.assertEqual(lines[0], "K_POINTS crystal_b")
        self.assertEqual(int(lines[1]), 2 * len(result["path"]))
        self.assertEqual(len(lines) - 2, 2 * len(result["path"]))
        explicit = format_qe_kpoints(result, explicit=True).splitlines()
        self.assertEqual(explicit[0], "K_POINTS crystal")
        self.assertEqual(int(explicit[1]), len(result["explicit"]))
        figure = build_bz_figure(result)
        self.assertGreaterEqual(len(figure.data), 3)

    def test_smaller_spacing_produces_more_points(self):
        atoms = bulk("Fe", "bcc", a=2.86)
        coarse = recommend_kpath(atoms, reference_distance=0.2)
        fine = recommend_kpath(atoms, reference_distance=0.05)
        self.assertGreater(len(fine["explicit"]), len(coarse["explicit"]))

    def test_primary_and_custom_path_variants(self):
        result = recommend_kpath(bulk("Cu", "sc", a=3.6),
                                 reference_distance=0.2)
        primary = primary_path(result["path"])
        self.assertLess(len(primary), len(result["path"]))
        primary_result = with_path(result, primary)
        self.assertEqual(primary_result["path"], primary)
        self.assertLess(len(primary_result["explicit"]), len(result["explicit"]))

        custom = parse_path_expression("Γ-X-M-Γ | R-X",
                                       result["point_coords"].keys())
        self.assertEqual(custom, [("GAMMA", "X"), ("X", "M"),
                                  ("M", "GAMMA"), ("R", "X")])
        custom_result = with_path(result, custom)
        self.assertEqual(custom_result["path"], custom)
        self.assertIn("K_POINTS crystal_b", format_qe_kpoints(custom_result))

    def test_custom_path_rejects_unknown_points(self):
        result = recommend_kpath(bulk("Cu", "sc", a=3.6))
        with self.assertRaisesRegex(ValueError, "Unknown point"):
            parse_path_expression("Γ-NOT_A_POINT",
                                  result["point_coords"].keys())


if __name__ == "__main__":
    unittest.main()
