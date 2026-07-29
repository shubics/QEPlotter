import unittest

import numpy as np
from ase import Atoms

from gui.page_structure import _formula_identity
from qeplotter.analysis.bilayer import (
    analyse_stacking, classify_stacking, detect_bilayer, ibrav2cell,
    parse_qe_block,
)
from qeplotter.core.utils import BOHR_TO_ANGSTROM
from qeplotter.core.utils import strip_number
from qeplotter.structure import analyse_angles, analyse_bonds, build_3dmol_html
from qeplotter.structure.viz import element_colors


class StructureAnalysisTests(unittest.TestCase):
    @staticmethod
    def _tmd_bilayer(upper_offset, upper_orientation=1, janus=False,
                     upper_metal="Mo"):
        """Small commensurate MX2-like bilayer for registry regression tests."""
        cell = np.array([
            [3.0, 0.0, 0.0],
            [-1.5, 3.0 * np.sqrt(3) / 2, 0.0],
            [0.0, 0.0, 20.0],
        ])
        q = np.array([1 / 3, 2 / 3])
        upper_q = q if upper_orientation == 1 else -q
        lower_symbols = ["S", "Mo", "Se" if janus else "S"]
        upper_symbols = ["S" if janus else "S", upper_metal, "Se" if janus else "S"]
        fractional = np.array([
            [*q, 0.20], [0.0, 0.0, 0.25], [*q, 0.30],
            [*(upper_offset + upper_q), 0.50],
            [*upper_offset, 0.55],
            [*(upper_offset + upper_q), 0.60],
        ], dtype=float)
        fractional[:, :2] %= 1.0
        return cell, lower_symbols + upper_symbols, fractional

    def test_all_six_tmd_stacking_registries_are_distinguished(self):
        q = np.array([1 / 3, 2 / 3])
        cases = {
            "AA": (np.array([0.0, 0.0]), 1),
            "AB": (-q, 1),
            "BA": (q, 1),
            "AA′": (q, -1),
            "AB′": (np.array([0.0, 0.0]), -1),
            "A′B": (2 * q, -1),
        }
        for expected, (offset, orientation) in cases.items():
            with self.subTest(stacking=expected):
                cell, species, fractional = self._tmd_bilayer(
                    offset, orientation)
                self.assertEqual(
                    classify_stacking(cell, species, fractional), expected)

    def test_general_registry_is_not_mislabelled_as_aa(self):
        cell, species, fractional = self._tmd_bilayer(
            np.array([0.17, 0.11]), 1)
        result = analyse_stacking(cell, species, fractional)
        self.assertEqual(result["label"], "General registry")
        self.assertEqual(result["confidence"], "safe fallback")
        self.assertIsNotNone(result["shift"])

    def test_janus_interface_and_layer_order_are_reported(self):
        q = np.array([1 / 3, 2 / 3])
        cell, species, fractional = self._tmd_bilayer(
            q, -1, janus=True, upper_metal="W")
        result = analyse_stacking(cell, species, fractional)
        self.assertEqual(result["label"], "AA′")
        self.assertEqual(result["bilayer_type"], "Janus heterobilayer")
        self.assertEqual(result["lower_formula"], "MoSSe")
        self.assertEqual(result["upper_formula"], "WSSe")
        self.assertEqual(result["interface"], "Se | S")

    def test_tolerance_changes_bonds_and_viewer(self):
        atoms = Atoms("CO", positions=[[0, 0, 0], [0, 0, 1.7]],
                      cell=[8, 8, 8], pbc=True)
        low, _ = analyse_bonds(atoms, tol=0.8, include_periodic=False)
        high, _ = analyse_bonds(atoms, tol=1.6, include_periodic=False)
        self.assertEqual((len(low), len(high)), (0, 1))
        self.assertNotEqual(
            build_3dmol_html(atoms, bond_tol=0.8),
            build_3dmol_html(atoms, bond_tol=1.6),
        )

    def test_periodic_bonds_are_optional_and_explicit(self):
        atoms = Atoms("Li", positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=True)
        local, _ = analyse_bonds(atoms, tol=1.0, include_periodic=False)
        periodic, geometry = analyse_bonds(
            atoms, tol=1.0, include_periodic=True)
        self.assertEqual(len(local), 0)
        self.assertEqual(len(periodic), 3)
        self.assertTrue(all(any(record["offset"]) for record in geometry))
        self.assertTrue(all("[" in label for label in periodic["atom_2"]))

    def test_bilayer_and_monolayer_are_distinguished(self):
        cell = np.diag([3.0, 3.0, 20.0])
        bilayer_frac = np.array([
            [0, 0, 0.20], [0.5, 0.5, 0.22],
            [0, 0, 0.40], [0.5, 0.5, 0.42],
        ])
        monolayer_frac = np.array([
            [0, 0, 0.45], [0.5, 0.5, 0.50], [0, 0, 0.55],
        ])
        self.assertTrue(
            detect_bilayer(cell, ["W", "S", "W", "Se"], bilayer_frac)["is_bilayer"])
        self.assertFalse(
            detect_bilayer(cell, ["S", "W", "S"], monolayer_frac)["is_bilayer"])

    def test_angle_rows_align_with_geometry(self):
        atoms = Atoms("HCO", positions=[[0, 0, 0], [0, 0, 1], [1, 0, 1]],
                      cell=[4, 4, 4], pbc=True)
        frame, geometry = analyse_angles(atoms, tol=1.6)
        self.assertEqual(len(frame), len(geometry))
        keys = list(zip(frame.vertex, frame.neighbor_1, frame.neighbor_2))
        self.assertEqual(len(keys), len(set(keys)))
        self.assertTrue(all(len(record["points"]) == 3 for record in geometry))

    def test_equivalent_angle_values_are_not_removed(self):
        atoms = Atoms("Li", positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=True)
        frame, _ = analyse_angles(atoms, tol=1.0, include_periodic=True)
        self.assertTrue(frame["angle (°)"].duplicated().any())

    def test_element_colours_never_repeat(self):
        symbols = ["H", "He", "Li", "Be", "B", "C", "N", "O",
                   "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S"]
        atoms = Atoms(symbols, positions=np.zeros((len(symbols), 3)),
                      cell=[10, 10, 10])
        colours = element_colors(atoms)
        self.assertEqual(len(colours), len(set(colours.values())))

    def test_qe_alat_positions_are_cartesian_and_expressions_are_supported(self):
        cell, species, fractional = parse_qe_block([
            "&SYSTEM",
            "  ibrav = 4, celldm(1) = 10, celldm(3) = 2, nat = 1",
            "/",
            "ATOMIC_POSITIONS (alat)",
            "H 0 1/2 0",
        ])
        self.assertEqual(species, ["H"])
        cartesian = fractional @ cell
        np.testing.assert_allclose(
            cartesian[0], [0, 5 * BOHR_TO_ANGSTROM, 0], atol=1e-10
        )
        self.assertFalse(np.allclose(fractional[0], [0, 0.5, 0]))

    def test_qe_parenthesized_cell_units_and_conventional_parameters(self):
        cell, species, fractional = parse_qe_block([
            "&SYSTEM",
            "  ibrav = 0, A = 4.0, nat = 1",
            "/",
            "CELL_PARAMETERS (angstrom)",
            "4 0 0",
            "0 5 0",
            "0 0 6",
            "ATOMIC_POSITIONS (angstrom)",
            "Si 2 2.5 3",
        ])
        np.testing.assert_allclose(cell, np.diag([4.0, 5.0, 6.0]))
        np.testing.assert_allclose(fractional, [[0.5, 0.5, 0.5]])
        self.assertEqual(species, ["Si"])

    def test_qe_monoclinic_celldm_is_a_cosine_not_an_angle(self):
        cell = ibrav2cell(12, {1: 10.0, 2: 1.2, 3: 1.4, 4: 0.5})
        cosine_ab = np.dot(cell[0], cell[1]) / (
            np.linalg.norm(cell[0]) * np.linalg.norm(cell[1])
        )
        self.assertAlmostEqual(cosine_ab, 0.5)

    def test_supported_qe_ibrav_variants_have_nonzero_cells(self):
        common = {1: 10.0, 2: 1.2, 3: 1.4, 4: 0.2, 5: 0.25, 6: 0.3}
        for ibrav in (1, 2, 3, -3, 4, 5, -5, 6, 7, 8, 9, -9, 91,
                      10, 11, 12, -12, 13, -13, 14):
            with self.subTest(ibrav=ibrav):
                self.assertGreater(abs(np.linalg.det(ibrav2cell(ibrav, common))), 1e-8)

    def test_qe_species_suffixes_resolve_to_elements(self):
        self.assertEqual(strip_number("Fe1"), "Fe")
        self.assertEqual(strip_number("Fe_up"), "Fe")
        self.assertEqual(strip_number("C-h"), "C")
        self.assertEqual(strip_number("se2"), "Se")

    def test_heterobilayer_formula_preserves_each_layer(self):
        atoms = Atoms(
            ["W", "Se", "W", "S"],
            scaled_positions=[
                [0, 0, 0.20], [0.3, 0.3, 0.22],
                [0, 0, 0.40], [0.3, 0.3, 0.42],
            ],
            cell=np.diag([3.0, 3.0, 20.0]),
            pbc=True,
        )
        identity = _formula_identity(atoms)
        self.assertEqual(identity["display"], "WSe / WS")
        self.assertEqual(identity["layers"], ("WSe", "WS"))


if __name__ == "__main__":
    unittest.main()
