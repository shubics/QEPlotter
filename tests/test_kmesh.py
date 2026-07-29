import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.build import bulk, graphene
from ase.lattice import all_variants

from qeplotter.kmesh import (
    build_kmesh_figure,
    format_qe_automatic,
    format_qe_ir_kpoints,
    full_grid_points,
    irreducible_kmesh,
    orbit_members,
)


class NativeIrreducibleKMeshTests(unittest.TestCase):
    def test_reducer_does_not_call_spglib_reciprocal_mesh(self):
        with patch(
            "spglib.get_ir_reciprocal_mesh",
            side_effect=AssertionError("forbidden reciprocal-mesh shortcut"),
        ):
            result = irreducible_kmesh(
                bulk("Cu", "sc", a=3.6), (2, 2, 2))
        self.assertEqual(result["irreducible_count"], 4)
        self.assertEqual(result["engine"],
                         "QEPlotter native integer-orbit reducer")

    def test_simple_cubic_gamma_mesh_has_known_orbits(self):
        result = irreducible_kmesh(
            bulk("Cu", "sc", a=3.6), (2, 2, 2))
        self.assertEqual(result["total_grid_points"], 8)
        self.assertEqual(result["irreducible_count"], 4)
        self.assertEqual(
            [point["multiplicity"] for point in result["points"]],
            [1, 3, 3, 1],
        )
        self.assertEqual(result["compatible_spatial_rotations"], 48)
        self.assertTrue(result["full_crystal_symmetry_preserved"])

    def test_shifted_cubic_mesh_coordinates_and_weight(self):
        result = irreducible_kmesh(
            bulk("Cu", "sc", a=3.6),
            (2, 2, 2),
            shift=(1, 1, 1),
        )
        self.assertEqual(result["irreducible_count"], 1)
        point = result["points"][0]
        self.assertEqual(point["multiplicity"], 8)
        np.testing.assert_allclose(point["frac"], (0.25, 0.25, 0.25))
        self.assertAlmostEqual(point["normalized_weight"], 1.0)

    def test_time_reversal_is_an_explicit_physical_choice(self):
        atoms = Atoms(
            "SiO",
            scaled_positions=[[0.123, 0.217, 0.319],
                              [0.371, 0.113, 0.429]],
            cell=[[3.1, 0.2, 0.1],
                  [0.4, 4.2, 0.3],
                  [0.2, 0.5, 5.3]],
            pbc=True,
        )
        without = irreducible_kmesh(
            atoms, (3, 3, 2), time_reversal=False, symprec=1e-7)
        with_reversal = irreducible_kmesh(
            atoms, (3, 3, 2), time_reversal=True, symprec=1e-7)
        self.assertEqual(without["spacegroup_number"], 1)
        self.assertEqual(without["irreducible_count"], 18)
        self.assertEqual(with_reversal["irreducible_count"], 10)

    def test_anisotropic_grid_uses_only_grid_stabiliser(self):
        result = irreducible_kmesh(
            bulk("Cu", "sc", a=3.6), (4, 6, 4))
        self.assertEqual(result["detected_unique_point_rotations"], 48)
        self.assertEqual(result["compatible_spatial_rotations"], 16)
        self.assertEqual(result["dropped_spatial_rotations"], 32)
        self.assertFalse(result["full_crystal_symmetry_preserved"])
        self.assertEqual(
            sum(point["multiplicity"] for point in result["points"]), 96)

    def test_2d_mesh_is_supported(self):
        atoms = graphene(formula="C2", a=2.46, thickness=0.0, vacuum=10.0)
        result = irreducible_kmesh(atoms, (6, 6, 1))
        self.assertEqual(result["total_grid_points"], 36)
        self.assertLess(result["irreducible_count"], 36)
        self.assertEqual(
            sum(point["multiplicity"] for point in result["points"]), 36)

    def test_global_translation_does_not_change_reduction(self):
        atoms = bulk("Si", "diamond", a=5.43)
        original = irreducible_kmesh(
            atoms, (4, 4, 4), shift=(1, 0, 1))
        translated = atoms.copy()
        translated.set_scaled_positions(
            translated.get_scaled_positions() + [0.173, 0.291, 0.407])
        moved = irreducible_kmesh(
            translated, (4, 4, 4), shift=(1, 0, 1))
        self.assertEqual(original["spacegroup_number"],
                         moved["spacegroup_number"])
        self.assertEqual(
            [point["multiplicity"] for point in original["points"]],
            [point["multiplicity"] for point in moved["points"]],
        )

    def test_every_full_point_maps_to_one_orbit(self):
        result = irreducible_kmesh(
            bulk("Fe", "bcc", a=2.86), (5, 4, 3), shift=(1, 0, 1))
        complete = full_grid_points(result)
        self.assertEqual(len(complete), 60)
        mapped = [point["irreducible_index"] for point in complete]
        self.assertEqual(min(mapped), 1)
        self.assertEqual(max(mapped), result["irreducible_count"])
        for index, point in enumerate(result["points"]):
            members = orbit_members(result, index)
            self.assertEqual(len(members), point["multiplicity"])

    def test_weights_and_qe_exports_are_consistent(self):
        result = irreducible_kmesh(
            bulk("Cu", "fcc", a=3.6), (4, 4, 4), shift=(1, 1, 1))
        self.assertAlmostEqual(
            sum(point["normalized_weight"] for point in result["points"]),
            1.0,
            places=14,
        )
        automatic = format_qe_automatic(result).splitlines()
        self.assertEqual(automatic[0], "K_POINTS automatic")
        self.assertEqual(automatic[1].split(), ["4", "4", "4", "1", "1", "1"])

        normalized = format_qe_ir_kpoints(result)
        multiplicity = format_qe_ir_kpoints(
            result, weight_mode="multiplicity")
        self.assertEqual(
            int(normalized.splitlines()[1]), result["irreducible_count"])
        normalized_sum = sum(
            float(line.split()[3])
            for line in normalized.splitlines()[2:])
        multiplicity_sum = sum(
            float(line.split()[3])
            for line in multiplicity.splitlines()[2:])
        self.assertAlmostEqual(normalized_sum, 1.0, places=10)
        self.assertEqual(multiplicity_sum, 64.0)

        without_reversal = irreducible_kmesh(
            bulk("Cu", "fcc", a=3.6), (2, 2, 2),
            time_reversal=False)
        self.assertIn(
            "noinv=.true.", format_qe_automatic(without_reversal))
        self.assertIn(
            "noinv=.true.", format_qe_ir_kpoints(without_reversal))

    def test_reciprocal_convention_and_plot_are_explicit(self):
        result = irreducible_kmesh(
            bulk("Mg", "hcp", a=3.2, c=5.2), (4, 4, 3))
        self.assertIn("2π", result["reciprocal_convention"])
        first = result["points"][0]
        reconstructed = first["frac"] @ result["reciprocal_input"]
        self.assertEqual(reconstructed.shape, (3,))
        figure = build_kmesh_figure(result, show_full_grid=True)
        self.assertGreaterEqual(len(figure.data), 3)

    def test_all_three_dimensional_bravais_variants_satisfy_invariants(self):
        ignored_2d = {"OBL", "RECT", "CRECT", "HEX2D", "SQR", "LINE"}
        checked = 0
        for lattice in all_variants():
            if lattice.name in ignored_2d:
                continue
            with self.subTest(lattice=lattice.name,
                              variant=lattice.variant):
                atoms = Atoms(
                    "H", scaled_positions=[[0, 0, 0]],
                    cell=lattice.tocell(), pbc=True)
                result = irreducible_kmesh(
                    atoms, (3, 2, 2), symprec=1e-5)
                self.assertEqual(
                    sum(point["multiplicity"]
                        for point in result["points"]),
                    12,
                )
                self.assertAlmostEqual(
                    sum(point["normalized_weight"]
                        for point in result["points"]),
                    1.0,
                    places=14,
                )
                checked += 1
        self.assertGreaterEqual(checked, 23)

    def test_invalid_inputs_and_safety_limit_are_rejected(self):
        atoms = bulk("Cu", "sc", a=3.6)
        invalid = [
            ((4, 0, 4), (0, 0, 0), "positive"),
            ((4.5, 4, 4), (0, 0, 0), "integers"),
            ((4, 4), (0, 0, 0), "exactly three"),
            ((4, 4, 4), (0, 2, 0), "0 or 1"),
        ]
        for mesh, shift, message in invalid:
            with self.subTest(mesh=mesh, shift=shift):
                with self.assertRaisesRegex(ValueError, message):
                    irreducible_kmesh(atoms, mesh, shift=shift)
        with self.assertRaisesRegex(ValueError, "safe limit"):
            irreducible_kmesh(
                atoms, (100, 100, 100), max_grid_points=1000)


if __name__ == "__main__":
    unittest.main()
