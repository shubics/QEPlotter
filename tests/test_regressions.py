import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gui.io_helpers import save_file
from qeplotter.api import launch_gui
from qeplotter.converters.soc import _cg_prob
from qeplotter.core.io import read_band_xdistances
from qeplotter.plotting.dos import plot_pdos_dir
from qeplotter.plotting.fatbands import (
    _layer_material_labels,
    _layer_plot_title,
)


ROOT = Path(__file__).resolve().parents[1]


class GeneralRegressionTests(unittest.TestCase):
    def test_soc_clebsch_gordan_probabilities_need_no_sympy(self):
        self.assertAlmostEqual(_cg_prob(1, 1.5, 1.5, 1), 1.0)
        self.assertAlmostEqual(_cg_prob(1, 1.5, 0.5, 0), 2 / 3)
        self.assertAlmostEqual(_cg_prob(1, 1.5, 0.5, 1), 1 / 3)
        self.assertAlmostEqual(
            sum(_cg_prob(2, 2.5, 0.5, ml) for ml in range(-2, 3)), 1.0
        )

    def test_launch_gui_targets_modular_app_with_current_python(self):
        with patch("subprocess.run") as run:
            launch_gui()
        command = run.call_args.args[0]
        self.assertEqual(command[:3], [sys.executable, "-m", "streamlit"])
        self.assertEqual(Path(command[-1]), ROOT / "gui_mod.py")

    def test_setup_version_works_outside_repository(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, str(ROOT / "setup.py"), "--version"],
                cwd=directory,
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.stdout.strip(), "2.0")

    def test_final_high_symmetry_label_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            band_file = root / "bands.gnu"
            band_file.write_text("0 -1\n1 -0.5\n2 -0.2\n3 -0.1\n")
            k_file = root / "K_POINTS"
            k_file.write_text(
                "K_POINTS crystal_b\n"
                "3\n"
                "0 0 0 2 ! G\n"
                "0.5 0 0 2 ! X\n"
                "0.5 0.5 0 1 ! M\n"
            )
            _x, _bands, ticks, labels, _segments = read_band_xdistances(
                band_file, k_file
            )
        self.assertEqual(labels, ["G", "X", "M"])
        self.assertEqual(ticks[-1], 3.0)

    def test_pdos_draws_one_fermi_line_on_its_axis(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.pdos_atm#1(S)_wfc#1(s)"
            np.savetxt(path, np.array([[-1.0, 0.1], [0.0, 0.2], [1.0, 0.1]]))
            with patch("matplotlib.pyplot.show"):
                plot_pdos_dir(directory, fermi_level=0.0)
            axis = plt.gca()
            self.assertEqual(len(axis.lines), 2)  # one PDOS curve + one Fermi line
            plt.close(axis.figure)

    def test_layer_legend_uses_material_formulas(self):
        atom_names = ["W1", "Mo2", "S3", "S4", "S5", "S6"]
        assignment = {
            "W1": "bottom", "S4": "bottom", "S6": "bottom",
            "Mo2": "top", "S3": "top", "S5": "top",
        }
        self.assertEqual(
            _layer_material_labels(atom_names, assignment),
            ("WS₂", "MoS₂"),
        )

    def test_identical_layer_formulas_keep_lower_upper_identity(self):
        atom_names = ["Mo1", "S2", "S3", "Mo4", "S5", "S6"]
        assignment = {
            "Mo1": "bottom", "S2": "bottom", "S3": "bottom",
            "Mo4": "top", "S5": "top", "S6": "top",
        }
        self.assertEqual(
            _layer_material_labels(atom_names, assignment),
            ("MoS₂ · lower", "MoS₂ · upper"),
        )

    def test_layer_plot_title_names_hetero_and_homobilayers(self):
        self.assertEqual(
            _layer_plot_title("WS₂", "MoS₂"),
            "WS₂ / MoS₂ · layer-resolved fatbands",
        )
        self.assertEqual(
            _layer_plot_title("MoS₂ · lower", "MoS₂ · upper"),
            "MoS₂ bilayer · layer-resolved fatbands",
        )

    def test_uploaded_filename_cannot_escape_session_directory(self):
        uploaded = Mock()
        uploaded.name = "../../outside.dat"
        uploaded.getbuffer.return_value = memoryview(b"data")
        with tempfile.TemporaryDirectory() as directory:
            with patch("gui.io_helpers.ensure_temp_dir", return_value=directory):
                path = save_file(uploaded)
            self.assertEqual(Path(path), Path(directory) / "outside.dat")
            self.assertEqual(Path(path).read_bytes(), b"data")


if __name__ == "__main__":
    unittest.main()
