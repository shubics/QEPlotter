import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qeplotter.api import plot_from_file
from qeplotter.plotting.dos import plot_dos
from qeplotter.plotting.fatbands import plot_fatbands
from qeplotter.plotting.style import figure_size


class PlotControlTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_dos_applies_common_text_layout_and_legend_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            dos_file = Path(directory) / "sample.dos"
            np.savetxt(
                dos_file,
                np.array([[-1.0, 0.1], [0.0, 0.4], [1.0, 0.2]]),
            )
            with patch("matplotlib.pyplot.show"):
                plot_dos(
                    dos_file,
                    figsize=(7.5, 4.5),
                    plot_title="Custom DOS",
                    x_label="Custom energy",
                    y_label="Custom density",
                    show_grid=False,
                    show_legend=False,
                )

        fig = plt.gcf()
        ax = fig.axes[0]
        np.testing.assert_allclose(fig.get_size_inches(), [7.5, 4.5])
        self.assertEqual(ax.get_title(), "Custom DOS")
        self.assertEqual(ax.get_xlabel(), "Custom energy")
        self.assertEqual(ax.get_ylabel(), "Custom density")
        self.assertIsNone(ax.get_legend())
        self.assertFalse(any(line.get_visible() for line in ax.get_xgridlines()))

    def test_empty_custom_text_preserves_automatic_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            dos_file = Path(directory) / "sample.dos"
            np.savetxt(dos_file, np.array([[-1.0, 0.1], [1.0, 0.2]]))
            with patch("matplotlib.pyplot.show"):
                plot_dos(
                    dos_file,
                    plot_title="  ",
                    x_label="",
                    y_label=None,
                    legend_title="DOS curves",
                    legend_location="lower left",
                )

        ax = plt.gca()
        self.assertEqual(ax.get_title(), "Total DOS")
        self.assertEqual(ax.get_xlabel(), "Energy (eV)")
        self.assertEqual(ax.get_ylabel(), "DOS")
        self.assertEqual(ax.get_legend().get_title().get_text(), "DOS curves")

    def test_high_level_api_forwards_common_style_options(self):
        options = {
            "figsize": (9, 5),
            "plot_title": "API title",
            "x_label": "kx",
            "y_label": "Energy",
            "show_title": False,
            "show_grid": False,
            "show_legend": False,
            "legend_location": "upper left",
            "legend_title": "Channels",
        }
        with patch("qeplotter.api.plot_dos") as mocked_plot_dos:
            plot_from_file(plot_type="dos", dos_file="sample.dos", **options)

        forwarded = mocked_plot_dos.call_args.kwargs
        for key, value in options.items():
            self.assertEqual(forwarded[key], value)

    @patch("qeplotter.plotting.fatbands.read_band_xdistances")
    @patch("qeplotter.plotting.fatbands.read_fatband_files")
    def test_fatband_legend_control_manages_continuous_colour_scale(
        self, mocked_read_fatbands, mocked_read_bands
    ):
        mocked_read_fatbands.return_value = (
            [("Mo1", "d")],
            np.array([1, 2]),
            np.array([[-1.0, 0.0], [-0.8, 0.2]]),
            [np.array([[0.2, 0.4], [0.3, 0.5]])],
        )
        mocked_read_bands.return_value = (
            np.array([0.0, 1.0]),
            np.array([[-1.0, -0.8]]),
            [0.0, 1.0],
            ["G", "X"],
            [(0, 1)],
        )

        with patch("matplotlib.pyplot.show"):
            plot_fatbands(
                "projections",
                "kpath.in",
                "bands.gnu",
                mode="heat_total",
                show_legend=False,
            )
        self.assertEqual(len(plt.gcf().axes), 1)

        plt.close("all")
        with patch("matplotlib.pyplot.show"):
            plot_fatbands(
                "projections",
                "kpath.in",
                "bands.gnu",
                mode="heat_total",
                show_legend=True,
                legend_title="Projection weight",
            )
        self.assertEqual(len(plt.gcf().axes), 2)
        self.assertEqual(plt.gcf().axes[1].get_ylabel(), "Projection weight")

    def test_figure_size_rejects_invalid_values(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            figure_size((8, 0), (6, 6))
        with self.assertRaisesRegex(ValueError, "exactly two"):
            figure_size((8,), (6, 6))


if __name__ == "__main__":
    unittest.main()
