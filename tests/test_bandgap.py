import tempfile
import unittest
from pathlib import Path

import numpy as np

from qeplotter.analysis.bandgap import _find_band_gap, detect_band_gap


class BandGapRegressionTests(unittest.TestCase):
    def _files(self, bands):
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        blocks = []
        x = np.arange(len(bands[0]), dtype=float)
        for band in bands:
            blocks.append("\n".join(
                f"{coordinate:.1f} {energy:.6f}"
                for coordinate, energy in zip(x, band)
            ))
        band_file = root / "bands.gnu"
        band_file.write_text("\n\n".join(blocks) + "\n")
        k_file = root / "K_POINTS"
        k_file.write_text(
            "K_POINTS crystal_b\n"
            "2\n"
            f"0 0 0 {len(x)} ! G\n"
            "0.5 0 0 1 ! X\n"
        )
        return directory, band_file, k_file

    def test_crossing_band_is_metallic_not_a_false_gap(self):
        directory, band_file, k_file = self._files([[-1.0, 1.0]])
        self.addCleanup(directory.cleanup)
        result = detect_band_gap(band_file, k_file, 0.0)
        self.assertTrue(result["metallic"])
        self.assertIsNone(
            _find_band_gap(
                np.array([0.0, 1.0]), np.array([[-1.0, 1.0]]),
                fermi_level=0.0,
            )
        )

    def test_direct_insulating_gap_is_reported(self):
        directory, band_file, k_file = self._files([
            [-1.0, -0.5, -1.0],
            [1.0, 0.5, 1.0],
        ])
        self.addCleanup(directory.cleanup)
        result = detect_band_gap(band_file, k_file, 0.0)
        self.assertFalse(result["metallic"])
        self.assertTrue(result["is_direct"])
        self.assertAlmostEqual(result["gap"], 1.0)

    def test_missing_fermi_is_explicitly_rejected(self):
        directory, band_file, k_file = self._files([
            [-1.0, -0.5],
            [0.5, 1.0],
        ])
        self.addCleanup(directory.cleanup)
        result = detect_band_gap(band_file, k_file)
        self.assertIsNone(result["metallic"])
        self.assertIn("Fermi level", result["error"])


if __name__ == "__main__":
    unittest.main()
