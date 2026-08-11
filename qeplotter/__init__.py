"""
QEPlotter – Quantum ESPRESSO plotting toolkit.

This package is the maintained modular implementation evolved from the archived
``archive/monolithic/qep.py`` script:
  - core/utils.py    : shared constants and helpers
  - core/io.py       : file parsers (kpath, bands, fatbands)
  - plotting/bands.py    : band structure plotting
  - plotting/dos.py      : DOS and PDOS plotting
  - plotting/fatbands.py : fatband plotting (bubble, line, heat modes)
  - analysis/bandgap.py  : band gap detection
  - analysis/bilayer.py  : bilayer stacking analysis
  - converters/fatbands.py : projwfc output standardizer
  - converters/soc.py     : SOC to (l, ml) basis converter
  - api.py               : high-level plot_from_file dispatcher
"""

from qeplotter.api import plot_from_file, launch_gui
from qeplotter.analysis.bandgap import detect_band_gap
from qeplotter.analysis.bilayer import (
    analyse_file,
    analyse_stacking,
    classify_stacking,
    detect_bilayer,
)
from qeplotter.converters.fatbands import convert_consistent
from qeplotter.converters.soc import convert_soc_proj_to_ml
from qeplotter.version import __version__
