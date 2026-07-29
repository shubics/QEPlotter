"""
Streamlit GUI package for QEPlotter.

Split out of the former monolithic ``gui_mod.py`` into focused modules:

  - theme.py          : page config + dark-theme CSS
  - io_helpers.py     : upload/temp-file helpers, Fermi auto-detect, channel scan
  - page_plot.py      : Visualization dashboard (band / fatbands / DOS / overlay)
  - page_structure.py : Crystal-structure viewer & analysis (CIF / POSCAR / QE)
  - page_tools.py     : Converters, gap detector, bilayer analysis

``gui_mod.py`` at the repo root is now a thin router over these pages.
"""
