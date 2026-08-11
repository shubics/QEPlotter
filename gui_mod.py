"""
QEPlotter — Streamlit dashboard (modular edition).

Thin router over the ``gui/`` package:
  • Structure      → CIF / POSCAR / QE 3D viewer + analysis   (gui.page_structure)
  • K-path         → native high-symmetry path + first BZ     (gui.page_kpath)
  • K-grid         → native uniform-grid IBZ reduction         (gui.page_kmesh)
  • Representations → Γ-point irreps, SALCs, orbital matching (gui.page_symmetry_representations)
  • Tools          → converters, gap detector, bilayer        (gui.page_tools)
  • Visualization  → band / fatbands / DOS / PDOS / overlay  (gui.page_plot)

Run with:  streamlit run gui_mod.py
"""
import os
import sys

import matplotlib

# Non-interactive backend (critical for headless / macOS) — set before pyplot import.
matplotlib.use("Agg")

import streamlit as st

# Make local packages importable regardless of the working directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gui.theme import configure_page, inject_css
from qeplotter.version import __version__

configure_page()
inject_css()

try:
    from gui.page_plot import render_dashboard
    from gui.page_structure import render_structure
    from gui.page_kpath import render_kpath
    from gui.page_kmesh import render_kmesh
    from gui.page_symmetry_representations import render_symmetry_representations
    from gui.page_tools import render_tools
except ImportError as _e:
    st.error(f"Critical error: could not import GUI modules. {_e}")
    st.stop()


def main():
    with st.sidebar:
        st.title("QEPlotter")
        st.caption(f"Scientific post-processing workspace · v{__version__}")
        st.divider()
        mode = st.radio(
            "Workspace",
            ["Crystal Structure", "K-path & Brillouin Zone",
             "K-grid & IBZ", "Symmetry Representations",
             "Utilities", "Plots & Data"],
        )

    if mode == "Crystal Structure":
        render_structure()
    elif mode == "K-path & Brillouin Zone":
        render_kpath()
    elif mode == "K-grid & IBZ":
        render_kmesh()
    elif mode == "Symmetry Representations":
        render_symmetry_representations()
    elif mode == "Utilities":
        render_tools()
    else:
        render_dashboard()


if __name__ == "__main__":
    main()
