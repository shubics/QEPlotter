"""Computational Utilities page: converters, gap detector, bilayer analysis."""
import os
import shutil
from io import StringIO
from contextlib import redirect_stdout

import streamlit as st

from qeplotter.analysis.bandgap import detect_band_gap
from qeplotter.analysis.bilayer import analyse_file
from qeplotter.converters.fatbands import convert_consistent
from qeplotter.converters.soc import convert_soc_proj_to_ml
from gui.io_helpers import save_file, ensure_temp_dir
from gui.theme import remember_tab, remembered_tabs


_TOOLS_TAB_STATE = "utilities_active_tab"


def render_tools():
    st.title("Computational Utilities")
    st.caption("Converters and focused analysis helpers for Quantum ESPRESSO data.")

    tab1, tab2, tab3, tab4 = remembered_tabs(
        ["Standard Converter", "SOC Converter", "Gap Detector", "Bilayer Analysis"],
        _TOOLS_TAB_STATE,
    )

    # --- 1a. STANDARD CONVERTER ---
    with tab1:
        st.markdown("#### Standard `proj.out` → `.pdos` Converter")
        st.info("Converts standard non-SOC QE `projwfc.x` output to plotting-friendly format.")

        f_std = st.file_uploader(
            "Upload proj.out", key="t_p_std_uploader",
            on_change=remember_tab,
            args=(_TOOLS_TAB_STATE, "Standard Converter"),
        )

        if f_std:
            p = save_file(f_std)
            out_d = os.path.join(ensure_temp_dir(), "converted_pdos")

            if st.button(
                "Convert (Standard)", key="btn_conv_std", type="primary",
                on_click=remember_tab,
                args=(_TOOLS_TAB_STATE, "Standard Converter"),
            ):
                try:
                    if os.path.exists(out_d):
                        shutil.rmtree(out_d)
                    log = StringIO()
                    with redirect_stdout(log):
                        convert_consistent(p, outdir=out_d)
                    st.success("Conversion Complete!")
                    create_download_button(out_d, "converted_pdos.zip")
                    with st.expander("Logs"):
                        st.text(log.getvalue())
                except Exception as e:
                    st.error(f"Error: {e}")

    # --- 1b. SOC CONVERTER ---
    with tab2:
        st.markdown("#### SOC `proj.out` → `.pdos` Converter")
        st.info("Converts Spin-Orbit Coupled (SOC) QE `projwfc.x` output to plotting-friendly format.")

        f_soc = st.file_uploader(
            "Upload soc proj.out", key="t_p_soc_uploader",
            on_change=remember_tab,
            args=(_TOOLS_TAB_STATE, "SOC Converter"),
        )

        if f_soc:
            p = save_file(f_soc)
            out_d = os.path.join(ensure_temp_dir(), "soc_pdos")

            if st.button(
                "Convert (SOC Mode)", key="btn_conv_soc", type="primary",
                on_click=remember_tab,
                args=(_TOOLS_TAB_STATE, "SOC Converter"),
            ):
                try:
                    if os.path.exists(out_d):
                        shutil.rmtree(out_d)
                    log = StringIO()
                    with redirect_stdout(log):
                        convert_soc_proj_to_ml(p, outdir=out_d)
                    st.success("SOC Conversion Complete!")
                    create_download_button(out_d, "soc_pdos.zip")
                    with st.expander("Logs"):
                        st.text(log.getvalue())
                except Exception as e:
                    st.error(f"Error: {e}")

    # --- 2. GAP DETECTOR ---
    with tab3:
        st.markdown("#### Band Gap Analysis")
        st.info("Detect properties like direct/indirect gap, VBM, and CBM.")
        c1, c2 = st.columns(2)
        gap_change = {
            "on_change": remember_tab,
            "args": (_TOOLS_TAB_STATE, "Gap Detector"),
        }
        fb = c1.file_uploader("Band File (.gnu)", key="t_bg_uploader", **gap_change)
        fk = c2.file_uploader("K-Path File", key="t_kg_uploader", **gap_change)
        fermi = st.number_input("Fermi Level (eV)", value=0.0, format="%.4f", key="t_fermi_input", **gap_change)

        if st.button(
            "Analyze Gap", key="btn_analyze_gap", type="primary",
            on_click=remember_tab,
            args=(_TOOLS_TAB_STATE, "Gap Detector"),
        ) and fb and fk:
            run_tool(detect_band_gap, save_file(fb), save_file(fk), fermi)

    # --- 3. BILAYER ---
    with tab4:
        st.markdown("#### Structure Analyzer")
        st.info("Analyze layer separation and atomic coordinates.")
        fs = st.file_uploader(
            "Input File (.in / .out)", key="t_s_uploader",
            on_change=remember_tab,
            args=(_TOOLS_TAB_STATE, "Bilayer Analysis"),
        )
        if st.button(
            "Analyze Structure", key="btn_analyze_struc", type="primary",
            on_click=remember_tab,
            args=(_TOOLS_TAB_STATE, "Bilayer Analysis"),
        ) and fs:
            run_tool(analyse_file, save_file(fs))


def create_download_button(folder_path, zip_name):
    """Zips folder and creates download button."""
    if os.path.exists(folder_path):
        shutil.make_archive(folder_path, 'zip', folder_path)
        zip_file = folder_path + ".zip"
        with open(zip_file, "rb") as f:
            st.download_button(
                label="Download result (ZIP)",
                data=f,
                file_name=zip_name,
                mime="application/zip",
                type="primary"
            )


def run_tool(func, *args, **kwargs):
    log = StringIO()
    try:
        with redirect_stdout(log):
            func(*args, **kwargs)
        st.success("Execution Complete")
        with st.expander("View Logs", expanded=True):
            st.code(log.getvalue())
    except Exception as e:
        st.error(f"Error: {e}")
