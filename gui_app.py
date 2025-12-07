import streamlit as st
import os
import tempfile
import matplotlib
import matplotlib.pyplot as plt
from io import StringIO, BytesIO
from contextlib import redirect_stdout
import ast
import re
import numpy as np

# --- 1. APP CONFIGURATION ---
st.set_page_config(
    page_title="QEPlotter Pro",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Force non-interactive backend (Critical for macOS)
matplotlib.use("Agg")

# Import Backend
try:
    import qep as qep5
except ImportError:
    st.error("❌ Critical Error: `qep.py` not found. Please put it in the same folder.")
    st.stop()
    


# --- 2. STYLING ---
st.markdown("""
<style>
    .stButton button { width: 100%; border-radius: 8px; font-weight: bold; padding: 0.6rem; }
    button[kind="primary"] { background: linear-gradient(90deg, #FF4B4B 0%, #D13333 100%); border: none; }
    .stExpander { border: 1px solid #e0e0e0; border-radius: 8px; }
    .stTextInput input, .stNumberInput input { border-radius: 5px; }
</style>
""", unsafe_allow_html=True)

# --- 3. HELPER FUNCTIONS & CACHING ---
if 'temp_dir' not in st.session_state:
    st.session_state.temp_dir = tempfile.mkdtemp()

def save_file(uploaded_file, subdir=None):
    if uploaded_file is None: return None
    target_dir = st.session_state.temp_dir
    if subdir:
        target_dir = os.path.join(target_dir, subdir)
        os.makedirs(target_dir, exist_ok=True)
    file_path = os.path.join(target_dir, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path

@st.cache_data
def get_fermi_from_scf(scf_path):
    """Parses scf.out to find the Fermi energy."""
    try:
        with open(scf_path, 'r', errors='ignore') as f:
            content = f.read()
        # Search for "Fermi energy is" or "highest occupied level"
        # Pattern: "the Fermi energy is     1.2345 eV"
        m = re.search(r"the Fermi energy is\s+([-+]?\d*\.\d+)\s+eV", content)
        if m:
            return float(m.group(1))
        
        # Pattern: "highest occupied level      1.2345 eV"
        m2 = re.search(r"highest occupied level\s+([-+]?\d*\.\d+)\s+eV", content)
        if m2:
            return float(m2.group(1))
            
        return None
    except Exception:
        return None

@st.cache_data
def generate_plot_image(args):
    """
    Cached wrapper for qep5.plot_from_file.
    Returns the BytesIO object of the image.
    args must be a hashable dictionary (or tuple of items).
    """
    # Convert args back to dict if needed, simpler to just run logic here but 
    # st.cache works best on pure data processing. 
    # Since qep5 plots to plt, we need to capture it.
    
    # NOTE: We can't easily cache the side-effect of plt.show/savefig inside qep5 
    # unless we rewrite qep5 to return Figure. 
    # Instead, we will cache the PARSING steps for interactive plots, 
    # but for static plots we might just let it run (Matplotlib is fast enough for single plots).
    pass 

# --- 4. MAIN LAYOUT ---
def main():
    with st.sidebar:
        st.title("⚛️ QEPlotter")
        st.caption("Full-Feature Interface v3.3 (Pro)")
        mode = st.radio("Navigation", ["📊 Visualization Dashboard", "🛠 Tools & Utilities"])

    if mode == "📊 Visualization Dashboard":
        render_dashboard()
    else:
        render_tools()

def render_dashboard():
    # --- STEP 1: PLOT TYPE SELECTION ---
    c_type, c_space = st.columns([2, 3])
    with c_type:
        plot_type_ui = st.selectbox("Select Plot Type",
                                    ["Band Structure", "Fatbands (Projected)", "Total DOS", "PDOS Only",
                                     "Overlay Comparison"])

    p_map = {
        "Band Structure": "band", "Fatbands (Projected)": "fatbands",
        "Total DOS": "dos", "PDOS Only": "pdos", "Overlay Comparison": "overlay_band"
    }
    pt = p_map[plot_type_ui]
    args = {'plot_type': pt, 'savefig': None}
    paths = {}

    col_inputs, col_preview = st.columns([1, 1.5])

    with col_inputs:
        # --- A. FILE INPUTS ---
        with st.expander("📂 1. Input Files", expanded=True):
            if pt in ["band", "fatbands", "overlay_band"]:
                f_band = st.file_uploader("Band File (.gnu)", type=["gnu", "dat"], key="u_band")
                paths['band_file'] = save_file(f_band)
                f_kpath = st.file_uploader("K-Path File", type=["kpath", "in", "txt"], key="u_kpath")
                paths['kpath_file'] = save_file(f_kpath)

            if pt in ["fatbands", "pdos"]:
                f_pdos = st.file_uploader("PDOS Files", accept_multiple_files=True, key="u_pdos")
                if f_pdos:
                    subdir = "pdos_data"
                    for f in f_pdos: save_file(f, subdir=subdir)
                    paths['fatband_dir'] = os.path.join(st.session_state.temp_dir, subdir)
                    paths['pdos_dir'] = paths['fatband_dir']
                else:
                    paths['fatband_dir'] = None

            f_dos = st.file_uploader("Total DOS File (Optional)", type=["dos", "dat"], key="u_dos")
            paths['dos_file'] = save_file(f_dos)
            
            if pt == "overlay_band":
                st.caption("Comparison Data")
                f_b2 = st.file_uploader("Band File 2", key="u_b2")
                paths['band_file2'] = save_file(f_b2)
                f_k2 = st.file_uploader("K-Path File 2", key="u_k2")
                paths['kpath_file2'] = save_file(f_k2)

        # --- B. DATA SETTINGS (Smart Automation) ---
        with st.expander("⚙️ 2. Data Settings", expanded=True):
            # Auto-Fermi
            st.markdown("##### Fermi Energy")
            c_f1, c_f2 = st.columns([1, 1])
            scf_file = c_f1.file_uploader("Auto-detect from `scf.out`", type=["out", "in"])
            
            detected_fermi = 0.0
            if scf_file:
                # Cache lookup
                scf_path = save_file(scf_file, subdir="scf")
                val = get_fermi_from_scf(scf_path)
                if val is not None:
                    detected_fermi = val
                    c_f1.success(f"Found E_F = {val} eV")
                else:
                    c_f1.warning("Not found in file")
            
            args['fermi_level'] = c_f2.number_input("Fermi Level (eV)", value=detected_fermi, format="%.4f")
            args['shift_fermi'] = c_f2.checkbox("Shift E_F to 0", value=True)

            # Limits
            st.markdown("##### Plot Limits")
            
            if pt in ["pdos", "dos"]:
                st.markdown("---")
                if st.checkbox("Set Custom Y-Limits", value=False):
                    c3, c4 = st.columns(2)
                    args['y_range'] = (c3.number_input("Y-Min", value=0.0, min_value=-100.0, max_value=100.0), 
                                       c4.number_input("Y-Max", value=10.0, min_value=-100.0, max_value=100.0))
                else:
                    args['y_range'] = None
            else:
                c3, c4 = st.columns(2)
                args['y_range'] = (c3.number_input("Y-Min", value=-3.0, min_value=-50.0, max_value=50.0), 
                                   c4.number_input("Y-Max", value=3.0, min_value=-50.0, max_value=50.0))

            if pt in ["band", "fatbands"]:
                args['spin'] = st.checkbox("Spin Polarized")
                args['sub_orb'] = st.checkbox("Sub-Orbital Analysis")

        # --- C. PLOT MODE ---
        with st.expander("🎨 3. Visualization Settings", expanded=True):
            # Matplotlib settings
            col_w, col_h = st.columns(2)
            fig_width = col_w.number_input("Width", 12)
            fig_height = col_h.number_input("Height", 6)
            args['dpi'] = st.number_input("DPI", 200)
            args['cmap_name'] = st.selectbox("Colormap", ["tab10", "magma", "viridis", "jet"])
            
            # Fatband specific (same as before)
            if pt == "fatbands":
                # ... existing fatband logic ...
                # (Simplifying for brevity, keeping core logic)
                fb_mode = st.selectbox("Mode", ["most", "atomic", "orbital", "element_orbital", "normal", "o_atomic", "o_orbital", "o_element_orbital", "heat_total", "heat_atomic", "heat_orbital", "heat_element_orbital", "layer"])
                args['fatbands_mode'] = fb_mode
                if "heat" in fb_mode or fb_mode in ['normal', 'most']:
                     args['highlight_channel'] = st.text_input("Highlight", "Mo")
                     args['overlay_bands_in_heat'] = st.checkbox("Overlay Lines", True)

            # Overlay specific
            if pt == "overlay_band":
                st.caption("Overlay Appearance")
                c_o1, c_o2 = st.columns(2)
                args['label1'] = c_o1.text_input("Label 1", "System A")
                args['color1'] = c_o1.color_picker("Color 1", "#FF0000")
                args['label2'] = c_o2.text_input("Label 2", "System B")
                args['color2'] = c_o2.color_picker("Color 2", "#0000FF")

    # ==========================================
    # RIGHT COLUMN: EXECUTION
    # ==========================================
    with col_preview:
        st.subheader("🚀 Execution")
        st.info("Interactive mode disabled. Adjust settings on the left and click Generate to update.")
        
        if st.button("GENERATE VISUALIZATION", type="primary"):
            args.update(paths)
            
            # Validation
            if pt == "band" and (not args.get('band_file') or not args.get('kpath_file')):
                st.error("Missing Band or K-Path file.")
                return

            with st.spinner("Processing..."):
                try:
                    # Static Path (Standard)
                    log_io = StringIO()
                    plt.close('all')
                    with redirect_stdout(log_io):
                        qep5.plot_from_file(**args)
                    
                    if plt.get_fignums():
                        fig = plt.gcf()
                        # Adjust size explicitly
                        fig.set_size_inches(fig_width, fig_height)
                        st.pyplot(fig, use_container_width=True)
                        
                        buf = BytesIO()
                        fig.savefig(buf, format="png", dpi=args.get('dpi', 200), bbox_inches='tight')
                        st.download_button("💾 Download PNG", buf, "plot.png", "image/png")
                    else:
                        st.warning("No plot generated.")

                except Exception as e:
                    st.error(f"Error: {e}")
                    import traceback
                    st.code(traceback.format_exc())

# ==========================================
# TOOLS PAGE (UTILITIES)
# ==========================================
import shutil 

def render_tools():
    st.title("🛠️ Computational Utilities")

    tab1, tab2, tab3 = st.tabs(["Fatband Converter", "Gap Detector", "Bilayer Analysis"])

    # --- 1. CONVERTER ---
    with tab1:
        st.markdown("#### `proj.out` → `.pdos` Converter")
        st.info("Converts QE `projwfc.x` output to plotting-friendly format.")

        f = st.file_uploader("Upload proj.out", key="t_p_uploader")

        if f:
            p = save_file(f)
            out_d = os.path.join(st.session_state.temp_dir, "converted_pdos")

            col1, col2 = st.columns(2)

            # STANDARD CONVERT
            if col1.button("Convert (Standard)", key="btn_conv_std", use_container_width=True):
                try:
                    # Clean previous run
                    if os.path.exists(out_d): shutil.rmtree(out_d)
                    
                    log = StringIO()
                    with redirect_stdout(log):
                        qep5.convert_consistent(p, outdir=out_d)
                    st.success("Conversion Complete!")
                    create_download_button(out_d, "converted_pdos.zip")
                    with st.expander("Logs"): st.text(log.getvalue())
                except Exception as e:
                    st.error(f"Error: {e}")

            # SOC CONVERT
            if col2.button("Convert (SOC Mode)", key="btn_conv_soc", use_container_width=True):
                try:
                    if os.path.exists(out_d): shutil.rmtree(out_d)
                    
                    log = StringIO()
                    with redirect_stdout(log):
                        qep5.convert_soc_proj_to_ml(p, outdir=out_d)
                    st.success("SOC Conversion Complete!")
                    create_download_button(out_d, "soc_pdos.zip")
                    with st.expander("Logs"): st.text(log.getvalue())
                except Exception as e:
                    st.error(f"Error: {e}")

    # --- 2. GAP DETECTOR ---
    with tab2:
        st.markdown("#### Band Gap Analysis")
        c1, c2 = st.columns(2)
        fb = c1.file_uploader("Band File", key="t_bg_uploader")
        fk = c2.file_uploader("K-Path File", key="t_kg_uploader")
        fermi = st.number_input("Fermi Level", key="t_fermi_input")

        if st.button("Analyze Gap", key="btn_analyze_gap", type="primary") and fb and fk:
            run_tool(qep5.detect_band_gap, save_file(fb), save_file(fk), fermi)

    # --- 3. BILAYER ---
    with tab3:
        st.markdown("#### Structure Analyzer")
        fs = st.file_uploader("Input File", key="t_s_uploader")
        if st.button("Analyze Structure", key="btn_analyze_struc", type="primary") and fs:
            run_tool(qep5.analyse_file, save_file(fs))

def create_download_button(folder_path, zip_name):
    """Zips folder and creates download button."""
    if os.path.exists(folder_path):
        shutil.make_archive(folder_path, 'zip', folder_path)
        zip_file = folder_path + ".zip"
        with open(zip_file, "rb") as f:
            st.download_button(
                label="⬇️ Download Result (ZIP)",
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
            
if __name__ == "__main__":
    main()
