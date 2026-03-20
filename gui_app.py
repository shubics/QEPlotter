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
    page_title="QEPlotter",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- CUSTOM CSS ---
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap');
    
    html, body, [class*="css"]  {
        font-family: 'Inter', sans-serif;
    }
    
    /* Sleek Header Gradient */
    h1 {
        background: -webkit-linear-gradient(45deg, #3B82F6, #8B5CF6);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-weight: 700 !important;
    }
    
    /* Modern Containers */
    .stApp {
        background-color: #0F172A;
    }
    
    /* Button Styling */
    .stButton>button {
        border-radius: 8px;
        font-weight: 600;
        border: none;
        transition: all 0.2s ease;
    }
    
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
    }
</style>
""", unsafe_allow_html=True)

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
        
        # Pattern 1: Metal - "the Fermi energy is     1.2345 eV"
        m = re.search(r"the Fermi energy is\s+([-+]?\d*\.?\d+)\s+eV", content)
        if m:
            return float(m.group(1))
        
        # Pattern 2: Insulator/Semiconductor - "highest occupied, lowest unoccupied level (ev):    -4.7993   -3.6536"
        m2 = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)", content)
        if m2:
            homo = float(m2.group(1))
            lumo = float(m2.group(2))
            return 0.5 * (homo + lumo)  # mid-gap as Fermi estimate
        
        # Pattern 3: "highest occupied level (ev):     1.2345"
        m3 = re.search(r"highest occupied level\s*\(ev\):\s+([-+]?\d*\.?\d+)", content)
        if m3:
            return float(m3.group(1))
            
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

def get_available_channels(pdos_dir):
    import glob, re, os
    if not pdos_dir or not os.path.exists(pdos_dir): return [], [], [], []
    files = glob.glob(os.path.join(pdos_dir, '*pdos*'))
    
    atoms = set()
    elements = set()
    orbitals = set()
    elem_orbs = set()
    
    pattern = re.compile(r'atm#(\d+)\(([A-Za-z]+)\)_wfc#\d+\(([a-zA-Z0-9_.]+)\)')
    for f in files:
        m = pattern.search(os.path.basename(f))
        if m:
            num, elem, orb = m.groups()
            base_orb = orb.split('_')[0] 
            
            atoms.add(f"{elem}{num}")
            elements.add(elem)
            orbitals.add(base_orb)
            elem_orbs.add(f"{elem}-{base_orb}")
    
    return sorted(list(atoms)), sorted(list(elements)), sorted(list(orbitals)), sorted(list(elem_orbs))

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
        tab_data, tab_settings, tab_style = st.tabs(["📂 Data & Files", "⚙️ Core Settings", "🎨 Plot Styling"])
        
        # --- A. FILE INPUTS ---
        with tab_data:
            if pt in ["band", "fatbands", "overlay_band"]:
                st.markdown("##### Band Structure Data")
                f_band = st.file_uploader("Band File (.gnu)", type=["gnu", "dat"], key="u_band", help="Quantum ESPRESSO bands.dat.gnu file")
                paths['band_file'] = save_file(f_band)
                f_kpath = st.file_uploader("K-Path File", type=["kpath", "in", "txt"], key="u_kpath", help="K_POINTS file in crystal_b format for symmetry labels")
                paths['kpath_file'] = save_file(f_kpath)

            if pt in ["fatbands", "pdos", "band"]:
                st.markdown("##### PDOS Projection Data")
                if pt == "band":
                    st.caption("Required for colored band modes (atomic, orbital, most, etc.)")
                f_pdos = st.file_uploader("PDOS Files", accept_multiple_files=True, key="u_pdos", help="Select all outdir/*pdos* files produced by projwfc.x")
                if f_pdos:
                    subdir = "pdos_data"
                    for f in f_pdos: save_file(f, subdir=subdir)
                    paths['fatband_dir'] = os.path.join(st.session_state.temp_dir, subdir)
                    paths['pdos_dir'] = paths['fatband_dir']
                else:
                    paths['fatband_dir'] = None

            st.markdown("##### DOS Data")
            f_dos = st.file_uploader("Total DOS File (Optional)", key="u_dos", help="File containing Total Density of States")
            paths['dos_file'] = save_file(f_dos)
            
            if pt == "overlay_band":
                st.markdown("##### Comparison Data")
                f_b2 = st.file_uploader("Band File 2", key="u_b2")
                paths['band_file2'] = save_file(f_b2)
                f_k2 = st.file_uploader("K-Path File 2", key="u_k2")
                paths['kpath_file2'] = save_file(f_k2)

        # --- B. DATA SETTINGS ---
        with tab_settings:
            # Auto-Fermi
            st.markdown("##### Fermi Energy")
            f_scf = st.file_uploader("SCF Output (auto-detect Fermi)", type=["out", "txt"], key="u_scf", help="Upload scf.out to auto-read Fermi energy")
            auto_fermi = 0.0
            if f_scf:
                scf_path = save_file(f_scf)
                paths['scf_file'] = scf_path
                detected = get_fermi_from_scf(scf_path)
                if detected is not None:
                    auto_fermi = detected
                    st.success(f"Fermi energy detected: **{detected:.4f} eV**")
                else:
                    st.warning("Could not find Fermi energy in this file.")
            c_f1, c_f2 = st.columns(2)
            args['fermi_level'] = c_f1.number_input("Fermi Level (eV)", value=auto_fermi, format="%.4f", help="Absolute Fermi energy to shift plots relative to")
            args['shift_fermi'] = c_f2.checkbox("Shift E_F to 0", value=True, help="Shift energy axis so Fermi level is at 0")

            if pt in ["band", "fatbands", "pdos", "dos"]:
                st.markdown("##### Calculation Properties")
                c_prop1, c_prop2 = st.columns(2)
                args['spin'] = c_prop1.checkbox("Spin Polarized", help="Check if calculation used nspin=2 or noncolin=true")
                args['sub_orb'] = c_prop2.checkbox("Sub-Orbital Analysis", help="Check if you want m-resolved or SOC states")

            if pt == "band":
                 bm = st.selectbox("Band Mode", ["normal", "atomic", "orbital", "element_orbital", "most"], help="Mode for coloring bands")
                 args['band_mode'] = bm
                 if bm != 'normal':
                     st.info("⚠️ For colored bands, you must upload Fatband/PDOS files in the Data tab.")
                 else:
                     args['band_mode'] = 'normal'

            if pt == "fatbands":
                 st.markdown("##### Fatband Projection")
                 
                 fb_style = st.selectbox("Plot Style", ["Scatter (Bubble)", "Lines (o_)", "Heatmap (heat_)"], help="Visual style for projecting orbital weights")
                 
                 proj_opts = []
                 if fb_style == "Scatter (Bubble)":
                     proj_opts = ["Most Dominant", "Atomic", "Orbital", "Element-Orbital"]
                 elif fb_style == "Lines (o_)":
                     proj_opts = ["Layer", "Normal", "Atomic", "Orbital", "Element-Orbital"]
                 elif fb_style == "Heatmap (heat_)":
                     proj_opts = ["Total", "Atomic", "Orbital", "Element-Orbital"]
                 
                 fb_proj = st.selectbox("Projection Type", proj_opts, help="How to group projected orbitals")
                 
                 mode_map = {
                     "Most Dominant": "most",
                     "Normal": "normal",
                     "Atomic": "atomic", 
                     "Orbital": "orbital",
                     "Element-Orbital": "element_orbital",
                     "Layer": "layer",
                     "Total": "total"
                 }
                 
                 base_m = mode_map[fb_proj]
                 
                 if fb_style == "Lines (o_)":
                     if base_m in ["normal", "layer"]:
                         fb_mode = base_m
                     else:
                         fb_mode = f"o_{base_m}"
                 elif fb_style == "Heatmap (heat_)":
                     fb_mode = f"heat_{base_m}"
                 else:
                     fb_mode = base_m
                 
                 args['fatbands_mode'] = fb_mode
                 st.caption(f"*(Internal mode: `{fb_mode}`)*")

                 # Layer Assignment Logic
                 if base_m == 'layer':
                     st.caption("Layer Mapping (Assign 'top' or 'bottom' to atoms)")
                     a_list, _, _, _ = get_available_channels(paths.get('fatband_dir'))
                     
                     st.markdown("**(Optional) Auto-Assign from Structure**")
                     f_struc = st.file_uploader("Upload .in or .out file", key="u_struc_layer", help="Upload a QE structure file to automatically detect layers based on Z-coordinates")
                     
                     auto_top, auto_bot = [], []
                     if f_struc:
                         try:
                             from qep import parse_qe_block, custom_labeling
                             import numpy as np
                             content = f_struc.getvalue().decode('utf-8').splitlines()
                             cell, species, frac = parse_qe_block(content)
                             if len(species) > 0:
                                 # Robust PBC-aware median split for slabs
                                 z = np.mod(frac[:, 2], 1.0)
                                 order = np.argsort(z)
                                 sorted_z = z[order]
                                 gaps = np.append(np.diff(sorted_z), 1.0 + sorted_z[0] - sorted_z[-1])
                                 
                                 max_gap_idx = np.argmax(gaps)
                                 shift = sorted_z[(max_gap_idx + 1) % len(z)]
                                 shifted_z = np.mod(z - shift, 1.0)
                                 median_z = np.median(shifted_z)
                                 
                                 # QE projwfc.x outputs global 1-based indexing for atoms (e.g. Se3, Se4)
                                 labels = [f"{sp}{i+1}" for i, sp in enumerate(species)]
                                 for i, sz in enumerate(shifted_z):
                                     if labels[i] in a_list:
                                         if sz > median_z:
                                             auto_top.append(labels[i])
                                         else:
                                             auto_bot.append(labels[i])
                                 st.success("Layers correctly detected via median splitting!")
                         except Exception as e:
                             st.warning(f"Could not auto-detect layers: {e}")
                             
                     if not a_list:
                         st.warning("Please upload PDOS files in the Data tab to enable layer mapping.")
                     else:
                         top_atoms = st.multiselect("Top Layer Atoms", a_list, default=auto_top)
                         bot_bot_options = [x for x in a_list if x not in top_atoms]
                         valid_auto_bot = [x for x in auto_bot if x in bot_bot_options]
                         bot_atoms = st.multiselect("Bottom Layer Atoms", bot_bot_options, default=valid_auto_bot)
                         
                         l_map = {}
                         for a in top_atoms: l_map[a] = 'top'
                         for a in bot_atoms: l_map[a] = 'bottom'
                         if l_map:
                             args['layer_assignment'] = l_map
                             
                 # Dynamic Highlight Channels List
                 atoms, elements, orbitals, exp_orbs = get_available_channels(paths.get('fatband_dir'))
                 hl_options = elements # default
                 if base_m == 'orbital':
                     hl_options = orbitals
                 elif base_m == 'element_orbital':
                     hl_options = exp_orbs
                 elif base_m == 'atomic':
                     hl_options = elements
                     
                 # Inject generic if lists are empty (e.g. before upload)
                 if not hl_options: hl_options = ["Mo", "S", "d", "p", "Mo-d"]

                 if fb_mode in ['o_orbital', 'o_atomic', 'o_element_orbital']:
                      args['dual'] = st.checkbox("Dual Channel Mode", help="Highlight two contrasting channels with a diverging colormap")
                 
                 if args.get('dual'):
                      c_h1, c_h2 = st.columns(2)
                      idx2 = 1 if len(hl_options) > 1 else 0
                      h1 = c_h1.selectbox("Channel 1", hl_options, index=0)
                      h2 = c_h2.selectbox("Channel 2", hl_options, index=idx2)
                      args['highlight_channel'] = (h1, h2)
                 elif "heat" in fb_mode or fb_mode in ['normal', 'most', 'o_orbital', 'o_atomic', 'o_element_orbital']:
                      args['highlight_channel'] = st.selectbox("Highlight Channel", hl_options, index=0, help="Specific element/orbital to highlight")
                 
                 if "heat" in fb_mode:
                     args['overlay_bands_in_heat'] = st.checkbox("Overlay Lines", True, help="Add line bands on top of heatmap")

            args['plot_total_dos'] = st.checkbox("Plot Total DOS side-by-side", value=False, help="Requires DOS file uploaded")

            if pt in ["band", "fatbands"]:
                args['show_band_gap'] = st.checkbox("📏 Show Band Gap Arrow", value=False, help="Detect and annotate the band gap (VBM → CBM) on the plot")

        # --- C. PLOT MODE & STYLE ---
        with tab_style:
            st.markdown("##### Plot Dimensions & Limits")
            col_w, col_h = st.columns(2)
            fig_width = col_w.number_input("Width", 12)
            fig_height = col_h.number_input("Height", 6)
            args['dpi'] = st.number_input("DPI", 200)

            c3, c4 = st.columns(2)
            if pt in ["pdos", "dos"]:
                if st.checkbox("Set Custom Y-Limits", value=False):
                    args['y_range'] = (c3.number_input("Y-Min", value=0.0, min_value=-100.0, max_value=100.0), 
                                       c4.number_input("Y-Max", value=10.0, min_value=-100.0, max_value=100.0))
                else:
                    args['y_range'] = None
            else:
                args['y_range'] = (c3.number_input("Y-Min", value=-3.0, min_value=-50.0, max_value=50.0), 
                                   c4.number_input("Y-Max", value=3.0, min_value=-50.0, max_value=50.0))
            
            if pt == 'dos' or args.get('plot_total_dos', False):
                 c5, c6 = st.columns(2)
                 use_x = c5.checkbox("Set Custom DOS Limits (X-Axis)", value=False)
                 if use_x:
                     args['x_range'] = (0.0, c6.number_input("Max DOS Value", value=10.0))
                 elif pt == 'dos':
                     use_x2 = st.checkbox("Set Custom Energy Limits", value=False)
                     if use_x2:
                        args['x_range'] = (c5.number_input("Energy-Min", value=-10.0), 
                                           c6.number_input("Energy-Max", value=10.0))
                     else:
                        args['x_range'] = None
                 else:
                     args['x_range'] = None

            if pt in ['dos', 'pdos']:
                args['vertical'] = st.checkbox("Vertical Orientation (Energy on Y)", value=True)
            
            st.markdown("##### Colors & Visuals")
            args['cmap_name'] = st.selectbox("Colormap", ["tab10", "magma", "viridis", "jet", "coolwarm", "bwr"], help="Matplotlib colormap")
            
            if pt == "fatbands":
                c_adv1, c_adv2 = st.columns(2)
                args['s_min'] = c_adv1.number_input("Min Bubble Size", 1.0, 100.0, 10.0)
                args['s_max'] = c_adv2.number_input("Max Bubble Size", 10.0, 500.0, 100.0)
                args['weight_threshold'] = c_adv1.number_input("Weight Threshold", 0.0, 1.0, 0.01)
                
                if "heat_" in args.get('fatbands_mode', ''):
                    args['heat_vmin'] = c_adv1.number_input("Heatmap Min Value", value=0.0)
                    args['heat_vmax'] = c_adv2.number_input("Heatmap Max Value", value=0.0)
                    if args['heat_vmax'] == 0.0: args['heat_vmax'] = None

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
        
        if st.button("GENERATE VISUALIZATION", type="primary"):
            args.update(paths)
            
            # Validation
            if pt == "band" and (not args.get('band_file') or not args.get('kpath_file')):
                st.error("Missing Band or K-Path file.")
                return
            if pt == "band" and args.get('band_mode', 'normal') != 'normal' and not args.get('fatband_dir'):
                st.error(f"Band mode '{args['band_mode']}' requires PDOS projection files. Upload them in the Data tab.")
                return
            if pt == "fatbands" and not args.get('fatband_dir'):
                st.error("Fatband mode requires PDOS projection files. Upload them in the Data tab.")
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
                        
                        buf = BytesIO()
                        fig.savefig(buf, format="png", dpi=args.get('dpi', 200), bbox_inches='tight')
                        st.image(buf.getvalue(), caption="Plot Preview", use_container_width=True)
                        
                        st.download_button("💾 Download PNG", buf, "plot.png", "image/png")
                    else:
                        st.warning("No plot generated.")

                    # Show backend logs
                    log_text = log_io.getvalue()
                    if log_text.strip():
                        with st.expander("📋 Backend Logs", expanded=False):
                            st.code(log_text)

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
        st.info("Detect properties like direct/indirect gap, VBM, and CBM.")
        c1, c2 = st.columns(2)
        fb = c1.file_uploader("Band File (.gnu)", key="t_bg_uploader")
        fk = c2.file_uploader("K-Path File", key="t_kg_uploader")
        fermi = st.number_input("Fermi Level (eV)", value=0.0, format="%.4f", key="t_fermi_input")

        if st.button("Analyze Gap", key="btn_analyze_gap", type="primary") and fb and fk:
            run_tool(qep5.detect_band_gap, save_file(fb), save_file(fk), fermi)

    # --- 3. BILAYER ---
    with tab3:
        st.markdown("#### Structure Analyzer")
        st.info("Analyze layer separation and atomic coordinates.")
        fs = st.file_uploader("Input File (.in / .out)", key="t_s_uploader")
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
