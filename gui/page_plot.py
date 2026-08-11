"""Visualization Dashboard page: band / fatbands / DOS / PDOS / overlay plots."""
import os
from io import StringIO, BytesIO
from contextlib import redirect_stdout

import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex

from qeplotter.api import plot_from_file
from gui.io_helpers import save_file, get_fermi_from_scf, get_available_channels
from gui.theme import remember_tab, remembered_tabs
from qeplotter.plotting.fatbands import _layer_material_labels


_COLORMAP_GROUPS = {
    "Categorical": (
        "tab10", "tab20", "Set1", "Set2", "Set3", "Paired", "Accent",
        "Dark2",
    ),
    "Perceptual sequential": (
        "viridis", "plasma", "inferno", "magma", "cividis", "turbo",
        "cubehelix",
    ),
    "Sequential": (
        "Blues", "Greens", "Oranges", "Purples", "Greys", "YlGnBu",
        "YlOrRd", "PuBuGn",
    ),
    "Diverging": (
        "coolwarm", "RdBu_r", "Spectral", "PiYG", "PRGn", "BrBG", "bwr",
        "seismic",
    ),
    "Cyclic": ("twilight", "twilight_shifted", "hsv"),
    "Legacy": ("jet",),
}
_COLORMAP_OPTIONS = tuple(
    name for names in _COLORMAP_GROUPS.values() for name in names
)
_CATEGORICAL_COLORMAPS = frozenset(_COLORMAP_GROUPS["Categorical"])
_PLOT_TAB_STATE = "plot_configuration_tab"


def _colormap_group(cmap_name):
    for group, names in _COLORMAP_GROUPS.items():
        if cmap_name in names:
            return group
    raise ValueError(f"Unknown GUI colormap: {cmap_name}")


def _colormap_option_label(cmap_name):
    return f"{_colormap_group(cmap_name)} · {cmap_name}"


def _colormap_preview_html(cmap_name):
    """Return an accessible HTML preview of a Matplotlib colormap."""
    cmap = plt.get_cmap(cmap_name)
    categorical = cmap_name in _CATEGORICAL_COLORMAPS

    if categorical:
        count = min(int(getattr(cmap, "N", 10)), 20)
        colours = [to_hex(cmap(index / max(count - 1, 1))) for index in range(count)]
        stops = []
        for index, colour in enumerate(colours):
            left = 100 * index / count
            right = 100 * (index + 1) / count
            stops.extend((f"{colour} {left:.3f}%", f"{colour} {right:.3f}%"))
        scale_text = f"{count} discrete colours"
    else:
        count = 17
        colours = [to_hex(cmap(index / (count - 1))) for index in range(count)]
        stops = [
            f"{colour} {100 * index / (count - 1):.3f}%"
            for index, colour in enumerate(colours)
        ]
        scale_text = "low → high"

    gradient = ", ".join(stops)
    return f"""
    <div aria-label="{cmap_name} colormap preview" style="margin:.25rem 0 .9rem;">
      <div style="height:30px;border:1px solid #505965;border-radius:4px;
                  background:linear-gradient(90deg,{gradient});"></div>
      <div style="display:flex;justify-content:space-between;margin-top:4px;
                  color:#A8B0BA;font-size:.76rem;">
        <span>{cmap_name}</span><span>{scale_text}</span>
      </div>
    </div>
    """


def render_dashboard():
    st.title("Plots & Data")
    st.caption(
        "Prepare band, fatband, DOS, and comparison figures from Quantum "
        "ESPRESSO outputs."
    )
    # --- STEP 1: PLOT TYPE SELECTION ---
    c_type, c_space = st.columns([2, 3])
    with c_type:
        plot_type_ui = st.selectbox(
            "Select Plot Type",
            ["Band Structure", "Fatbands (Projected)", "Total DOS", "PDOS Only",
             "Overlay Comparison"],
            key="plot_type_ui",
            on_change=remember_tab,
            args=(_PLOT_TAB_STATE, "Data & Files"),
        )

    p_map = {
        "Band Structure": "band", "Fatbands (Projected)": "fatbands",
        "Total DOS": "dos", "PDOS Only": "pdos", "Overlay Comparison": "overlay_band"
    }
    pt = p_map[plot_type_ui]
    args = {'plot_type': pt, 'savefig': None}
    paths = {}

    col_inputs, col_preview = st.columns([1, 1.5])

    with col_inputs:
        tab_data, tab_settings, tab_style = remembered_tabs(
            ["Data & Files", "Core Settings", "Plot Styling"],
            _PLOT_TAB_STATE,
        )
        settings_change = {
            "on_change": remember_tab,
            "args": (_PLOT_TAB_STATE, "Core Settings"),
        }
        style_change = {
            "on_change": remember_tab,
            "args": (_PLOT_TAB_STATE, "Plot Styling"),
        }

        # --- A. FILE INPUTS ---
        with tab_data:
            if pt in ["band", "fatbands", "overlay_band"]:
                st.markdown("##### Band Structure Data")
                f_band = st.file_uploader("Band File (.gnu)", type=["gnu", "dat"], key="u_band", help="Quantum ESPRESSO bands.dat.gnu file")
                paths['band_file'] = save_file(f_band)
                f_kpath = st.file_uploader("K-Path File", type=["kpath", "in", "txt"], key="u_kpath", help="K_POINTS file in crystal_b format for symmetry labels")
                paths['kpath_file'] = save_file(f_kpath)

            if pt == "dos":
                st.markdown("##### DOS Data (Required)")
                f_dos = st.file_uploader("Total DOS File", key="u_dos", help="File from dos.x (e.g. system.dos)")
                paths['dos_file'] = save_file(f_dos)

            if pt in ["fatbands", "pdos", "band"]:
                st.markdown("##### PDOS Projection Data")
                if pt == "band":
                    st.caption("Required for colored band modes (atomic, orbital, most, etc.)")
                elif pt == "pdos":
                    st.caption("Upload all *pdos* files from projwfc.x output")
                f_pdos = st.file_uploader("PDOS Files", accept_multiple_files=True, key="u_pdos", help="Select all outdir/*pdos* files produced by projwfc.x")
                if f_pdos:
                    subdir = "pdos_data"
                    for f in f_pdos:
                        save_file(f, subdir=subdir)
                    paths['fatband_dir'] = os.path.join(st.session_state.temp_dir, subdir)
                    paths['pdos_dir'] = paths['fatband_dir']
                else:
                    paths['fatband_dir'] = None

            if pt != "dos":
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
            f_scf = st.file_uploader("SCF Output (auto-detect Fermi)", type=["out", "txt"], key="u_scf", help="Upload scf.out to auto-read Fermi energy", **settings_change)
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
            args['fermi_level'] = c_f1.number_input("Fermi Level (eV)", value=auto_fermi, format="%.4f", help="Absolute Fermi energy to shift plots relative to", **settings_change)
            args['shift_fermi'] = c_f2.checkbox("Shift E_F to 0", value=True, help="Shift energy axis so Fermi level is at 0", **settings_change)

            if pt in ["band", "fatbands"]:
                st.markdown("##### Calculation Properties")
                c_prop1, c_prop2 = st.columns(2)
                args['spin'] = c_prop1.checkbox("Spin Polarized", help="Check if calculation used nspin=2 or noncolin=true", **settings_change)
                args['sub_orb'] = c_prop2.checkbox("Sub-Orbital Analysis", help="Check if you want m-resolved or SOC states", **settings_change)

            if pt == "pdos":
                st.markdown("##### PDOS Settings")
                args['pdos_mode'] = st.selectbox("PDOS Grouping Mode", ["atomic", "orbital", "element_orbital"],
                    help="How to group projected orbitals: by atom element, orbital type, or element-orbital pair", **settings_change)

            if pt == "band":
                bm = st.selectbox("Band Mode", ["normal", "atomic", "orbital", "element_orbital", "most"], help="Mode for coloring bands", **settings_change)
                args['band_mode'] = bm
                if bm != 'normal':
                    st.info("Colored bands require Fatband/PDOS files in the Data tab.")
                else:
                    args['band_mode'] = 'normal'

            if pt == "fatbands":
                st.markdown("##### Fatband Projection")

                fb_style = st.selectbox("Plot Style", ["Scatter (Bubble)", "Lines (o_)", "Heatmap (heat_)"], help="Visual style for projecting orbital weights", **settings_change)

                proj_opts = []
                if fb_style == "Scatter (Bubble)":
                    proj_opts = ["Most Dominant", "Atomic", "Orbital", "Element-Orbital"]
                elif fb_style == "Lines (o_)":
                    proj_opts = ["Layer", "Normal", "Atomic", "Orbital", "Element-Orbital"]
                elif fb_style == "Heatmap (heat_)":
                    proj_opts = ["Total", "Atomic", "Orbital", "Element-Orbital"]

                fb_proj = st.selectbox("Projection Type", proj_opts, help="How to group projected orbitals", **settings_change)

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
                    st.caption(
                        "Assign atoms to the two layers. Their material formulas "
                        "are inferred automatically for the plot legend."
                    )
                    a_list, _, _, _ = get_available_channels(paths.get('fatband_dir'))

                    st.markdown("**(Optional) Auto-Assign from Structure**")
                    f_struc = st.file_uploader("Upload .in or .out file", key="u_struc_layer", help="Upload a QE structure file to automatically detect layers based on Z-coordinates", **settings_change)

                    auto_top, auto_bot = [], []
                    auto_bottom_name, auto_top_name = None, None
                    if f_struc:
                        try:
                            from qeplotter.analysis.bilayer import parse_qe_block
                            from qeplotter.core.utils import strip_number
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
                                labels = [
                                    f"{strip_number(sp)}{i+1}"
                                    for i, sp in enumerate(species)
                                ]
                                for i, sz in enumerate(shifted_z):
                                    if labels[i] in a_list:
                                        if sz > median_z:
                                            auto_top.append(labels[i])
                                        else:
                                            auto_bot.append(labels[i])
                                auto_map = {
                                    **{atom: 'top' for atom in auto_top},
                                    **{atom: 'bottom' for atom in auto_bot},
                                }
                                auto_bottom_name, auto_top_name = _layer_material_labels(
                                    auto_top + auto_bot, auto_map
                                )
                                st.success(
                                    f"Layers detected: {auto_bottom_name} → "
                                    f"{auto_top_name}"
                                )
                        except Exception as e:
                            st.warning(f"Could not auto-detect layers: {e}")

                    if not a_list:
                        st.warning("Please upload PDOS files in the Data tab to enable layer mapping.")
                    else:
                        top_title = "Upper layer atoms"
                        if auto_top_name:
                            top_title += f" · {auto_top_name}"
                        top_atoms = st.multiselect(top_title, a_list, default=auto_top, **settings_change)
                        bot_bot_options = [x for x in a_list if x not in top_atoms]
                        valid_auto_bot = [x for x in auto_bot if x in bot_bot_options]
                        bottom_title = "Lower layer atoms"
                        if auto_bottom_name:
                            bottom_title += f" · {auto_bottom_name}"
                        bot_atoms = st.multiselect(
                            bottom_title, bot_bot_options, default=valid_auto_bot,
                            **settings_change,
                        )

                        l_map = {}
                        for a in top_atoms:
                            l_map[a] = 'top'
                        for a in bot_atoms:
                            l_map[a] = 'bottom'
                        if l_map:
                            args['layer_assignment'] = l_map
                            if len(l_map) == len(a_list):
                                bottom_name, top_name = _layer_material_labels(
                                    a_list, l_map
                                )
                                st.info(
                                    f"Plot colour scale: **{bottom_name}** → "
                                    f"mixed → **{top_name}**"
                                )
                        args['data_note'] = st.text_input(
                            "Data source note (optional)",
                            help=(
                                "Printed below the figure. Identify synthetic, "
                                "demonstration, or unverified data explicitly."
                            ),
                            placeholder="e.g. QE calculation: PBE, 12×12×1 mesh",
                            **settings_change,
                        )

                # Dynamic Highlight Channels List
                atoms, elements, orbitals, exp_orbs = get_available_channels(paths.get('fatband_dir'))
                hl_options = elements  # default
                if base_m == 'orbital':
                    hl_options = orbitals
                elif base_m == 'element_orbital':
                    hl_options = exp_orbs
                elif base_m == 'atomic':
                    hl_options = elements

                # Inject generic if lists are empty (e.g. before upload)
                if not hl_options:
                    hl_options = ["Mo", "S", "d", "p", "Mo-d"]

                if fb_mode in ['o_orbital', 'o_atomic', 'o_element_orbital']:
                    args['dual'] = st.checkbox("Dual Channel Mode", help="Highlight two contrasting channels with a diverging colormap", **settings_change)

                if args.get('dual'):
                    c_h1, c_h2 = st.columns(2)
                    idx2 = 1 if len(hl_options) > 1 else 0
                    h1 = c_h1.selectbox("Channel 1", hl_options, index=0, **settings_change)
                    h2 = c_h2.selectbox("Channel 2", hl_options, index=idx2, **settings_change)
                    args['highlight_channel'] = (h1, h2)
                elif "heat" in fb_mode or fb_mode in ['normal', 'most', 'o_orbital', 'o_atomic', 'o_element_orbital']:
                    args['highlight_channel'] = st.selectbox("Highlight Channel", hl_options, index=0, help="Specific element/orbital to highlight", **settings_change)

                if "heat" in fb_mode:
                    args['overlay_bands_in_heat'] = st.checkbox("Overlay Lines", True, help="Add line bands on top of heatmap", **settings_change)

            if pt in ["band", "fatbands"]:
                args['plot_total_dos'] = st.checkbox("Plot Total DOS side-by-side", value=False, help="Requires DOS file uploaded", **settings_change)

            if pt in ["band", "fatbands"]:
                args['show_band_gap'] = st.checkbox("Show Band Gap Arrow", value=False, help="Detect and annotate the band gap (VBM → CBM) on the plot", **settings_change)

        # --- C. PLOT MODE & STYLE ---
        with tab_style:
            st.markdown("##### Plot Dimensions & Limits")
            col_w, col_h = st.columns(2)
            fig_width = col_w.number_input(
                "Width (in)", min_value=4.0, max_value=24.0,
                value=12.0, step=0.5,
                **style_change,
            )
            fig_height = col_h.number_input(
                "Height (in)", min_value=3.0, max_value=18.0,
                value=6.0, step=0.5,
                **style_change,
            )
            args['figsize'] = (float(fig_width), float(fig_height))
            args['dpi'] = st.number_input(
                "DPI", min_value=72, max_value=1200, value=200, step=25,
                **style_change,
            )

            c3, c4 = st.columns(2)
            if pt in ["pdos", "dos"]:
                if st.checkbox("Set Custom Y-Limits", value=False, **style_change):
                    args['y_range'] = (c3.number_input("Y-Min", value=0.0, min_value=-100.0, max_value=100.0, **style_change),
                                       c4.number_input("Y-Max", value=10.0, min_value=-100.0, max_value=100.0, **style_change))
                else:
                    args['y_range'] = None
            else:
                args['y_range'] = (c3.number_input("Y-Min", value=-3.0, min_value=-50.0, max_value=50.0, **style_change),
                                   c4.number_input("Y-Max", value=3.0, min_value=-50.0, max_value=50.0, **style_change))

            if pt == 'dos' or args.get('plot_total_dos', False):
                c5, c6 = st.columns(2)
                use_x = c5.checkbox("Set Custom DOS Limits (X-Axis)", value=False, **style_change)
                if use_x:
                    args['x_range'] = (0.0, c6.number_input("Max DOS Value", value=10.0, **style_change))
                elif pt == 'dos':
                    use_x2 = st.checkbox("Set Custom Energy Limits", value=False, **style_change)
                    if use_x2:
                        args['x_range'] = (c5.number_input("Energy-Min", value=-10.0, **style_change),
                                           c6.number_input("Energy-Max", value=10.0, **style_change))
                    else:
                        args['x_range'] = None
                else:
                    args['x_range'] = None

            st.markdown("##### Text & Layout")
            title_col, title_toggle_col = st.columns([3, 1])
            args['show_title'] = title_toggle_col.checkbox(
                "Show title", value=True,
                **style_change,
            )
            custom_title = title_col.text_input(
                "Plot title",
                placeholder="Leave blank to use the automatic title",
                disabled=not args['show_title'],
                **style_change,
            )
            args['plot_title'] = custom_title.strip() or None

            label_x_col, label_y_col = st.columns(2)
            custom_x_label = label_x_col.text_input(
                "X-axis label", placeholder="Automatic",
                **style_change,
            )
            custom_y_label = label_y_col.text_input(
                "Y-axis label", placeholder="Automatic",
                **style_change,
            )
            args['x_label'] = custom_x_label.strip() or None
            args['y_label'] = custom_y_label.strip() or None

            grid_col, legend_col = st.columns(2)
            args['show_grid'] = grid_col.checkbox("Show grid", value=True, **style_change)
            args['show_legend'] = legend_col.checkbox(
                "Show legend / colour scale",
                value=True,
                help=(
                    "Controls categorical legends and the continuous colour "
                    "scale used by line, layer, and heatmap fatbands."
                ),
                **style_change,
            )

            legend_location_col, legend_title_col = st.columns(2)
            fatband_mode = args.get('fatbands_mode', '')
            uses_colour_scale = pt == 'fatbands' and (
                fatband_mode == 'normal'
                or fatband_mode == 'layer'
                or fatband_mode.startswith('o_')
                or fatband_mode.startswith('heat_')
            )
            args['legend_location'] = legend_location_col.selectbox(
                "Legend position",
                [
                    "best", "upper right", "upper left", "lower left",
                    "lower right", "center right", "center left",
                    "lower center", "upper center", "center",
                ],
                disabled=not args['show_legend'] or uses_colour_scale,
                help=(
                    "Applies to categorical legends. Continuous colour "
                    "scales remain beside the plot."
                ),
                **style_change,
            )
            custom_legend_title = legend_title_col.text_input(
                "Legend / scale title",
                placeholder="Optional",
                disabled=not args['show_legend'],
                **style_change,
            )
            args['legend_title'] = custom_legend_title.strip() or None

            st.markdown("##### Colors & Visuals")
            uses_cmap = pt == 'fatbands' or (
                pt == 'band' and args.get('band_mode', 'normal') != 'normal'
            )
            if uses_cmap:
                if uses_colour_scale:
                    cmap_options = tuple(
                        name for name in _COLORMAP_OPTIONS
                        if name not in _CATEGORICAL_COLORMAPS
                    )
                    default_cmap = (
                        'coolwarm' if fatband_mode == 'layer' else 'viridis'
                    )
                else:
                    cmap_options = _COLORMAP_OPTIONS
                    default_cmap = 'tab10'

                args['cmap_name'] = st.selectbox(
                    "Colormap",
                    cmap_options,
                    index=cmap_options.index(default_cmap),
                    format_func=_colormap_option_label,
                    help=(
                        "The preview below shows the exact colour order used "
                        "by the generated plot."
                    ),
                    **style_change,
                )
                st.markdown(
                    _colormap_preview_html(args['cmap_name']),
                    unsafe_allow_html=True,
                )
            else:
                args['cmap_name'] = 'tab10'
                st.caption(
                    "This plot mode uses fixed line colours; switch to a "
                    "projected band or fatband mode to select a colormap."
                )

            if pt == "fatbands":
                c_adv1, c_adv2 = st.columns(2)
                args['s_min'] = c_adv1.number_input("Min Bubble Size", 1.0, 100.0, 10.0, **style_change)
                args['s_max'] = c_adv2.number_input("Max Bubble Size", 10.0, 500.0, 100.0, **style_change)
                args['weight_threshold'] = c_adv1.number_input("Weight Threshold", 0.0, 1.0, 0.01, **style_change)

                if "heat_" in args.get('fatbands_mode', ''):
                    args['heat_vmin'] = c_adv1.number_input("Heatmap Min Value", value=0.0, **style_change)
                    args['heat_vmax'] = c_adv2.number_input("Heatmap Max Value", value=0.0, **style_change)
                    if args['heat_vmax'] == 0.0:
                        args['heat_vmax'] = None

            if pt == "overlay_band":
                st.caption("Overlay Appearance")
                c_o1, c_o2 = st.columns(2)
                args['label1'] = c_o1.text_input("Label 1", "System A", **style_change)
                args['color1'] = c_o1.color_picker("Color 1", "#557A9E", **style_change)
                args['label2'] = c_o2.text_input("Label 2", "System B", **style_change)
                args['color2'] = c_o2.color_picker("Color 2", "#B06A63", **style_change)

    # ==========================================
    # RIGHT COLUMN: EXECUTION
    # ==========================================
    with col_preview:
        st.subheader("Generate")

        if st.button("Generate visualization", type="primary"):
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
            if pt == "dos" and not args.get('dos_file'):
                st.error("Total DOS file is required for DOS plotting. Upload it in the Data tab.")
                return
            if pt == "pdos" and not args.get('pdos_dir'):
                st.error("PDOS files are required for Projected DOS plotting. Upload them in the Data tab.")
                return
            if pt == "overlay_band" and (not args.get('band_file') or not args.get('kpath_file') or not args.get('band_file2') or not args.get('kpath_file2')):
                st.error("Both Band files and both K-Path files are required for overlay comparison.")
                return

            with st.spinner("Processing..."):
                try:
                    # Static Path (Standard)
                    log_io = StringIO()
                    plt.close('all')
                    with redirect_stdout(log_io):
                        plot_from_file(**args)

                    if plt.get_fignums():
                        fig = plt.gcf()

                        buf = BytesIO()
                        fig.savefig(buf, format="png", dpi=args.get('dpi', 200), bbox_inches='tight')
                        st.image(buf.getvalue(), caption="Plot Preview", use_container_width=True)

                        st.download_button("Download PNG", buf, "plot.png", "image/png")
                    else:
                        st.warning("No plot generated.")

                    # Show backend logs
                    log_text = log_io.getvalue()
                    if log_text.strip():
                        with st.expander("Backend logs", expanded=False):
                            st.code(log_text)

                except Exception as e:
                    st.error(f"Error: {e}")
                    import traceback
                    st.code(traceback.format_exc())
