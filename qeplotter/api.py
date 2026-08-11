"""
High-level API: plot_from_file dispatcher.
Extracted verbatim from qep.py.
"""
from qeplotter.plotting.bands import plot_band, overlay_band_plot
from qeplotter.plotting.dos import plot_dos, plot_pdos_dir
from qeplotter.plotting.fatbands import plot_fatbands

def plot_from_file(
    plot_type='band',
    pdos_dir=None,
    fatband_dir=None,
    kpath_file=None,
    band_file=None,
    dos_file=None,
    pdos_mode='atomic',
    fatbands_mode='most',
    highlight_channel=None,
    dual=False,
    band_mode='normal',
    cmap_name='tab10',
    band_file2=None,
    kpath_file2=None,
    color1='red',
    color2='blue',
    label1='Band 1',
    label2='Band 2',
    s_min=10,
    s_max=100,
    weight_threshold=0.01,
    y_range=None,
    x_range=None,
    fermi_level=None,
    shift_fermi=False,
    plot_total_dos=False,
    overlay_bands_in_heat=False,
    heat_vmin=None,
    heat_vmax=None,
    dpi=None,
    layer_assignment=None,
    save_dir="saved",
    savefig=None,
    spin=False,
    sub_orb=False,
    vertical=False,
    show_band_gap=False,
    scf_file=None,
    data_note=None,
    figsize=None,
    plot_title=None,
    x_label=None,
    y_label=None,
    show_title=True,
    show_grid=True,
    show_legend=True,
    legend_location='best',
    legend_title=None,
):
    """
    High-level wrapper for plotting Quantum ESPRESSO band, DOS, PDOS, and fatbands in one function.

    Parameters
    ----------
    plot_type : str
        'band', 'dos', 'pdos', 'overlay_band', or 'fatbands'.
    (See individual plot functions for full parameter documentation.)
    """
    pt = plot_type.lower()
    common_style = {
        'figsize': figsize,
        'plot_title': plot_title,
        'x_label': x_label,
        'y_label': y_label,
        'show_title': show_title,
        'show_grid': show_grid,
        'show_legend': show_legend,
        'legend_location': legend_location,
        'legend_title': legend_title,
    }
    if pt == 'band':
        plot_band(
            band_file=band_file,
            kpath_file=kpath_file,
            fermi_level=fermi_level,
            shift_fermi=shift_fermi,
            y_range=y_range,
            dpi=dpi,
            band_mode=band_mode,
            fatband_dir=fatband_dir,
            cmap_name=cmap_name,
            save_dir=save_dir,
            savefig=savefig,
            spin=spin, sub_orb=sub_orb,
            plot_total_dos=plot_total_dos,
            dos_file=dos_file,
            x_range=x_range,
            show_band_gap=show_band_gap,
            scf_file=scf_file,
            **common_style,
        )
    elif pt == 'dos':
        plot_dos(
            dos_file, fermi_level, shift_fermi, y_range, x_range=x_range, dpi=dpi,
            save_dir=save_dir, savefig=savefig, vertical=vertical,
            **common_style,
        )
    elif pt == 'overlay_band':
        overlay_band_plot(
            band_file, kpath_file,
            band_file2, kpath_file2,
            fermi_level=fermi_level,
            shift_fermi=shift_fermi,
            y_range=y_range,
            dpi=dpi,
            color1=color1,
            color2=color2,
            label1=label1,
            label2=label2,
            save_dir=save_dir,
            savefig=savefig,
            **common_style,
        )
    elif pt == 'pdos':
        plot_pdos_dir(
            pdos_dir, fermi_level, shift_fermi, y_range, dpi=dpi, pdos_mode=pdos_mode,
            save_dir=save_dir, savefig=savefig,
            **common_style,
        )
    elif pt == 'fatbands':
        fb_dir = fatband_dir if fatband_dir is not None else pdos_dir
        if fb_dir is None or band_file is None or kpath_file is None:
            raise ValueError("fatband_dir (or file_path), band_file, and kpath_file are all required for 'fatbands'")
        plot_fatbands(
            fatband_dir=fb_dir,
            kpath_file=kpath_file,
            band_file=band_file,
            mode=fatbands_mode,
            highlight_channel=highlight_channel,
            dual=dual,
            fermi_level=fermi_level,
            shift_fermi=shift_fermi,
            y_range=y_range,
            cmap_name=cmap_name,
            s_min=s_min,
            s_max=s_max,
            weight_threshold=weight_threshold,
            plot_total_dos=plot_total_dos,
            dos_file=dos_file,
            overlay_bands_in_heat=overlay_bands_in_heat,
            heat_vmin=heat_vmin,
            heat_vmax=heat_vmax,
            dpi=dpi,
            layer_assignment=layer_assignment,
            save_dir=save_dir,
            savefig=savefig,
            spin=spin,
            sub_orb=sub_orb,
            x_range=x_range,
            show_band_gap=show_band_gap,
            scf_file=scf_file,
            data_note=data_note,
            **common_style,
        )
    else:
        raise ValueError("Use 'band','dos','pdos', 'overlay_band', or 'fatbands' for plot_type")

def launch_gui():
    """Launch the recommended modular Streamlit application."""
    import importlib.util
    import subprocess
    import sys
    from pathlib import Path

    app_spec = importlib.util.find_spec("gui_mod")
    if app_spec is None or app_spec.origin is None:
        raise FileNotFoundError(
            "Could not find the installed QEPlotter GUI module `gui_mod`."
        )
    app_path = Path(app_spec.origin)

    print(f"Launching QEPlotter from {app_path}...")
    try:
        subprocess.run(
            [sys.executable, "-m", "streamlit", "run", str(app_path)],
            check=True,
        )
    except Exception as e:
        print(f"Failed to launch GUI: {e}")
        print(f"Try running manually: streamlit run {app_path}")
