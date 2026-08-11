"""
Band structure plotting functions.
Extracted verbatim from qep.py (plot_band, overlay_band_plot).
"""
import os
import re
import matplotlib.pyplot as plt
import numpy as np

from qeplotter.core.io import read_band_xdistances, read_fatband_files
from qeplotter.core.utils import strip_number
from qeplotter.analysis.bandgap import _find_band_gap, _annotate_band_gap
from qeplotter.plotting.style import (
    apply_axis_style,
    apply_legend,
    colormap_colors,
    figure_size,
)


def overlay_band_plot(
    band_file1, kpath_file1, band_file2, kpath_file2,
    fermi_level=None, shift_fermi=False,
    y_range=None, dpi=None, color1='red', color2='blue',
    label1='Band1', label2='Band2',
    save_dir="saved", savefig=None,
    figsize=None, plot_title=None, x_label=None, y_label=None,
    show_title=True, show_grid=True, show_legend=True,
    legend_location='best', legend_title=None,
    cmap_name='tab10',
):
    """
    Overlay two band structures on the same plot.
    """
    x1, bands1, ticks1, labels1, segs1 = read_band_xdistances(band_file1, kpath_file1)
    x2, bands2, ticks2, labels2, segs2 = read_band_xdistances(band_file2, kpath_file2)

    if shift_fermi and fermi_level is not None:
        bands1 = bands1 - fermi_level
        bands2 = bands2 - fermi_level

    palette = colormap_colors(cmap_name, 2)
    color1 = palette[0] if color1 is None else color1
    color2 = palette[1] if color2 is None else color2

    plt.figure(figsize=figure_size(figsize, (8, 6)), dpi=dpi)

    for i, band in enumerate(bands1):
        for (s, e) in segs1:
            plt.plot(x1[s:e+1], band[s:e+1], color=color1, lw=1, alpha=0.85, label=label1 if i == 0 else None)
    for i, band in enumerate(bands2):
        for (s, e) in segs2:
            plt.plot(x2[s:e+1], band[s:e+1], color=color2, lw=1, alpha=0.85, label=label2 if i == 0 else None)

    for tx in ticks1:
        plt.axvline(tx, color='gray', ls='--', alpha=0.5)
    plt.xticks(ticks1, labels1)
    ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'

    if y_range:
        plt.ylim(y_range)
    if fermi_level is not None:
        y0 = 0.0 if shift_fermi else fermi_level
        plt.axhline(y0, color='r', ls='--', lw=1.2, label=f'Fermi = {fermi_level:.2f} eV')
    ax = plt.gca()
    apply_axis_style(
        ax,
        default_title="Overlaid Band Structures",
        default_xlabel="K-point Path",
        default_ylabel=ylabel,
        plot_title=plot_title,
        x_label=x_label,
        y_label=y_label,
        show_title=show_title,
        show_grid=show_grid,
        grid_alpha=0.4,
    )
    apply_legend(
        ax,
        show_legend=show_legend,
        location=legend_location,
        title=legend_title,
    )
    plt.tight_layout()

    # --- Save ---
    def sanitize(s):
        return re.sub(r'\W+', '_', s).strip('_')
    
    if savefig:
        filename = savefig
    else:
        filename = f"BandStructure_{sanitize(label1)}_vs_{sanitize(label2)}.png"
        
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    full_path = os.path.join(save_dir, filename)
    plt.savefig(full_path, dpi=dpi, bbox_inches='tight')
    print(f"Plot saved as {full_path}")

    plt.show()


def plot_band(
    band_file,
    kpath_file,
    fermi_level=None,
    shift_fermi=False,
    y_range=None,
    dpi=None,
    band_mode='normal',
    fatband_dir=None,
    cmap_name='tab10',
    save_dir="saved",
    savefig=None, 
    spin=False,
    sub_orb=False,
    plot_total_dos=False,
    dos_file=None,
    x_range=None,
    show_band_gap=False,
    scf_file=None,
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
    Plot the electronic band structure from Quantum ESPRESSO.
    """
    import matplotlib.cm as cm
    # 1) Read band data
    x_dist, band_energies, tick_positions, tick_labels, seg_ranges = read_band_xdistances(band_file, kpath_file)
    N_k = x_dist.shape[0]

    if shift_fermi and fermi_level is not None:
        band_energies = band_energies - fermi_level

    # --- PREPARE DOS DATA IF REQUESTED ---
    if plot_total_dos:
        if dos_file is None:
            raise ValueError("dos_file must be provided when plot_total_dos=True")
        dos_data = np.loadtxt(dos_file)
        if dos_data.ndim != 2 or dos_data.shape[1] < 2:
            raise ValueError(f"Unexpected DOS file format: {dos_file}")
        E_dos = dos_data[:, 0]
        DOS = dos_data[:, 1]
        
        if shift_fermi and fermi_level is not None:
            E_dos = E_dos - fermi_level

    # --- SETUP FIGURE ---
    resolved_figsize = figure_size(
        figsize, (10, 6) if plot_total_dos else (8, 6)
    )
    if plot_total_dos:
        if dpi is not None:
            fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=resolved_figsize, dpi=dpi, sharey=True)
        else:
            fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=resolved_figsize, sharey=True)
    else:
        if dpi is not None:
            fig, ax1 = plt.subplots(1, 1, figsize=resolved_figsize, dpi=dpi)
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=resolved_figsize)
        ax2 = None

    if band_mode == 'normal' or band_mode is None:
        line_color = colormap_colors(cmap_name, 1)[0]
        for band in band_energies:
            for (s,e) in seg_ranges:
                if e > s:
                    ax1.plot(x_dist[s:e+1], band[s:e+1], color=line_color, lw=1)
                else:
                    ax1.plot(x_dist[s], band[s], '.', color=line_color, markersize=2)
        title = 'Band Structure'
    else:

        if fatband_dir is None:
            raise ValueError(f"band_mode='{band_mode}' requires fatband_dir with projection files")
        # 2) Read fatband projection grids
        labels, uniq_ik, E_grid, W_grids = read_fatband_files(fatband_dir,spin,sub_orb)

        if len(uniq_ik) != N_k:
            print(f"Warning: fatband N_k={len(uniq_ik)} vs band file N_k={N_k}. They should match for correct coloring.")

          # strip off the atom numbers so all 'Bi1','Bi2',... become just 'Bi'
        ch_labels = [f"{a}-{o}" for (a, o) in labels]
        elems = [a for (a, _) in labels]
        orbs  = [o for (_,o) in labels]

        if band_mode == 'atomic':
            unique_keys = sorted(set(elems))

            group_indices = {key: [i for i,a in enumerate(elems) if a==key] for key in unique_keys}
        elif band_mode == 'orbital':
            unique_keys = sorted(set(orbs))
            group_indices = {key: [i for i,o in enumerate(orbs) if o==key] for key in unique_keys}
        elif band_mode in ('element_orbital', 'most'):
            unique_keys = sorted(set(ch_labels))
            group_indices = {key: [i for i,lab in enumerate(ch_labels) if lab==key] for key in unique_keys}
        else:
            raise ValueError(f"Unknown band_mode: {band_mode}")

        # 3) For each band, determine dominant group:
        nbands = band_energies.shape[0]
        band_colors = [None]*nbands

        cmap = cm.get_cmap(cmap_name, len(unique_keys))

        for b in range(nbands):

            group_sums = {key: 0.0 for key in unique_keys}

            k_iter = min(N_k, E_grid.shape[0])
            for i in range(k_iter):
                Eb = band_energies[b,i]

                if shift_fermi and fermi_level is not None:
                    row = E_grid[i,:] - fermi_level
                else:
                    row = E_grid[i,:]

                if row.size == 0:
                    continue
                j = np.argmin(np.abs(row - Eb))

                for key, idx_list in group_indices.items():

                    try:
                        wvals = [W_grids[idx][i,j] for idx in idx_list]
                    except Exception:

                        continue
                    group_sums[key] += np.sum(wvals)

            sums = np.array([group_sums[k] for k in unique_keys])
            if sums.sum() <= 0:
                color = 'k'
            else:
                imax = np.argmax(sums)
                color = cmap(imax)
            
            # Plot colored band
            for (s,e) in seg_ranges:
                if e > s:
                    ax1.plot(x_dist[s:e+1], band_energies[b, s:e+1], color=color, lw=1.5)
                else:
                    ax1.plot(x_dist[s], band_energies[b, s], '.', color=color, markersize=3)
        
        # Legend for colored bands
        for i, key in enumerate(unique_keys):
            ax1.plot([], [], c=cmap(i), label=key, lw=2)
        title = f'Band Structure ({band_mode})'

    # --- COMMON PLOTTING ---
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels)
    # Dynamic Y-label based on shift
    ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'

    if y_range:
        ax1.set_ylim(y_range)
    apply_axis_style(
        ax1,
        default_title=title,
        default_xlabel='K-point Path',
        default_ylabel=ylabel,
        plot_title=plot_title,
        x_label=x_label,
        y_label=y_label,
        show_title=show_title,
        show_grid=show_grid,
    )
    ax1.autoscale(enable=True, axis='x', tight=True)
    
    # Add Fermi line to Band Plot
    if fermi_level is not None:
        y0 = 0.0 if shift_fermi else fermi_level
        ax1.axhline(y0, color='r', ls='--', lw=1.0)

    # --- TOTAL DOS PLOTTING ---
    if plot_total_dos:
        # Side-by-side means Energy on Y, DOS on X (Standard Vertical Layout)
        dos_color = colormap_colors(cmap_name, 1)[0]
        ax2.plot(DOS, E_dos, color=dos_color, lw=1, label='Total DOS')
        ax2.set_xlabel('DOS')
        ax2.set_title("Total DOS")
        
        # Share Y axis with band plot
        if y_range:
            ax2.set_ylim(y_range)
        if x_range:
            ax2.set_xlim(x_range)
            
        # Draw Fermi Line
        if fermi_level is not None:
             y0 = 0.0 if shift_fermi else fermi_level
             # Show original Fermi value in legend, consistent with plot_dos
             label_f = f'Fermi = {fermi_level:.2f} eV' 
             ax2.axhline(y0, color='r', ls='--', lw=1.2, label=label_f)
        else:
             # Just a zero line if no fermi info
             ax2.axhline(0, color='gray', ls='--', lw=0.8)

        if show_grid:
            ax2.grid(True, ls='--', alpha=0.4)
        else:
            ax2.grid(False)
        
        # Hide Y-axis labels on DOS plot since they are shared
        plt.setp(ax2.get_yticklabels(), visible=False)
        
        if fermi_level is not None:
             apply_legend(
                 ax2,
                 show_legend=show_legend,
                 location=legend_location,
                 title=legend_title,
                 fontsize='small',
             )

    apply_legend(
        ax1,
        show_legend=show_legend,
        location=legend_location,
        title=legend_title,
        fontsize='small',
        ncol=2,
    )

    # --- BAND GAP ANNOTATION ---
    if show_band_gap:
        gap_info = _find_band_gap(x_dist, band_energies, fermi_level, shift_fermi, scf_file=scf_file)
        _annotate_band_gap(ax1, gap_info)

    plt.tight_layout()
    if savefig:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        out = os.path.join(save_dir, os.path.basename(savefig))
        plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
        print(f"Saved figure to {out}")

    plt.show()
