"""
Fatband plotting functions.
Extracted verbatim from qep.py (plot_fatbands — all bubble, line/layer, heat modes).
"""
import os
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import numpy as np
from ase import Atoms

from qeplotter.core.io import read_band_xdistances, read_fatband_files
from qeplotter.core.utils import strip_number
from qeplotter.analysis.bandgap import _find_band_gap, _annotate_band_gap
from qeplotter.plotting.style import (
    apply_axis_style,
    apply_legend,
    display_text,
    figure_size,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _layer_material_labels(atom_names, layer_assignment):
    """Return readable lower/upper material names from assigned QE atom labels."""
    grouped = {"bottom": [], "top": []}
    for atom in atom_names:
        layer = layer_assignment.get(atom)
        if layer in grouped:
            grouped[layer].append(strip_number(atom))

    formulas = {}
    for layer, symbols in grouped.items():
        try:
            formula = Atoms(symbols=symbols).get_chemical_formula(
                mode="metal", empirical=True
            )
        except (KeyError, ValueError):
            formula = ""
        formulas[layer] = formula.translate(_SUBSCRIPT) if formula else layer.title()

    bottom, top = formulas["bottom"], formulas["top"]
    if bottom == top:
        return f"{bottom} · lower", f"{top} · upper"
    return bottom, top


def _layer_plot_title(bottom_material, top_material):
    """Build a concise, material-aware title for a bilayer band plot."""
    bottom_formula = bottom_material.split(" · ", 1)[0]
    top_formula = top_material.split(" · ", 1)[0]
    if bottom_formula == top_formula:
        material_name = f"{bottom_formula} bilayer"
    else:
        material_name = f"{bottom_formula} / {top_formula}"
    return f"{material_name} · layer-resolved fatbands"


def _create_plot_axes(plot_total_dos, dpi, figsize):
    """Create the shared fatband layout with an optional DOS side panel."""
    resolved_figsize = figure_size(
        figsize, (10, 6) if plot_total_dos else (8, 6)
    )
    if plot_total_dos:
        return plt.subplots(
            1,
            2,
            gridspec_kw={'width_ratios': [3, 1]},
            figsize=resolved_figsize,
            dpi=dpi,
            sharey=True,
        )
    fig, ax = plt.subplots(1, 1, figsize=resolved_figsize, dpi=dpi)
    return fig, (ax, None)

def plot_fatbands(
    fatband_dir,
    kpath_file,
    band_file,
    mode='most',
    highlight_channel=None,
    dual=False,
    fermi_level=None,
    shift_fermi=False,
    y_range=None,
    cmap_name='tab10',
    s_min=10,
    s_max=100,
    weight_threshold=0.01,
    plot_total_dos=False,
    dos_file=None,
    overlay_bands_in_heat=False,
    heat_vmin=None,
    heat_vmax=None,
    dpi=None,
    layer_assignment=None ,
        save_dir="saved",
        savefig=None,
        spin=False,
        sub_orb=False,
        x_range=None,
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
      Plot fatbands from Quantum ESPRESSO data.

      This function visualizes band structures with atomic/orbital-resolved "fatbands"
      from QE projwfc files, supporting several visualization modes.


    Parameters
    ----------
    fatband_dir : str
        Directory containing fatband (PDOS/projwfc) files produced by Quantum ESPRESSO.
    kpath_file : str
        Path to the QE K_POINTS file in "crystal_b" format (defines k-point path and high-symmetry labels).
    band_file : str
        Path to the QE band structure file (usually ends with .bands.dat.gnu).
    mode : str, optional
        Visualization mode for fatbands. Controls how channel-resolved information is displayed:
          - 'most', 'atomic', 'orbital', 'element_orbital': Bubble modes (dominant channel per (k,E) as colored bubbles)
          - 'normal', 'o_atomic', 'o_orbital', 'o_element_orbital': Line modes (band color encodes channel fraction)
          - 'heat_total', 'heat_atomic', 'heat_orbital', 'heat_element_orbital': Heatmap modes (background color = channel weight)
          - 'layer': Color bands by layer (requires `layer_assignment`).
    highlight_channel : str or list, optional
        Channel(s) to highlight in line/heatmap modes.
    dual : bool, optional
        If True (with a list of two highlight channels), uses a colorbar to interpolate between the two channels along each band.
    fermi_level : float, optional
        Value of the Fermi energy (in eV).
    shift_fermi : bool, optional
        If True, shifts all band energies and DOS so that the Fermi level appears at 0 eV on the plot.
    y_range : tuple or list, optional
        (ymin, ymax) values for the energy axis (eV).
    cmap_name : str, optional
        Name of the matplotlib colormap to use for coloring channels/bands (default: 'tab10').
    s_min : float, optional
        Minimum marker size for bubbles or heatmap points (default: 10).
    s_max : float, optional
        Maximum marker size for bubbles or heatmap points (default: 100).
    weight_threshold : float, optional
        Fraction (0–1) of the global max channel weight to be plotted.
    plot_total_dos : bool, optional
        If True, plots the total DOS (from `dos_file`) alongside the fatbands.
    dos_file : str, optional
        Path to the total DOS file.
    overlay_bands_in_heat : bool, optional
        If True, overlays plain band structure lines on top of heatmap plots.
    heat_vmin : float, optional
        Minimum value for heatmap color normalization.
    heat_vmax : float, optional
        Maximum value for heatmap color normalization.
    dpi : int, optional
        Output resolution (dots per inch).
    layer_assignment : dict, optional
        Only for 'layer' mode. Dictionary mapping atom names to 'top' or 'bottom' layer.
    data_note : str, optional
        Provenance or status note printed below the plot. Use this to identify
        demonstration, synthetic, or otherwise unverified data explicitly.
    save_dir : str, optional
        Directory to save the plot. Default "saved".
    savefig : str, optional
        Filename for saving the generated plot.

      Returns
      -------
      None. Plots and saves the figure.
      """


    labels, uniq_ik, E_grid, W_grids = read_fatband_files(fatband_dir,spin,sub_orb)
    N_k, N_e = E_grid.shape
    if shift_fermi and fermi_level is not None:
        E_grid = E_grid - fermi_level
    x_dist, band_energies, tick_positions, tick_labels, seg_ranges = read_band_xdistances(band_file, kpath_file)
    if len(x_dist) != N_k:
        print(f"Warning: fatband N_k={N_k} vs band file x_dist length={len(x_dist)}. They should match.")
    if plot_total_dos:
        if dos_file is None:
            raise ValueError("dos_file must be provided when plot_total_dos=True")
        dos_data = np.loadtxt(dos_file)
        if dos_data.ndim != 2 or dos_data.shape[1] < 2:
            raise ValueError(f"Unexpected DOS file format: {dos_file}")
        E_dos = dos_data[:,0]
        DOS   = dos_data[:,1]
        if shift_fermi and fermi_level is not None:
            E_dos = E_dos - fermi_level
    atom_labels = [a for (a, _) in labels]
    element_labels = [strip_number(a) for (a, _) in labels]
    orbs = [o for (_, o) in labels]
    element_orbital_labels = [
        f"{element}-{orb}"
        for element, orb in zip(element_labels, orbs)
    ]
    atom_orbital_labels = [
        f"{atom}-{orb}" for atom, orb in zip(atom_labels, orbs)
    ]
    bubble_modes = {'most','atomic','orbital','element_orbital'}
    line_modes = {'normal','o_atomic','o_orbital','o_element_orbital'}
    heat_modes = {'heat_total','heat_atomic','heat_orbital','heat_element_orbital'}


    if mode in bubble_modes:
        # Build grouped weights Wg
        if mode == 'atomic':
            unique_keys = sorted(set(atom_labels))
            Wg = np.zeros((len(unique_keys), N_k, N_e))
            for i, key in enumerate(unique_keys):
                for idx, a in enumerate(atom_labels):
                    if a == key:
                        Wg[i] += W_grids[idx]
        elif mode == 'orbital':
            unique_keys = sorted(set(orbs))
            Wg = np.zeros((len(unique_keys), N_k, N_e))
            for i, key in enumerate(unique_keys):
                for idx, o in enumerate(orbs):
                    if o == key:
                        Wg[i] += W_grids[idx]
        elif mode == 'element_orbital':
            unique_keys = sorted(set(element_orbital_labels))
            Wg = np.zeros((len(unique_keys), N_k, N_e))
            for i, key in enumerate(unique_keys):
                for idx, lab in enumerate(element_orbital_labels):
                    if lab == key:
                        Wg[i] += W_grids[idx]
        else:  # 'most'
            unique_keys = sorted(set(atom_orbital_labels))
            Wg = np.zeros((len(unique_keys), N_k, N_e))
            for i, key in enumerate(unique_keys):
                for idx, lab in enumerate(atom_orbital_labels):
                    if lab == key:
                        Wg[i] += W_grids[idx]
        # Determine dominant channel and weight
        idx_max = np.argmax(Wg, axis=0)  # shape (N_k, N_e)
        val_max = np.max(Wg, axis=0)
        # Flatten
        X_flat = np.repeat(x_dist, N_e)
        E_flat = E_grid.flatten()
        idx_flat = idx_max.flatten()
        val_flat = val_max.flatten()
        # Threshold
        global_max = np.nanmax(val_flat)
        thr = weight_threshold * global_max
        mask = val_flat > thr
        X_plot = X_flat[mask]
        E_plot = E_flat[mask]
        idx_plot = idx_flat[mask]
        val_plot = val_flat[mask]
        # Colors/sizes
        cmap = plt.get_cmap(cmap_name, len(unique_keys))
        colors = [cmap(i) for i in idx_plot]
        sizes = s_min + (s_max - s_min) * (val_plot / global_max if global_max>0 else 0)
        # Setup figure
        fig, (ax1, ax2) = _create_plot_axes(plot_total_dos, dpi, figsize)
        # Scatter
        ax1.scatter(X_plot, E_plot, s=sizes, c=colors, edgecolor='k', lw=0.3, alpha=0.8, zorder=1)
        # Overlay band lines (split segments)
        for band in band_energies:
            for (s,e) in seg_ranges:
                if e > s:
                    y = band[s:e+1]
                    x = x_dist[s:e+1]
                    if shift_fermi and fermi_level is not None:
                        y = y - fermi_level
                    ax1.plot(x, y, color='gray', lw=0.5, zorder=0)
                else:
                    x = x_dist[s:s+1]
                    y = band[s:s+1]
                    if shift_fermi and fermi_level is not None:
                        y = y - fermi_level
                    ax1.plot(x, y, 'o', color='gray', markersize=2, zorder=0)
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels(tick_labels)
        if y_range:
            ax1.set_ylim(y_range)
        title_mode = mode.capitalize() if mode!='most' else 'Most'
        ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'
        apply_axis_style(
            ax1,
            default_title=f'Fatbands ({title_mode})',
            default_xlabel='K-point Path',
            default_ylabel=ylabel,
            plot_title=plot_title,
            x_label=x_label,
            y_label=y_label,
            show_title=show_title,
            show_grid=show_grid,
        )
        # Legend
        for i, key in enumerate(unique_keys):
            ax1.scatter([], [], c=[cmap(i)], label=key, edgecolor='k', lw=0.3)
        apply_legend(
            ax1,
            show_legend=show_legend,
            location=legend_location,
            title=legend_title,
            fontsize='small',
            ncol=2,
        )
        # Total DOS panel
        if plot_total_dos:
            ax2.plot(DOS, E_dos, 'k-', lw=1)
            ax2.set_xlabel('DOS')
            if y_range:
                ax2.set_ylim(y_range)
            if x_range:
                ax2.set_xlim(x_range)
            ax2.axvline(0, color='gray', ls='--', lw=0.8)
            if show_grid:
                ax2.grid(True, ls='--', alpha=0.3)
            else:
                ax2.grid(False)

        # --- BAND GAP ANNOTATION ---
        if show_band_gap:
            _be = band_energies - fermi_level if (shift_fermi and fermi_level is not None) else band_energies
            gap_info = _find_band_gap(x_dist, _be, fermi_level, shift_fermi, scf_file=scf_file)
            _annotate_band_gap(ax1, gap_info)

        plt.tight_layout()
        if savefig:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            out = os.path.join(save_dir, os.path.basename(savefig))
            plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
            print(f"Saved figure to {out}")

        plt.show()
        return


    elif mode in line_modes or mode == 'layer':
        if mode == 'layer':
            # 1. Unique atom names from labels (now like 'W1', 'Mo2', 'S3' ...)
            atom_names = sorted(set(a for (a, _) in labels))

            # 2. Use provided assignment, raise if not provided or incomplete
            if layer_assignment is None:
                raise ValueError("You must provide layer_assignment as a dict when using mode='layer'.")

            # 3. Validate all atoms are in assignment
            for atom in atom_names:
                if atom not in layer_assignment:
                    raise ValueError(f"layer_assignment is missing entry for atom '{atom}'. "
                                     f"Please supply all atoms: {atom_names}")

            # 4. Sum weights by layer
            W_top = np.zeros((N_k, N_e))
            W_bottom = np.zeros((N_k, N_e))
            for i, (a, o) in enumerate(labels):
                if layer_assignment[a] == 'top':
                    W_top += W_grids[i]
                elif layer_assignment[a] == 'bottom':
                    W_bottom += W_grids[i]
                else:
                    raise ValueError(
                        f"layer_assignment for atom '{a}' must be 'top' or 'bottom', not '{layer_assignment[a]}'.")

            W_sum = W_top + W_bottom
            W_sum[W_sum <= 0] = np.nan
            frac = W_top / W_sum  # 1 = all top, 0 = all bottom
            bottom_material, top_material = _layer_material_labels(
                atom_names, layer_assignment
            )
            colorbar_label = 'Layer character'

        elif mode in line_modes:

            if dual:

                if isinstance(highlight_channel, str):

                    groups = [g.strip() for g in highlight_channel.split(',')]

                elif isinstance(highlight_channel, (list, tuple)):

                    groups = list(highlight_channel)

                else:

                    raise ValueError(

                        "For dual=True, highlight_channel must be 'g1,g2' or a two-element list/tuple"

                    )

                if len(groups) != 2:
                    raise ValueError(f"dual mode needs exactly two groups, got {groups!r}")

                key1, key2 = groups

                if mode == 'o_atomic':

                    valid = sorted(set(atom_labels))



                elif mode == 'o_orbital':

                    valid = sorted(set(orbs))

                else:

                    valid = sorted(set(element_orbital_labels))

                if key1 not in valid or key2 not in valid:
                    raise ValueError(f"dual keys {groups!r} must be among {valid}")

                W1 = np.zeros((N_k, N_e))

                W2 = np.zeros((N_k, N_e))

                if mode == 'o_atomic':
                    for i, a in enumerate(atom_labels):
                        if a == key1:
                            W1 += W_grids[i]
                        elif a == key2:
                            W2 += W_grids[i]

                elif mode == 'o_orbital':

                    for i, o in enumerate(orbs):

                        if o == key1:
                            W1 += W_grids[i]

                        elif o == key2:
                            W2 += W_grids[i]

                else:

                    for i, lab in enumerate(element_orbital_labels):

                        if lab == key1:
                            W1 += W_grids[i]

                        elif lab == key2:
                            W2 += W_grids[i]

                W12 = W1 + W2

                W12_safe = W12.copy()

                W12_safe[W12_safe <= 0] = np.nan

                frac = W2 / W12_safe

                colorbar_label = f'Fraction of {key2}   (0={key1}, 1={key2})'

            else:

                if mode == 'o_atomic':

                    if highlight_channel is None:
                        raise ValueError("highlight_channel must be provided for o_atomic mode")

                    unique_atoms = sorted(set(atom_labels))

                    if highlight_channel not in unique_atoms:
                        raise ValueError(f"highlight_channel '{highlight_channel}' not in atomic keys {unique_atoms}")

                    Wtot = np.zeros((N_k, N_e));
                    Whigh = np.zeros((N_k, N_e))

                    for idx, a in enumerate(atom_labels):

                        Wtot += W_grids[idx]

                        if a == highlight_channel:
                            Whigh += W_grids[idx]

                elif mode == 'o_orbital':

                    if highlight_channel is None:
                        raise ValueError("highlight_channel must be provided for o_orbital mode")

                    unique_orbs = sorted(set(orbs))

                    if highlight_channel not in unique_orbs:
                        raise ValueError(f"highlight_channel '{highlight_channel}' not in orbital keys {unique_orbs}")

                    Wtot = np.zeros((N_k, N_e));
                    Whigh = np.zeros((N_k, N_e))

                    for idx, o in enumerate(orbs):

                        Wtot += W_grids[idx]

                        if o == highlight_channel:
                            Whigh += W_grids[idx]

                elif mode == 'o_element_orbital':

                    if highlight_channel is None:
                        raise ValueError("highlight_channel must be provided for o_element_orbital mode")

                    unique_eo = sorted(set(element_orbital_labels))

                    if highlight_channel not in unique_eo:
                        raise ValueError(
                            f"highlight_channel '{highlight_channel}' not in element-orbital keys {unique_eo}")

                    Wtot = np.zeros((N_k, N_e));
                    Whigh = np.zeros((N_k, N_e))

                    for idx, lab in enumerate(element_orbital_labels):

                        Wtot += W_grids[idx]

                        if lab == highlight_channel:
                            Whigh += W_grids[idx]

                else:  # normal

                    Wtot = np.ones((N_k, N_e))

                    Whigh = np.zeros((N_k, N_e))

                Wtot_safe = Wtot.copy()

                Wtot_safe[Wtot_safe <= 0] = np.nan

                frac = Whigh / Wtot_safe

                colorbar_label = f'Fraction of {highlight_channel}'

        # ------ PLOTTING PART (shared for all line/layer modes) ------

        cmap = plt.get_cmap(cmap_name)

        norm = plt.Normalize(0.0, 1.0)

        nbands = band_energies.shape[0]

        fig, (ax1, ax2) = _create_plot_axes(plot_total_dos, dpi, figsize)

        for b in range(nbands):

            y_line = band_energies[b].copy()

            if shift_fermi and fermi_level is not None:
                y_line = y_line - fermi_level

            x_line = x_dist

            for (s, e) in seg_ranges:

                if e <= s:
                    continue

                xs = x_line[s:e + 1];
                ys = y_line[s:e + 1]

                points = np.array([xs, ys]).T.reshape(-1, 1, 2)

                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                frac_vals = []

                for i_k in range(s, e + 1):

                    Eb = band_energies[b, i_k]

                    Eb0 = Eb - fermi_level if (shift_fermi and fermi_level is not None) else Eb

                    row = E_grid[i_k, :]

                    j = np.argmin(np.abs(row - Eb0))

                    fv = frac[i_k, j]

                    if np.isnan(fv): fv = 0.0

                    frac_vals.append(fv)

                frac_seg = 0.5 * (np.array(frac_vals[:-1]) + np.array(frac_vals[1:]))

                lc = mcoll.LineCollection(
                    segments, array=frac_seg, cmap=cmap, norm=norm,
                    linewidth=2.2, capstyle='round', zorder=2,
                )

                ax1.add_collection(lc)

        ax1.set_xticks(tick_positions)

        display_tick_labels = [
            'Γ' if str(label).strip().upper() in {'G', 'GAMMA'} else label
            for label in tick_labels
        ]
        ax1.set_xticklabels(display_tick_labels)

        ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'

        if y_range:
            ax1.set_ylim(y_range)

        if mode == 'layer':
            default_title = _layer_plot_title(bottom_material, top_material)

        else:

            title_mode = mode if mode != 'normal' else f"Highlight {highlight_channel}"

            default_title = f'Fatbands ({title_mode})'

        apply_axis_style(
            ax1,
            default_title=default_title,
            default_xlabel='K-point Path',
            default_ylabel=ylabel,
            plot_title=plot_title,
            x_label=x_label,
            y_label=y_label,
            show_title=show_title,
            show_grid=show_grid,
            grid_alpha=0.22,
        )

        if show_legend:
            scale_title = display_text(legend_title, colorbar_label)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

            sm.set_array([])

            cbar = plt.colorbar(sm, ax=ax1, pad=0.02)

            cbar.set_label(scale_title)

            if mode == 'layer':
                cbar.set_ticks([0.0, 0.5, 1.0])
                cbar.set_ticklabels([bottom_material, 'Mixed', top_material])
                cbar.set_label('')
                cbar.ax.set_title(scale_title, fontsize=10, pad=8)
                cbar.ax.tick_params(length=0, pad=6, labelsize=10)
                cbar.outline.set_linewidth(0.7)

        if overlay_bands_in_heat:

            for band in band_energies:

                for (s, e) in seg_ranges:

                    if e > s:

                        y = band[s:e + 1];
                        x = x_dist[s:e + 1]

                        if shift_fermi and fermi_level is not None:
                            y = y - fermi_level

                        ax1.plot(x, y, color='lightgray', lw=0.5, zorder=0)

                    else:

                        x = x_dist[s];
                        y = band[s]

                        if shift_fermi and fermi_level is not None:
                            y = y - fermi_level

                        ax1.plot(x, y, color='lightgray', lw=0.5, zorder=0)

        ax1.set_axisbelow(True)
        ax1.margins(x=0)
        if shift_fermi and fermi_level is not None:
            ax1.axhline(0.0, color='#555555', ls='--', lw=0.9, alpha=0.8, zorder=0)



        # Total DOS panel for line/layer modes
        if plot_total_dos:
            if shift_fermi and fermi_level is not None:
                E_dos_plot = E_dos
            else:
                E_dos_plot = E_dos
            ax2.plot(DOS, E_dos_plot, 'k-', lw=1)
            ax2.set_xlabel('DOS')
            if y_range:
                ax2.set_ylim(y_range)
            if x_range:
                ax2.set_xlim(x_range)
            ax2.axvline(0, color='gray', ls='--', lw=0.8)
            if show_grid:
                ax2.grid(True, ls='--', alpha=0.3)
            else:
                ax2.grid(False)

        # --- BAND GAP ANNOTATION ---
        if show_band_gap:
            _be = band_energies - fermi_level if (shift_fermi and fermi_level is not None) else band_energies
            gap_info = _find_band_gap(x_dist, _be, fermi_level, shift_fermi, scf_file=scf_file)
            _annotate_band_gap(ax1, gap_info)

        if data_note:
            fig.text(
                0.01, 0.01, str(data_note), ha='left', va='bottom',
                fontsize=7.5, color='#666666',
            )
            fig.tight_layout(rect=(0, 0.035, 1, 1))
        else:
            fig.tight_layout()
        if savefig:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            out = os.path.join(save_dir, os.path.basename(savefig))
            plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
            print(f"Saved figure to {out}")

        plt.show()

    elif mode in heat_modes:
        # Heatmap modes: show intensity (weight) as colored background along bands
        if mode == 'heat_total':
            # sum all channels
            Wtot = np.zeros((N_k, N_e))
            for arr in W_grids:
                Wtot += arr
            heat_grid = Wtot
            heat_label = 'Total weight'
        else:
            if highlight_channel is None:
                raise ValueError(f"highlight_channel must be provided for {mode}")
            if mode == 'heat_atomic':
                unique_atoms = sorted(set(atom_labels))
                if highlight_channel not in unique_atoms:
                    raise ValueError(f"highlight_channel '{highlight_channel}' not in atomic keys {unique_atoms}")
                Whigh = np.zeros((N_k, N_e))
                for idx, a in enumerate(atom_labels):
                    if a == highlight_channel:
                        Whigh += W_grids[idx]
                heat_grid = Whigh
                heat_label = f"Weight of atom {highlight_channel}"
            elif mode == 'heat_orbital':
                unique_orbs = sorted(set(orbs))
                if highlight_channel not in unique_orbs:
                    raise ValueError(f"highlight_channel '{highlight_channel}' not in orbital keys {unique_orbs}")
                Whigh = np.zeros((N_k, N_e))
                for idx, o in enumerate(orbs):
                    if o == highlight_channel:
                        Whigh += W_grids[idx]
                heat_grid = Whigh
                heat_label = f"Weight of orbital {highlight_channel}"
            elif mode == 'heat_element_orbital':
                unique_eo = sorted(set(element_orbital_labels))
                if highlight_channel not in unique_eo:
                    raise ValueError(f"highlight_channel '{highlight_channel}' not in element-orbital keys {unique_eo}")
                Whigh = np.zeros((N_k, N_e))
                for idx, lab in enumerate(element_orbital_labels):
                    if lab == highlight_channel:
                        Whigh += W_grids[idx]
                heat_grid = Whigh
                heat_label = f"Weight of {highlight_channel}"
            else:
                raise ValueError(f"Unknown heat mode: {mode}")
        # Setup figure
        fig, (ax1, ax2) = _create_plot_axes(plot_total_dos, dpi, figsize)
        # Flatten grid
        X_flat = np.repeat(x_dist, N_e)
        E_flat = E_grid.flatten()
        W_flat = heat_grid.flatten()
        # Color normalization
        vmin = heat_vmin if heat_vmin is not None else np.nanmin(W_flat)
        vmax = heat_vmax if heat_vmax is not None else np.nanmax(W_flat)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(cmap_name)
        # Plot each point as small square marker
        mask = ~np.isnan(W_flat)
        Xp = X_flat[mask]
        Ep = E_flat[mask]
        Wp = W_flat[mask]
        colors = cmap(norm(Wp))
        s_heat = (s_min + s_max) / 2.0
        ax1.scatter(Xp, Ep, s=s_heat, c=colors, marker='s', edgecolor='none', alpha=1.0, zorder=0)
        # Optionally overlay band lines
        if overlay_bands_in_heat:
            for band in band_energies:
                for (s,e) in seg_ranges:
                    if e > s:
                        y = band[s:e+1]
                        x = x_dist[s:e+1]
                        if shift_fermi and fermi_level is not None:
                            y = y - fermi_level
                        ax1.plot(x, y, color='lightgray', lw=0.5, zorder=1)
                    else:
                        x = x_dist[s]
                        y = band[s]
                        if shift_fermi and fermi_level is not None:
                            y = y - fermi_level
                        ax1.plot(x, y, 'o', color='lightgray', markersize=2, zorder=1)
        # X-ticks
        ax1.set_xticks(tick_positions)
        ax1.set_xticklabels(tick_labels)
        ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'
        if y_range:
            ax1.set_ylim(y_range)
        # Title and colorbar
        mode_title = mode.replace('_',' ').title()
        apply_axis_style(
            ax1,
            default_title=f'Fatbands ({mode_title})',
            default_xlabel='K-point Path',
            default_ylabel=ylabel,
            plot_title=plot_title,
            x_label=x_label,
            y_label=y_label,
            show_title=show_title,
            show_grid=show_grid,
        )
        if show_legend:
            scale_title = display_text(legend_title, heat_label)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax1, pad=0.02)
            cbar.set_label(scale_title)
        # Total DOS panel
        if plot_total_dos:
            if shift_fermi and fermi_level is not None:
                E_dos_plot = E_dos
            else:
                E_dos_plot = E_dos
            ax2.plot(DOS, E_dos_plot, 'k-', lw=1)
            ax2.set_xlabel('DOS')
            if y_range:
                ax2.set_ylim(y_range)
            if x_range:
                ax2.set_xlim(x_range)
            ax2.axvline(0, color='gray', ls='--', lw=0.8)
            if show_grid:
                ax2.grid(True, ls='--', alpha=0.3)
            else:
                ax2.grid(False)

        # --- BAND GAP ANNOTATION ---
        if show_band_gap:
            _be = band_energies - fermi_level if (shift_fermi and fermi_level is not None) else band_energies
            gap_info = _find_band_gap(x_dist, _be, fermi_level, shift_fermi, scf_file=scf_file)
            _annotate_band_gap(ax1, gap_info)

        plt.tight_layout()
        if savefig:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            out = os.path.join(save_dir, os.path.basename(savefig))
            plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
            print(f"Saved figure to {out}")

        plt.show()

    else:
        raise ValueError(f"Unknown fatbands mode: {mode}")
