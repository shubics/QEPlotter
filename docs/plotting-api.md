# Plotting API reference

`qeplotter.plot_from_file()` is the maintained high-level entry point for band,
DOS/PDOS, overlay, and fatband plots.

```python
from qeplotter import plot_from_file
```

## Core parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `plot_type` | `str` | `'band'` | Plot family: `'band'`, `'dos'`, `'pdos'`, `'fatbands'`, or `'overlay_band'`. |
| `band_file` | `str` | `None` | QE band file (`.gnu`); required for band and fatband plots. |
| `kpath_file` | `str` | `None` | QE `K_POINTS crystal_b` path; required for band and fatband plots. |
| `fermi_level` | `float` | `None` | Fermi energy in eV. |
| `shift_fermi` | `bool` | `False` | Shift the supplied Fermi level to 0 eV. |
| `y_range` | `tuple` | `None` | Energy-axis limits, for example `(-5, 5)`. |
| `dpi` | `int` | `None` | Resolution for saved figures. |
| `save_dir` | `str` | `'saved'` | Output directory. |
| `savefig` | `str` | `None` | Output filename; if omitted, the plot is displayed only. |

## Figure styling

The same styling options are accepted by every plot family and are also
available in the GUI's **Plot Styling** tab. Empty custom text values preserve
QEPlotter's automatic labels.

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `figsize` | `tuple` | `None` | Figure width and height in inches; plot-specific defaults are used when omitted. |
| `plot_title` | `str` | `None` | Custom plot title. |
| `x_label` | `str` | `None` | Custom label for the main X axis. |
| `y_label` | `str` | `None` | Custom label for the main Y axis. |
| `show_title` | `bool` | `True` | Show or hide the plot title. |
| `show_grid` | `bool` | `True` | Show or hide grid lines, including the optional DOS panel. |
| `show_legend` | `bool` | `True` | Show categorical legends and continuous fatband colour scales. |
| `legend_location` | `str` | `'best'` | Matplotlib placement for categorical legends. |
| `legend_title` | `str` | `None` | Optional title for categorical legends or continuous fatband colour scales. |

## Band and fatband parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `band_mode` | `str` | `'normal'` | Band colouring: `'normal'`, `'atomic'`, `'orbital'`, `'element_orbital'`, or `'most'`. |
| `fatband_dir` | `str` | `None` | Directory containing `projwfc.x` projection files. |
| `fatbands_mode` | `str` | `'most'` | Fatband mode: bubble modes, `heat_*`, `o_*`, or `'layer'`. |
| `spin` | `bool` | `False` | Enable spin-polarized projection parsing. |
| `sub_orb` | `bool` | `False` | Split orbitals into available sub-components. |
| `highlight_channel` | `str/tuple` | `None` | Atom or orbital channel to highlight. |
| `dual` | `bool` | `False` | Compare two channels in a line fatband plot. |
| `layer_assignment` | `dict` | `None` | Map every projected atom to `'top'` or `'bottom'`. |
| `data_note` | `str` | `None` | Fatband provenance/status footer; use it to identify demonstration or unverified data. |
| `s_min` | `float` | `10` | Minimum bubble marker size. |
| `s_max` | `float` | `100` | Maximum bubble marker size. |
| `weight_threshold` | `float` | `0.01` | Ignore projection weights below this threshold. |
| `cmap_name` | `str` | `'tab10'` | Matplotlib colormap used by projected Band, DOS, PDOS, Overlay, and Fatbands plots. Normal bands retain solid black lines. |
| `show_band_gap` | `bool` | `False` | Detect and annotate VBM/CBM. |
| `scf_file` | `str` | `None` | Optional `scf.out` reference for gap detection. |

## DOS parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `dos_file` | `str` | `None` | `dos.x` output; required for total DOS. |
| `pdos_dir` | `str` | `None` | Directory containing PDOS files. |
| `pdos_mode` | `str` | `'atomic'` | Group PDOS by `'atomic'`, `'orbital'`, or `'element_orbital'`. |
| `plot_total_dos` | `bool` | `False` | Add a total-DOS side panel to band/fatband plots. |
| `x_range` | `tuple` | `None` | DOS-axis range. |
| `vertical` | `bool` | `False` | Put energy on the horizontal axis for DOS-only plots. |

## Heatmap parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `overlay_bands_in_heat` | `bool` | `False` | Draw band lines over a fatband heatmap. |
| `heat_vmin` | `float` | `None` | Lower heatmap normalization bound. |
| `heat_vmax` | `float` | `None` | Upper heatmap normalization bound. |

## Overlay parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| `band_file2` | `str` | `None` | Second band file. |
| `kpath_file2` | `str` | `None` | Second K-path file. |
| `label1` | `str` | `'Band 1'` | Legend label for the first band set. |
| `label2` | `str` | `'Band 2'` | Legend label for the second band set. |
| `color1` | `str` | `'red'` | First band colour. |
| `color2` | `str` | `'blue'` | Second band colour. |

The examples in the main README show typical calls. For scientific figures,
keep the corresponding QE inputs and record calculation provenance alongside
the generated output.
