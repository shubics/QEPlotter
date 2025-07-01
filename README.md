# QEPlotter

**Quantum ESPRESSO Band & Fatband Visualization Toolkit**

QEPlotter is a Python library designed to streamline the post-processing and visualization of Quantum ESPRESSO electronic structure outputs. With a unified `plot_from_file` API, you can generate publication‑quality plots for:

* **Band structures** (plain or projection‑colored)
* **Density of states (DOS/PDOS)**
* **Fatband projections** (bubble, line, heatmap, layer modes)

---

## 🚀 Installation

```bash
pip install qeplotter
```

Or to work from the source:

```bash
git clone https://github.com/shubics/QEPlotter.git
cd QEPlotter
pip install -e .
```

Requirements:

* Python 3.8 or newer
* `numpy`, `matplotlib`
  - (Optional) TeX Live for rendering LaTeX in axis labels
  - (Optional) Mathematica/Matlab/PDF tools for export

---

## 🎉 Quickstart

All plotting is driven by a single function: `plot_from_file(**kwargs)`.
By specifying `plot_type`, you can switch between band, DOS, PDOS, or fatband modes.

```python
from qeplotter.qep import plot_from_file
# Ensure output directory exists
o import os
os.makedirs('outputs', exist_ok=True)

# 1) Band structure (plain black lines)
plot_from_file(
    plot_type   = 'band',
    band_file   = 'inputs/bso.bands.dat.gnu',   # QE bands output
    kpath_file  = 'inputs/bso.kpath',          # crystal_b K_POINTS
    savefig     = 'outputs/bso_band.png'
)

# 2) Total DOS
plot_from_file(
    plot_type   = 'dos',
    dos_file    = 'inputs/bso.dos',            # energy & DOS columns
    fermi_level = -5.3,                        # draw horizontal Fermi line
    shift_fermi = True,                        # shift energies so Fermi=0
    savefig     = 'outputs/bso_dos.png'
)

# 3) Projected DOS (atomic contributions)
plot_from_file(
    plot_type   = 'pdos',
    pdos_dir    = 'inputs/bso_pdos',           # folder of projwfc files
    pdos_mode   = 'atomic',                    # group by atom type
    savefig     = 'outputs/bso_pdos_atomic.png'
)

# 4) Fatband 'bubble' view: atomic
plot_from_file(
    plot_type      = 'fatbands',
    band_file      = 'inputs/bso.bands.dat.gnu',
    kpath_file     = 'inputs/bso.kpath',
    fatband_dir    = 'inputs/bso_pdos',
    fatbands_mode  = 'atomic',                  # bubble sized by orbital weight
    fermi_level    = 10.83,
    shift_fermi    = True,
    weight_threshold = 0.02,                    # show only >2% weights
    savefig        = 'outputs/bso_fatbands_atomic.png'
)
```

---

## 🔍 Parameter Reference

| Name                   | Type                 | Default    | Required for       | Description                                                                                                                                                                                                                                                                                                                         |         |          |                                      |
| ---------------------- | -------------------- | ---------- | ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- | -------- | ------------------------------------ |
| **plot\_type**         | `str`                | ―          | always             | `'band'`                                                                                                                                                                                                                                                                                                                            | `'dos'` | `'pdos'` | `'fatbands'`. Selects the plot mode. |
| **band\_file**         | `str`                | None       | `band`, `fatbands` | Path to QE `.bands.dat.gnu` file.                                                                                                                                                                                                                                                                                                   |         |          |                                      |
| **kpath\_file**        | `str`                | None       | `band`, `fatbands` | QE `K_POINTS crystal_b` file (defines k‑point path and labels).                                                                                                                                                                                                                                                                     |         |          |                                      |
| **dos\_file**          | `str`                | None       | `dos`, optional    | Two‑column DOS: energy \[eV], DOS. Required if `plot_type='dos'`; optional for fatbands DOS panel.                                                                                                                                                                                                                                  |         |          |                                      |
| **pdos\_dir**          | `str`                | None       | `pdos`             | Directory of `pdos_atm#...` files (projwfc output).                                                                                                                                                                                                                                                                                 |         |          |                                      |
| **fatband\_dir**       | `str`                | None       | `fatbands`         | Directory of fatband PDOS files (for projection coloring).                                                                                                                                                                                                                                                                          |         |          |                                      |
| **band\_mode**         | `str`                | `'normal'` | `band`             | Band coloring: `'normal'` (black), or `'atomic'`, `'orbital'`, `'element_orbital'`, `'most'`.                                                                                                                                                                                                                                       |         |          |                                      |
| **pdos\_mode**         | `str`                | `'atomic'` | `pdos`             | Grouping of PDOS: `'atomic'`, `'orbital'`, `'element_orbital'`.                                                                                                                                                                                                                                                                     |         |          |                                      |
| **fatbands\_mode**     | `str`                | `'most'`   | `fatbands`         | 🎨 Visualization mode:<br>• **Bubble**: `'most'`, `'atomic'`, `'orbital'`, `'element_orbital'`<br>• **Line**: `'normal'`, `'o_atomic'`, `'o_orbital'`, `'o_element_orbital'`<br>• **Heatmap**: `'heat_total'`, `'heat_atomic'`, `'heat_orbital'`, `'heat_element_orbital'`<br>• **Layer**: `'layer'` (requires `layer_assignment`). |         |          |                                      |
| **highlight\_channel** | `str` or `list`      | None       | line, heatmap      | Which channel to emphasize (e.g. `'Mo'`, `'d'`, `'Mo-d'`). Provide a list of two for dual color.                                                                                                                                                                                                                                    |         |          |                                      |
| **dual**               | `bool`               | `False`    | line               | If `True` with two channels, interpolates colors between them.                                                                                                                                                                                                                                                                      |         |          |                                      |
| **fermi\_level**       | `float`              | None       | any                | Draws a horizontal line at this energy; if `shift_fermi=True`, shifts all data by `-fermi_level`.                                                                                                                                                                                                                                   |         |          |                                      |
| **shift\_fermi**       | `bool`               | `False`    | any                | If `True`, subtracts `fermi_level` from all energies (bands & DOS) so Fermi sits at 0 eV.                                                                                                                                                                                                                                           |         |          |                                      |
| **y\_range**           | `tuple[float,float]` | None       | any                | Axis limits `(min, max)` for energy (bands) or DOS.                                                                                                                                                                                                                                                                                 |         |          |                                      |
| **cmap\_name**         | `str`                | `'tab10'`  | any                | Matplotlib colormap for all colored plots.                                                                                                                                                                                                                                                                                          |         |          |                                      |
| **s\_min**, **s\_max** | `float`              | `10`,`100` | bubble, heatmap    | Min & max marker sizes for bubble & heat modes.                                                                                                                                                                                                                                                                                     |         |          |                                      |
| **weight\_threshold**  | `float`              | `0.01`     | bubble, heatmap    | Minimum fraction of maximum weight to display.                                                                                                                                                                                                                                                                                      |         |          |                                      |
| **dpi**                | `int`                | None       | any                | Resolution of saved figures (dots per inch).                                                                                                                                                                                                                                                                                        |         |          |                                      |
| **layer\_assignment**  | `dict[str,str]`      | None       | `'layer'` mode     | Map atoms (`'Mo2'`,`'S3'`, etc.) to `'top'` or `'bottom'`.                                                                                                                                                                                                                                                                          |         |          |                                      |
| **savefig**            | `str`                | None       | any                | Path & filename for saving the figure. If omitted, just displays on screen.                                                                                                                                                                                                                                                         |         |          |                                      |

> **Tip:** you only need to supply the parameters relevant to your plot—others will be ignored.

---

## 📁 Project Structure & Examples

```
QEPlotter/
├── qeplotter/               # Core library files
│   ├── __init__.py
│   └── qep.py
│
├── saved/                   # Output PNGs generated by scripts
│   ├── bso_band.png
│   ├── bso_fatbands_atomic.png
│   ├── bso_heatmap_atomic.png
│   ├── ice_band.png
│   ├── ice_dos.png
│   ├── nsp_band.png
│   ├── nsp_band_atomic.png
│   └── ... (other generated figures)
│
├── requirements.txt         # Python dependencies
└── README.md                # Project description and usage


---

## 📄 License

MIT © shubics

