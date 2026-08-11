# ⚛️ QEPlotter
> **Quantum ESPRESSO Band, DOS/PDOS, Fatband, and Analysis Toolkit**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![QEPlotter 2.0](https://img.shields.io/badge/QEPlotter-2.0-7893AE.svg)](https://github.com/shubics/QEPlotter)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Streamlit Cloud](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://qepweb.streamlit.app)
[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B?style=flat&logo=Streamlit&logoColor=white)](https://streamlit.io)

**QEPlotter** is a user-friendly Python library and GUI for post-processing Quantum ESPRESSO (QE) outputs.  
It provides a unified API to generate publication-ready plots (Bands, Fatbands, DOS) and robust tools for structure analysis and projection conversion.

---

## 🚀 Quick Start

### 1. Installation

```bash
git clone https://github.com/shubics/QEPlotter.git
cd QEPlotter
pip install -r requirements.txt
```

### 2. Run the GUI Dashboard

The easiest way to use QEPlotter is via the interactive web dashboard:

```bash
streamlit run gui_mod.py
```

> **Note:** `gui_app.py` is also available for legacy monolithic `qep.py` compatibility, but `gui_mod.py` powered by the `qeplotter` package is the recommended interface.

---

## ✨ Features at a Glance

### 🖥️ **Web Dashboard (`gui_mod.py`)**
- **Instant Visualization**: Drag & drop standard QE output files (`bands.dat.gnu`, `scf.out`, `pdos`).
- **High-Res Static Plots**: Publication-quality images via Matplotlib.
- **Auto-Fermi Detection**: Automatically reads Fermi energy from `scf.out` (supports metals, insulators, and semiconductors).
- **Band Gap Annotation**: Detects VBM/CBM and draws an arrow with the gap value directly on the plot.
- **Smart Analysis Tools**: Built-in band gap detector, bilayer structure analyzer, and projection converters.
- **Crystal Structure Explorer**: Browser-rendered 3D structures with symmetry, bonds, angles, bilayer recognition, and ordered six-site TMD stacking analysis.
- **Native K-path Engine**: Built-in Bravais classification, first-Brillouin-zone construction, and QE-ready paths without SeeK-path.
- **Irreducible K-grid Engine**: Exact integer-orbit reduction, weights, full-grid mapping, 3D IBZ inspection, and QE export from structure symmetry.
- **Symmetry Representations**: Γ-point reducible characters, numerical irreducible decomposition, SALCs, and symmetry-allowed `s/p/d` orbital matching.

### 📦 **Python API (Two Ways to Use)**

You have two options for using the Python API. **Both produce the exact same results with identical parameters.**

#### **Option 1: Modular Package ✅ (Recommended)**
The core is available as a clean, structured package (`qeplotter/`).
```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='fatbands',
    band_file='bands.dat.gnu',
    kpath_file='KPOINTS',
    fatband_dir='pdos_files/',
    fatbands_mode='atomic',
    spin=True
)
```

#### **Option 2: Standalone Script (`qep.py`)**
If you only need the legacy plotting and converter API in a portable,
single-file form, copy `qep.py` to your project directory.
```python
import qep

qep.plot_from_file(
    plot_type='fatbands',
    band_file='bands.dat.gnu',
    kpath_file='KPOINTS',
    fatband_dir='pdos_files/',
    fatbands_mode='atomic',
    spin=True
)
```

> **Note:** All examples below use the maintained modular package. The
> standalone script retains the plotting API, but the Structure Explorer,
> native K-path/K-grid engines, symmetry representations, and full v2 web
> interface are available only through `qeplotter` and `gui_mod.py`.

---

## 🖥️ Interface Preview

The views below were generated from the included
[2H monolayer MoS₂ structure](examples/structures/MoS2_monolayer.vasp).
QEPlotter identifies its `P-6m2` symmetry, constructs the conventional
high-symmetry band path, and reduces a `6 × 6 × 1` mesh to seven irreducible
k-points.

<p align="center">
  <img src="docs/images/mos2-structure-explorer.jpg" width="49%" alt="MoS2 monolayer in the QEPlotter crystal structure explorer">
  <img src="docs/images/mos2-kpath.jpg" width="49%" alt="MoS2 high-symmetry path and Brillouin zone in QEPlotter">
</p>

<p align="center">
  <img src="docs/images/mos2-irreducible-kgrid.jpg" width="760" alt="MoS2 irreducible 6 by 6 by 1 k-grid in QEPlotter">
</p>

## 📁 Project Structure

```
QEPlotter/
├── gui_mod.py              # 🖥  Main Streamlit GUI (uses modular package)
├── gui_app.py              # 🖥  Legacy GUI (uses monolithic qep.py)
│
├── qeplotter/              # 📦 Core modular Python package
│   ├── __init__.py         #     Package exports
│   ├── api.py              #     plot_from_file() dispatcher
│   ├── core/
│   │   ├── io.py           #     Parsers: K-path, band files, fatband files
│   │   └── utils.py        #     Helpers (strip_number, etc.)
│   ├── plotting/
│   │   ├── bands.py        #     Band structure + colored band modes
│   │   ├── dos.py          #     Total DOS + PDOS
│   │   └── fatbands.py     #     Bubble, Line, Heatmap, Layer fatbands
│   ├── analysis/
│   │   ├── bandgap.py      #     Band gap detection + plot annotation
│   │   └── bilayer.py      #     Structure analysis (stacking, interlayer distance)
│   ├── structure/           #     Multi-format structures, symmetry, bonds, browser 3D
│   ├── kpath/               #     Native K-path recipes + Wigner–Seitz BZ geometry
│   ├── kmesh/               #     Native uniform-grid → irreducible-BZ reducer
│   ├── symmetry/            #     Γ representations, character decomposition + SALCs
│   └── converters/
│       ├── fatbands.py     #     proj.out → pdos converter (consistent grid)
│       └── soc.py          #     SOC (j,mj) → (l,ml) basis converter
│
├── qep.py                  # 📄 Standalone monolithic script (all-in-one)
├── requirements.txt        # Python dependencies
├── pyproject.toml          # Modern Python build configuration
├── setup.py                # pip install support
└── examples/               # Sample QE outputs for testing
```

---

## 🔧 Core Capabilities & Technical Details

### 1. Band Structure Visualization
The library reads standard Quantum ESPRESSO `bands.x` output (`.gnu` format).

*   **Spin-Polarized Support**: When `spin=True`, fatband projection files with SOC/spin structure are parsed correctly. The `spin` flag is passed to the file reader for proper filename matching.
*   **Overlay Plotting**: The `overlay_band` mode superimposes two distinct band structures (e.g., *Bulk vs. Monolayer* or *PBE vs. HSE*). Handles disparate Fermi levels by aligning to 0 eV when `shift_fermi=True`.
*   **Band Gap Annotation**: When `show_band_gap=True`, automatically detects VBM/CBM and draws an arrow with the gap value. Supports SCF-based detection via `scf_file` for higher accuracy.

### 2. Fatbands (Projected Band Structure)
This is the library's most advanced feature, visualizing the contribution of specific atomic orbitals to the electronic bands. It parses filenames like `atm#1(Mo)_wfc#2(d)` to map weights to bands.

#### A. Visualization Modes
| Mode | Mechanism | Best For |
|:---|:---|:---|
| **Bubble** | Plots scatter markers where size $S \propto w_{orbital}$. Scaling is linear: $S = S_{min} + w(S_{max}-S_{min})$. | Quantitative analysis of contributions at specific k-points. |
| **Line** | Colors band segments based on weight using a `LineCollection`. | Clean, publication-quality plots where bubble clutter is undesirable. |
| **Heatmap** | Interpolates weights onto a dense grid using `imshow` (gaussian/bilinear interpolation). | Complex bands with high mixing (e.g., topological insulators). |

#### B. Grouping Logic
*   **`atomic`**: Sums all orbitals ($s, p, d, f$) belonging to a specific Atom ID (e.g., "Mo1").
*   **`orbital`**: Sums contributions by orbital type across *all* atoms (e.g., "Total d-character").
*   **`element_orbital`**: Granular breakdown (e.g., "Mo-d" vs "S-p").
*   **`sub_orb`**: If enabled, splits $d$ orbitals into $d_{z^2}, d_{x^2-y^2}, d_{xy}, d_{xz}, d_{yz}$. (Requires $m_j/m_l$ resolved data).

#### C. Layer-Resolved Plots
Useful for **Van der Waals heterostructures**.
*   **Input**: A dictionary mapping Atom Labels to Layer Names.
    *   *Example*: `{'Mo1': 'L1', 'S2': 'L1', 'W3': 'L2'}`
*   **Process**: The code aggregates weights for all atoms in "L1" vs "L2" and plots them with distinct colors (e.g., Blue for L1, Red for L2).
*   **Auto-Detection**: The GUI can automatically assign layers from a structure file using PBC-aware median splitting.

### 3. Density of States (DOS)
*   **Total DOS**: Reads standard `dos.x` output. Supports vertical orientation (Energy on Y-axis) for side-by-side comparison with bands.
*   **PDOS (Projected)**: Sums specific columns from `projwfc.x` output files grouped by `atomic`, `orbital`, or `element_orbital` mode.
*   **Side-by-Side DOS**: Band/fatband plots can include a Total DOS panel with `plot_total_dos=True`.

### 4. Advanced Analysis Tools

#### 🧭 Native K-path & Brillouin Zone
The modular GUI recommends band paths without a SeeK-path dependency. `spglib`
extracts the primitive cell and ASE supplies only the unimodular mapping to the
AFLOW canonical primitive basis. QEPlotter's internal Setyawan–Curtarolo registry
contains the named paths, extended Bravais cases, and lattice-dependent parameters
for all 14 three-dimensional Bravais lattices. The separate geometry engine builds
the first Brillouin zone as a reciprocal-lattice Wigner–Seitz cell using a Voronoi
construction; BZ vertices are never treated as band-path points. Both compact
`K_POINTS crystal_b` and explicit `K_POINTS crystal` cards can be downloaded.
Users can choose the full recommendation, its primary continuous branch, or
compose a custom route such as `Γ-X-M-Γ | R-X`; all exports and the 3D path update
to match the selected alternative.

#### Irreducible K-grid
The separate **K-grid & IBZ** page reduces uniform SCF/NSCF meshes to
irreducible Brillouin-zone representatives with multiplicities, normalized
weights, complete orbit mapping, and QE-ready output. spglib supplies only the
crystal's spatial symmetry operations; QEPlotter performs the shifted-grid
compatibility checks and exact integer-address orbit reduction itself.

See [Irreducible k-grid engine](docs/irreducible-k-grid.md) for the algorithm,
coordinate convention, magnetic/time-reversal considerations, and validation
invariants.

#### Symmetry & Orbital Representations
The dedicated representations page standardises an uploaded structure to its
primitive cell, identifies complete Wyckoff/site-symmetry orbits, and constructs
the transformation matrices for `s`, `p`, `d`, or atomic-displacement bases at
Γ. QEPlotter derives the finite point-group irreducible characters numerically,
so every detected symmetry class is included. It then decomposes the selected
reducible representation, generates normalised symmetry-adapted linear
combinations (SALCs), and lists `s/p/d` targets that share an irrep and are
therefore allowed to interact by symmetry. The result is structure-derived;
assigning irreps to individual electronic bands additionally requires QE
wavefunctions and is intentionally reported as a separate future capability.

Example: for the four H `1s` orbitals in the tetrahedral CH₄ validation
structure, the complete 24-operation `T_d` analysis gives:

$$
\Gamma_{\mathrm{H\,1s}}=A_1\oplus T_2 .
$$

The automated suite also covers diamond Si, rocksalt NaCl, hcp Mg,
wurtzite ZnO, 2H monolayer MoS₂, and the finite-group engine for all 32
crystallographic point groups.

See [Symmetry and Orbital Representations](docs/symmetry-representations.md)
for the equations, basis conventions, SALC construction, validation table, and
scientific limits.

#### 🔍 Band Gap Detector
Two detection methods with automatic fallback:

1.  **SCF-Based** (`scf_file`): Parses `scf.out` to extract HOMO/LUMO levels directly. Most accurate for insulators and semiconductors.
2.  **Fermi-Level Based** (fallback): Scans band energies relative to $E_F$ to find VBM/CBM.

Detection output:
*   **VBM/CBM positions** (energy, k-point index)
*   **Gap type**: Direct ($k_{VBM} = k_{CBM}$) or Indirect ($k_{VBM} \ne k_{CBM}$)
*   **Plot annotation**: Arrow from VBM→CBM with $E_g$ label

#### 🧱 Structure Analyzer (`analyse_file`)
Parses `scf.in` or `scf.out` to determine 2D material properties:
*   **Stacking Order**: Safe detection of R-type AA/AB/BA and H-type AA′/AB′/A′B registries for commensurate trigonal-prismatic bilayers. Janus/heterobilayer interfaces remain explicit; unsupported or twisted structures are reported as `General registry`.
*   **Interlayer Distance**: Calculates the vertical ($\Delta z$) distance between defined layers.

Method, naming convention, applicability, and references:
[Bilayer stacking analysis](docs/stacking-analysis.md).

#### 🔄 Consistency Converters
Quantum ESPRESSO's `projwfc.x` has known quirks.
*   **`convert_consistent`**: Fixes the "varying rows" issue where projections with zero weight are omitted, breaking standard plotters. Forces a dense rectangular data grid.
*   **`convert_soc_proj_to_ml`**: For **SOC** runs. Converts from the $(j, m_j)$ coupled basis to the standard $(l, m_l)$ basis (e.g., $d_{xy}, d_{z^2}$) using Clebsch-Gordan coefficients.

---

## 📁 Supported Input Files

| File | Description | Source |
|------|-------------|--------|
| `*.bands.dat.gnu` | Band structure data | `bands.x` |
| `*.kpath` / `*.in` | High-symmetry k-points path (crystal_b) | Custom |
| `scf.out` | Self-consistent field output | `pw.x` |
| `proj.out` | Projections log | `projwfc.x` |
| `*_pdos*` | Projected DOS / fatband files | `projwfc.x` |
| `*.dos` | Total DOS file | `dos.x` |

---

## 🎛️ Complete Parameter Reference

The `plot_from_file` function is the main entry point. Below is the **exhaustive list** of all available arguments.

### Core Parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`plot_type`** | `str` | `'band'` | **Required**. Options: `'band'`, `'dos'`, `'pdos'`, `'fatbands'`, `'overlay_band'`. |
| **`band_file`** | `str` | `None` | Path to QE band file (`.gnu`). Required for Band/Fatbands. |
| **`kpath_file`** | `str` | `None` | Path to K-path file (`crystal_b`). Required for Band/Fatbands. |
| **`fermi_level`** | `float` | `None` | Fermi energy (eV). Essential for energy alignment. |
| **`shift_fermi`** | `bool` | `False` | Moves Fermi level to 0 eV if `fermi_level` is provided. |
| **`y_range`** | `tuple` | `None` | Min/Max Y-axis limits. Example: `(-5, 5)`. |
| **`dpi`** | `int` | `None` | Resolution for saved figures. |
| **`save_dir`** | `str` | `'saved'` | Directory to save output images. |
| **`savefig`** | `str` | `None` | Filename to save plot (e.g., `'fig1.png'`). If `None`, shows window. |

### Band & Fatband Parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`band_mode`** | `str` | `'normal'` | Band coloring: `'normal'`, `'atomic'`, `'orbital'`, `'element_orbital'`, `'most'`. |
| **`fatband_dir`** | `str` | `None` | Folder with `projwfc.x` outputs. Required for Fatbands. |
| **`fatbands_mode`** | `str` | `'most'` | Fatband style: `'most'`, `'atomic'`, `'orbital'`, `'heat_*'`, `'layer'`, `'o_*'`. |
| **`spin`** | `bool` | `False` | Enable spin-polarized plotting. |
| **`sub_orb`** | `bool` | `False` | Split orbitals into sub-components ($d_{z^2}$, etc.). |
| **`highlight_channel`** | `str/tuple` | `None` | Atoms/orbitals to highlight. Ex: `'Mo'` or `('Mo', 'S')`. |
| **`dual`** | `bool` | `False` | Two-channel comparison in Line mode. |
| **`layer_assignment`** | `dict` | `None` | Layer Mode: Map atoms to layers. Ex: `{'Mo1': 'top', 'S2': 'bottom'}`. |
| **`s_min`** | `float` | `10` | Minimum marker size (bubble plots). |
| **`s_max`** | `float` | `100` | Maximum marker size (bubble plots). |
| **`weight_threshold`** | `float` | `0.01` | Ignore orbital weights below this value. |
| **`cmap_name`** | `str` | `'tab10'` | Matplotlib colormap name. |
| **`show_band_gap`** | `bool` | `False` | Detect and annotate VBM/CBM band gap on the plot. |
| **`scf_file`** | `str` | `None` | Path to `scf.out` for accurate SCF-based gap detection. |

### DOS Parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`dos_file`** | `str` | `None` | Path to `dos.x` output. Required for DOS plots. |
| **`pdos_dir`** | `str` | `None` | Folder with `pdos` files. Required for PDOS plots. |
| **`pdos_mode`** | `str` | `'atomic'` | PDOS grouping: `'atomic'`, `'orbital'`, `'element_orbital'`. |
| **`plot_total_dos`** | `bool` | `False` | Plot Total DOS side panel with Band/Fatband plots. |
| **`x_range`** | `tuple` | `None` | DOS axis range for side panel. |
| **`vertical`** | `bool` | `False` | Rotate DOS plot (Energy on Y-axis). |

### Heatmap Parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`overlay_bands_in_heat`** | `bool` | `False` | Draw band lines on top of heatmap. |
| **`heat_vmin`** | `float` | `None` | Minimum for heatmap color normalization. |
| **`heat_vmax`** | `float` | `None` | Maximum for heatmap color normalization. |

### Overlay Parameters

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`band_file2`** | `str` | `None` | Path to second band file. |
| **`kpath_file2`** | `str` | `None` | Path to second K-path file. |
| **`label1`** | `str` | `'Band 1'` | Legend label for first band. |
| **`label2`** | `str` | `'Band 2'` | Legend label for second band. |
| **`color1`** | `str` | `'red'` | Color for first band. |
| **`color2`** | `str` | `'blue'` | Color for second band. |

---

## 🖼️ Plot Examples

### 1. Basic Band Structure

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='band',
    band_file='molybden.bands.dat.gnu',
    kpath_file='molybden.kpath',
    fermi_level=6.234,
    shift_fermi=True,
    y_range=(-4, 6),
    band_mode='normal',
    spin=True
)
```

<p align="center">
  <img src="examples/example_outputs/BSO_Monolayer_spin_band.png" width="600" alt="Spin-polarized band structure">
</p>

### 2. Band Structure with Band Gap Annotation

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='band',
    band_file='Bi2Se3.bands.dat.gnu',
    kpath_file='Bi2Se3.kpath',
    fermi_level=-4.22,
    shift_fermi=True,
    y_range=(-3, 3),
    show_band_gap=True,                  # Annotate VBM/CBM with arrow
    scf_file='scf.out'                   # Use SCF data for accurate gap
)
```
<p align="center">
  <img src="examples/example_outputs/plot-4.png" width="600" alt="band plot">
</p>
<p align="center">
  <img src="examples/example_outputs/plot-5.png" width="600" alt="plotbyside">
</p>


### 3. Fatbands: Atomic Contribution (Bubble Plot)

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='fatbands',
    band_file='molybden.bands.dat.gnu',
    kpath_file='molybden.kpath',
    fatband_dir='./pdos_files',
    fatbands_mode='atomic',
    fermi_level=6.234,
    shift_fermi=True,
    s_min=1, s_max=150,
    spin=True
)
```

<p align="center">
  <img src="examples/example_outputs/fatbands_atomic.jpeg" width="600" alt="Fatbands atomic bubble plot">
</p>

### 4. Fatbands: Element Heatmap

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='fatbands',
    band_file='Bi2Se3.bands.dat.gnu',
    kpath_file='Bi2Se3.kpath',
    fatband_dir='./Bi2Se3_pdos',
    fatbands_mode='heat_orbital',
    highlight_channel='p',
    fermi_level=4.11,
    shift_fermi=True,
    cmap_name='inferno',
    heat_vmax=10.0,
    overlay_bands_in_heat=True
)
```

<p align="center">
  <img src="examples/example_outputs/fatbands_heat_atomic.jpeg" width="600" alt="Fatbands heatmap">
</p>

### 5. Element-Orbital Decomposition

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='fatbands',
    band_file='MoS2.bands.dat.gnu',
    kpath_file='MoS2.kpath',
    fatband_dir='./pdos',
    fatbands_mode='element_orbital',
    fermi_level=1.5,
    weight_threshold=0.1,
    spin=False
)
```



### 6. Overlay Band Comparison

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='overlay_band',
    band_file='MoS2_ML.bands.gnu',
    kpath_file='MoS2.kpath',
    label1='Monolayer', color1='red',
    band_file2='MoS2_Bulk.bands.gnu',
    kpath_file2='MoS2.kpath',
    label2='Bulk', color2='blue',
    fermi_level=5.0,
    shift_fermi=True
)
```

<p align="center">
  <img src="examples/example_outputs/BandStructure_MoS2_Monolayer_Spin_vs_MoS2_Monolayer_Non_Spin.png" width="600" alt="Overlay band comparison">
</p>

### 7. Layer-Resolved Fatbands

```python
from qeplotter import plot_from_file

layer_map = {
    'Mo1': 'top', 'S2': 'top', 'S3': 'top',
    'W4':  'bottom', 'Se5': 'bottom', 'Se6': 'bottom'
}

plot_from_file(
    plot_type='fatbands',
    band_file='hetero.bands.gnu',
    kpath_file='hetero.kpath',
    fatband_dir='./pdos',
    fatbands_mode='layer',
    layer_assignment=layer_map,
    fermi_level=3.5
)
```

<p align="center">
  <img src="examples/example_outputs/fatbands_layer.jpeg" width="600" alt="Layer-resolved fatbands">
</p>

### 8. Total DOS

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='dos',
    dos_file='system.dos',
    fermi_level=5.4,
    shift_fermi=True,
    vertical=True,
    savefig='dos_plot.pdf'
)
```

### 9. Projected DOS (PDOS)

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='pdos',
    pdos_dir='./pdos',
    pdos_mode='element_orbital',
    fermi_level=5.4,
    shift_fermi=True
)
```

### 10. Band + DOS Side-by-Side with Gap

```python
from qeplotter import plot_from_file

plot_from_file(
    plot_type='band',
    band_file='data.gnu',
    kpath_file='path.kpath',
    fermi_level=0.0,
    shift_fermi=True,
    plot_total_dos=True,
    dos_file='system.dos',
    show_band_gap=True,
    scf_file='scf.out',
    dpi=300,
    savefig='figure_1.png'
)
```

---

## 📜 License

MIT © shubics
