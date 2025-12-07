# ⚛️ QEPlotter
> **Quantum ESPRESSO Band, DOS/PDOS, Fatband, and Analysis Toolkit**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
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
streamlit run qui_app.py
```

---

## � Features at a Glance

### 🖥️ **QUI-Web Dashboard**
- **Instant Visualization**: Drag & drop standard QE output files (`bands.dat.gnu`, `scf.out`, `pdos`).
- **High-Res Static Plots**: Generates publication-quality images using Matplotlib.
- **Auto-Fermi**: Automatic detection of Fermi energy from `scf.out`.
- **Smart Analysis**: Built-in band gap detector and bilayer structure analyzer.

### � **Python Library (`qep.py`)**

If you prefer scripting, copy `qep.py` to your project and use it directly:

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

---

## 📁 Project Structure

The project uses a simple flat structure for easy usage:

```
QEPlotter/
├── qui_app.py         # Main Streamlit GUI Application
├── qep.py             # Core plotting & analysis library
├── requirements.txt   # Python dependencies
├── README.md          # Documentation
└── example_data/      # (Optional) Sample QE outputs
```

---

## 🔧 Core Capabilities & Technical Details

### 1. Band Structure Visualization
The library reads standard Quantum ESPRESSO `bands.x` output (`.gnu` format).

*   **Spin-Polarized Support**: Automatically detects spin channels. If `spin=True`, it plots Spin Up (red/blue) and Spin Down (blue/red) dashed lines.
*   **Overlay Plotting**: The `overlay_band` mode allows you to superimpose two distinct band structures (e.g., *Bulk vs. Monolayer* or *PBE vs. HSE*). It handles disparate Fermi levels by aligning them to 0 eV automatically if `shift_fermi=True`.

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

### 3. Density of States (DOS)
*   **Total DOS**: reads standard `dos.x` output. Supports automatic rotation of the plot (Energy on Y-axis) for side-by-side comparison with bands.
*   **PDOS (Projected)**: Sums specific columns from `projwfc.x` output files. It validates that the energy grid matches the Total DOS file (if provided) to ensure alignment.

### 4. Advanced Analysis Tools

#### 🔍 Band Gap Detector (`detect_band_gap`)
Calculates the electronic band gap by scanning the eigenenergies:
1.  **VBM Identification**: Finds the highest occupied state (energy $\le E_F$).
2.  **CBM Identification**: Finds the lowest unoccupied state (energy $> E_F$).
3.  **Nature Detection**: Compares the k-point indices of VBM and CBM.
    *   If $k_{VBM} == k_{CBM} \rightarrow$ **Direct Gap**.
    *   If $k_{VBM} \ne k_{CBM} \rightarrow$ **Indirect Gap**.

#### 🧱 Structure Analyzer (`analyse_file`)
Parses `scf.in` or `scf.out` to determine 2D material properties:
*   **Stacking Order**: Heuristic detection of AA, AB, or AA' stacking in bilayers based on relative coordinates.
*   **Interlayer Distance**: Calculates the vertical ($\Delta z$) distance between defined layers.

#### 🔄 Consistency Converters
Quantum ESPRESSO's `projwfc.x` has known quirks.
*   **`convert_consistent`**: Fixes the "varying rows" issue where projections with zero weight are sometimes omitted by QE, breaking standard plotters. It forces a dense rectangular data grid.
*   **`convert_soc_proj_to_ml`**: For **Spin-Orbit Coupling (SOC)** runs. QE outputs projections in the $(j, m_j)$ coupled basis. This tool uses Clebsch-Gordan coefficients to uncouple them back to the standard $(l, m_l)$ basis (e.g., $d_{xy}, d_{z^2}$) familiar to chemists.

---
---

## 📁 Supported Input Files

| File | Description | Source |
|------|-------------|--------|
| `*.bands.dat.gnu` | Band structure data | `bands.x` |
| `*.kpath` | High-symmetry k-points path | Custom / `qep` format |
| `scf.out` | Self-consistent field output | `pw.x` |
| `proj.out` | Projections log | `projwfc.x` |
| `*_pdos*` | Projected DOS files | `projwfc.x` |
| `*.dos` | Total DOS file | `dos.x` |

---

## 🎛️ Complete Parameter Reference

The `plot_from_file` function is the main entry point. Below is the **exhaustive list** of all 33 available arguments.

| Parameter | Type | Default | Description |
|:---|:---|:---|:---|
| **`plot_type`** | `str` | `'band'` | **Required**. Options: `'band'`, `'dos'`, `'pdos'`, `'fatbands'`, `'overlay_band'`. |
| **`band_file`** | `str` | `None` | Path to QE band file (`.gnu`). Required for Band/Fatbands. |
| **`kpath_file`** | `str` | `None` | Path to K-path file (`crystal_b`). Required for Band/Fatbands. |
| **`fermi_level`** | `float` | `None` | Fermi energy (eV). Essential for energy alignment. |
| **`shift_fermi`** | `bool` | `False` | Moves Fermi level to 0 eV if `fermi_level` is provided. |
| **`fatband_dir`** | `str` | `None` | Folder with `projwfc.x` outputs. Required for Fatbands. |
| **`pdos_dir`** | `str` | `None` | Folder with `pdos` files. Required for PDOS plots. |
| **`dos_file`** | `str` | `None` | Path to `dos.x` output file. Required for total DOS plots. |
| **`y_range`** | `tuple` | `None` | Min/Max Y-axis limits. Example: `(-5, 5)`. |
| **`spin`** | `bool` | `False` | Enable spin-polarized plotting (Up/Down channels). |
| **`sub_orb`** | `bool` | `False` | Split orbitals into sub-components ($d_{z^2}$ etc.) instead of total ($l$). |
| **`band_mode`** | `str` | `'normal'` | Band plot style: `'normal'`, `'atomic'`, `'orbital'`, `'element_orbital'`. |
| **`fatbands_mode`** | `str` | `'most'` | Fatband style: `'atomic'`, `'orbital'`, `'element_orbital'`, `'heat_*'`, `'layer'`. |
| **`pdos_mode`** | `str` | `'atomic'` | PDOS grouping: `'atomic'`, `'orbital'`, `'element_orbital'`. |
| **`highlight_channel`** | `str/tuple` | `None` | Atoms/orbitals to highlight in Line/Heatmap modes. Ex: `('Mo',)` |
| **`dual`** | `bool` | `False` | Compare two highlight channels in Line mode. |
| **`layer_assignment`** | `dict` | `None` | **Layer Mode**: Map atoms to layers. Ex: `{'Mo1': 'L1', 'S2': 'L2'}`. |
| **`s_min`** | `float` | `10` | Minimum marker size for bubble plots. |
| **`s_max`** | `float` | `100` | Maximum marker size for bubble plots. |
| **`weight_threshold`** | `float` | `0.01` | Ignore orbital weights below this value to reduce clutter. |
| **`cmap_name`** | `str` | `'tab10'` | Matplotlib colormap name (e.g., `'viridis'`, `'magma'`). |
| **`overlay_bands_in_heat`** | `bool` | `False` | Draw thin gray band lines on top of heatmap fatbands. |
| **`heat_vmin`** | `float` | `None` | Minimum value for heatmap color normalization. |
| **`heat_vmax`** | `float` | `None` | Maximum value for heatmap color normalization. |
| **`plot_total_dos`** | `bool` | `False` | Plot Total DOS side-by-side with Fatbands. |
| **`band_file2`** | `str` | `None` | **Overlay Mode**: Path to second band file. |
| **`kpath_file2`** | `str` | `None` | **Overlay Mode**: Path to second K-path file. |
| **`label1`** | `str` | `'Band 1'` | **Overlay Mode**: Legend label for first band. |
| **`label2`** | `str` | `'Band 2'` | **Overlay Mode**: Legend label for second band. |
| **`color1`** | `str` | `'red'` | **Overlay Mode**: Color for first band. |
| **`color2`** | `str` | `'blue'` | **Overlay Mode**: Color for second band. |
| **`save_dir`** | `str` | `'saved'` | Directory to save output images. |
| **`savefig`** | `str` | `None` | Filename to save plot (e.g., `'fig1.png'`). If `None`, shows window. |
| **`dpi`** | `int` | `None` | Resolution for saved figures (default: matplotlibrc). |

---

## 🖼️ Detailed Plot Examples

Here are comprehensive examples for various scenarios. You can run these in a Python script (e.g., `run_plot.py`) by importing `qep`.

### 1. Basic Band Structure (Spin-Polarized)
*Standard band plot with spin up/down.*

```python
import qep

qep.plot_from_file(
    plot_type='band',
    band_file='molybden.bands.dat.gnu',  # Your QE bands file
    kpath_file='molybden.kpath',         # Your K-points path
    fermi_level=6.234,                   # E_F from scf.out
    shift_fermi=True,                    # Shift so 0 eV is Fermi level
    y_range=(-4, 6),                     # Show 4eV below and 6eV above Fermi
    band_mode='normal',                  # Simple black lines
    spin=True                            # Enable spin-polarized plotting
)
```

### 2. Fatbands: Atomic Contribution (Bubble Plot)
*Bubble size indicates how much a specific atom contributes to the band.*

```python
qep.plot_from_file(
    plot_type='fatbands',
    band_file='molybden.bands.dat.gnu',
    kpath_file='molybden.kpath',
    fatband_dir='./pdos_files',          # Folder with projwfc.x outputs
    fatbands_mode='atomic',              # Color bubbles by Atom ID (Mo1, S2...)
    fermi_level=6.234,
    shift_fermi=True,
    s_min=1, s_max=150,                  # Min/Max size of bubbles
    spin=True
)
```

### 3. Fatbands: Orbital Heatmap
*Best for complex systems. Heatmap intensity shows orbital character (e.g., d-orbitals).*

```python
qep.plot_from_file(
    plot_type='fatbands',
    band_file='Bi2Se3.bands.dat.gnu',
    kpath_file='Bi2Se3.kpath',
    fatband_dir='./Bi2Se3_pdos', 
    fatbands_mode='heat_orbital',        # Heatmap based on s/p/d/f orbitals
    highlight_channel=('p',),            # Highlight 'p' orbitals specifically
    fermi_level=4.11,
    shift_fermi=True,
    cmap_name='inferno',                 # 'magma', 'viridis', 'plasma', etc.
    heat_vmax=10.0,                      # Contrast adjustment for heatmap
    overlay_bands_in_heat=True           # Draw thin band lines over the heatmap
)
```

### 4. Element-Orbital Decomposition (Detailed)
*Distinguish between Mo-d states and S-p states with different colors.*

```python
qep.plot_from_file(
    plot_type='fatbands',
    band_file='MoS2.bands.dat.gnu',
    kpath_file='MoS2.kpath',
    fatband_dir='./pdos',
    fatbands_mode='element_orbital',     # e.g., Mo-s, Mo-p, Mo-d, S-s, S-p...
    fermi_level=1.5,
    weight_threshold=0.1,                # Ignore contributions < 10% to clean up plot
    spin=False
)
```

### 5. Overlay Band Structure (Comparison)
*Compare two different calculations (e.g., Bulk vs Monolayer, or Theory vs Theory).*

```python
qep.plot_from_file(
    plot_type='overlay_band',
    # --- System 1 (e.g., Monolayer) ---
    band_file='MoS2_ML.bands.gnu',
    kpath_file='MoS2.kpath',
    label1='Monolayer',
    color1='red',
    
    # --- System 2 (e.g., Bulk) ---
    band_file2='MoS2_Bulk.bands.gnu',
    kpath_file2='MoS2.kpath',
    label2='Bulk',
    color2='blue',
    
    fermi_level=5.0,
    shift_fermi=True
)
```

### 6. Sub-Orbital Analysis ($d_{z^2}$, $d_{xy}$, etc.)
*Deep dive into crystal field splitting.*

```python
qep.plot_from_file(
    plot_type='fatbands',
    band_file='oxide.bands.gnu',
    kpath_file='oxide.kpath',
    fatband_dir='./pdos',
    fatbands_mode='orbital',
    sub_orb=True,                        # <--- Enables splitting d -> dz2, dxy...
    fermi_level=7.0,
    shift_fermi=True
)
```

### 7. Layer-Resolved Fatbands
*For 2D heterostructures: See which layer contributes to which band.*

```python
# Define which atoms belong to which layer
layer_map = {
    'Mo1': 'Layer1', 'S2': 'Layer1', 'S3': 'Layer1',
    'W4':  'Layer2', 'Se5': 'Layer2', 'Se6': 'Layer2'
}

qep.plot_from_file(
    plot_type='fatbands',
    band_file='hetero.bands.gnu',
    kpath_file='hetero.kpath',
    fatband_dir='./pdos',
    fatbands_mode='layer',               # Activate Layer mode
    layer_assignment=layer_map,          # Pass the mapping
    fermi_level=3.5
)
```

### 8. Total DOS Plot
*Simple Density of States.*

```python
qep.plot_from_file(
    plot_type='dos',
    dos_file='system.dos',               # Output from dos.x
    fermi_level=5.4,
    y_range=(-10, 10),                   # DOS value range (x-axis typically)
    shift_fermi=True,
    savefig='dos_plot.pdf'               # Save as PDF
)
```

### 9. Projected DOS (PDOS)
*DOS resolved by atom or orbital.*

```python
qep.plot_from_file(
    plot_type='pdos',
    fatband_dir='./pdos',                # Directory with projwfc.x files
    pdos_mode='element_orbital',         # Group by Element+Orbital (Mo-d, S-p)
    fermi_level=5.4,
    shift_fermi=True
)
```

### 10. High-Quality Customization
*Control resolution and saving.*

```python
qep.plot_from_file(
    plot_type='band',
    band_file='data.gnu',
    kpath_file='path.kpath',
    dpi=300,                             # High resolution for print
    savefig='figure_1.tiff',             # TIFF format
    fermi_level=0.0,
    font_family='serif'                  # (Optional internal customization)
)
```

---

## 📜 License

MIT © shubics
