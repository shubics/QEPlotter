# QEPlotter

**Quantum ESPRESSO Band, DOS/PDOS, Fatband, and Analysis Toolkit**

QEPlotter is a Python library for post-processing and visualization of Quantum ESPRESSO (QE) outputs.  
With a unified high-level API, you can generate publication-ready plots and also perform structure analysis, projection conversion, and band gap detection.

---

## 🚀 Installation

```bash
pip install git+https://github.com/shubics/QEPlotter.git
```

Or work from the source:

```bash
git clone https://github.com/shubics/QEPlotter.git
cd QEPlotter
pip install -e .
```

**Requirements:**
- Python 3.8+
- `numpy`, `matplotlib`  
  *Optional:* TeX Live (LaTeX labels), PDF/Matlab/Mathematica export tools

---

## 🔧 Key Features

### 🎯 Band Structure Visualization
- **Modes**:
  - `normal`: classic band plot
  - `atomic`, `orbital`, `element_orbital`: colored bands based on contribution
  - `overlay_band`: compare two band structures

### 🎨 Fatbands Visualization
- **Bubble Modes**:
  - `atomic`, `orbital`, `element_orbital`
- **Line Modes**:
  - `o_atomic`, `o_orbital`, `o_element_orbital`
- **Heatmap Modes**:
  - `heat_total`, `heat_atomic`, `heat_orbital`, `heat_element_orbital`
- **Layer Mode**:
  - Visualize interlayer contributions in bilayers

### 📈 DOS/PDOS Plots
- Total and projected DOS plots

### 🧰 Side Tools
- `detect_band_gap(band_file, kpt_file, fermi=None)`: Automatically detects direct or indirect band gaps and reports VBM/CBM.
- `analyse_file(path)`: Analyzes bilayer stacking and calculates vertical/3D interatomic distances.
- `convert_soc_proj_to_ml(...)`: Converts `proj.out` files into QEPlotter-compatible fatband files with atomic/orbital resolution.

---

## 📁 Input Files
| File Type | Description |
|----------|-------------|
| `*.bands.dat.gnu` | Band structure file (from `bands.x`) |
| `*.kpath` | High-symmetry k-point file |
| `proj.out` | Output from `projwfc.x` |
| `*_pdos` | Directory of projected DOS files |
| `*.dos` | Total DOS file from `dos.x` |

---
## 🎛️ Key Parameters

| Name                     | Type                 | Applies to         | Description |
|--------------------------|----------------------|--------------------|-------------|
| **plot_type**            | `str`                | all                | `'band'`, `'dos'`, `'pdos'`, `'fatbands'` |
| **band_file**            | `str`                | band, fatbands     | QE `.bands.dat.gnu` file |
| **kpath_file**           | `str`                | band, fatbands     | QE `K_POINTS crystal_b` file |
| **dos_file**             | `str`                | dos (+fatbands)    | Two-column file with energy and DOS |
| **pdos_dir**             | `str`                | pdos               | Folder with `pdos_atm#` files |
| **fatband_dir**          | `str`                | fatbands           | Folder with fatband projection files |
| **band_mode**            | `str`                | band               | `'normal'`, `'atomic'`, `'orbital'`, `'element_orbital'`, `'most'` |
| **pdos_mode**            | `str`                | pdos               | `'atomic'`, `'orbital'`, `'element_orbital'` |
| **fatbands_mode**        | `str`                | fatbands           | Various bubble, line, heatmap, or layer modes |
| **highlight_channel**    | `str` / `list`       | line, heatmap      | Emphasize specific atom/orbital/channel |
| **dual**                 | `bool`               | line               | Blend colors between two highlight channels |
| **spin**                 | `bool`               | band, fatbands     | Plot spin-resolved results |
| **sub_orb**              | `bool`               | band, fatbands     | Show sub-orbital separation |
| **fermi_level**          | `float`              | all                | Reference Fermi energy |
| **shift_fermi**          | `bool`               | all                | Shift Fermi level to 0 eV |
| **y_range**              | `(float,float)`      | all                | Set y-axis limits |
| **cmap_name**            | `str`                | all                | Matplotlib colormap |
| **s_min**, **s_max**     | `float`              | bubble/heatmap     | Marker size range |
| **weight_threshold**     | `float`              | bubble/heatmap     | Ignore small weights |
| **dpi**                  | `int`                | all                | Image resolution |
| **layer_assignment**     | `dict[str,str]`      | layer              | Atom-to-layer mapping |
| **plot_total_dos**       | `bool`               | fatbands           | Show total DOS beside fatbands |
| **overlay_bands_in_heat**| `bool`               | heatmap            | Overlay band lines on heatmap |
| **heat_vmin**, **heat_vmax** | `float`          | heatmap            | Color scale normalization |
| **savefig**              | `str`                | all                | Output path (omit to show interactively) |

---



## 🖼️ Plot Examples

### 1. **Normal Band Structure**
```python
plot_from_file(
    plot_type='band',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fermi_level=10.7431,
    shift_fermi=True,
    y_range=(-3, 3),
    band_mode='normal',
    spin=True
)
```

### 2. **Fatbands - Atomic Bubble Mode**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='Bandx.dat.gnu',
    kpath_file='proj.kpath',
    fatband_dir='BBN_pdos',
    fatbands_mode='atomic',
    fermi_level=-5.3321,
    shift_fermi=True,
    y_range=(-3,3),
    cmap_name='tab10',
    spin=True
)
```

### 3. **Element-Orbital Bubble Mode**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fatband_dir='BBS_pdos',
    fatbands_mode='element_orbital',
    fermi_level=10.7431,
    shift_fermi=True,
    cmap_name='magma',
    weight_threshold=0.05,
    s_max=300,
    s_min=2,
    spin=True
)
```

### 4. **O-Orbital Line Mode**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fatband_dir='BBS_pdos',
    fatbands_mode='o_orbital',
    fermi_level=10.7431,
    shift_fermi=True,
    dual=False,
    highlight_channel=('p',),
    spin=True
)
```

### 5. **Dual O-Element-Orbital Line Mode**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fatband_dir='BBS_pdos',
    fatbands_mode='o_element_orbital',
    fermi_level=10.7431,
    shift_fermi=True,
    dual=True,
    highlight_channel=('Bi1-p','O8-p'),
    spin=True
)
```

### 6. **Atomic Colored Band Mode**
```python
plot_from_file(
    plot_type='band',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fatband_dir='BBS_pdos',
    band_mode='atomic',
    spin=True
)
```

### 7. **Overlay Band Mode**
```python
plot_from_file(
    plot_type='overlay_band',
    band_file='MMS.bands.dat.gnu',
    kpath_file='MMS.kpath',
    band_file2='MMN.bands.dat.gnu',
    kpath_file2='MMN.kpath',
    label1='MoS2 Bulk Spin',
    label2='MoS2 Bulk Non-Spin',
    shift_fermi=True,
    color1='blue',
    color2='red'
)
```

### 8. **Layer Fatbands Mode**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='Bandx.dat.gnu',
    kpath_file='proj.kpath',
    fatband_dir='BBN_pdos',
    fatbands_mode='layer',
    layer_assignment={'Mo2':'bottom', 'S3':'top', 'S4':'bottom', 'S5':'top'},
    shift_fermi=True,
    spin=True
)
```

### 9. **Heatmap Mode - heat_atomic**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='BBS.bands.dat.gnu',
    kpath_file='BBS.kpath',
    fatband_dir='BBS_pdos',
    fatbands_mode='heat_atomic',
    shift_fermi=True,
    heat_vmax=40,
    heat_vmin=0,
    highlight_channel=('Bi',),
    spin=True
)
```

### 10. **Sub-Orbital Decomposition**
```python
plot_from_file(
    plot_type='fatbands',
    band_file='W-W.AA.bands.dat.gnu',
    kpath_file='W-W.AA.kpath',
    fatband_dir='W-W.AA_pdos',
    fatbands_mode='element_orbital',
    sub_orb=True,
    spin=False
)
```

---

## 🛠 Side Tools

### 🔍 `detect_band_gap(band_file, kpt_file, fermi=None)`
- Computes direct/indirect band gap
- Locates VBM and CBM high-symmetry points

### 🧪 `analyse_file(path)`
- Automatically detects stacking type (AA, AB, AA′)
- Calculates all vertical and 3D distances from QE `ibrav` input

### 🧰 `convert_soc_proj_to_ml(proj_out_path, ...)`
- Converts `proj.out` to QEPlotter-compatible `.pdos` files
- Supports j-split orbitals
- Fully supports SOC and spin

---

## ✅ Advantages

- ⚡ **Automated**: No manual editing of QE files
- 🎯 **Precise**: Detects and visualizes fine structural details
- 🌍 **Flexible**: Works for bulk, 2D, and heterostructures
- 🧩 **Modular**: Easy integration with ML pipelines

---
## 📁 Project Structure

```
QEPlotter/
├── qeplotter/
│   ├── __init__.py
│  
│   └── qep.py
├── saved/
├── requirements.txt
└── README.md
```

---

## 📜 License

MIT © shubics
