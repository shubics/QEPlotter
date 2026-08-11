# Project layout

QEPlotter 2.0 uses the modular package and `gui_mod.py` as its canonical
implementation. The former monolithic implementation is retained only as an
archive.

```text
QEPlotter/
├── gui_mod.py                  # Main Streamlit entry point
├── gui/                        # Streamlit pages and shared GUI helpers
├── qeplotter/                  # Maintained Python package
│   ├── api.py                  # plot_from_file() dispatcher
│   ├── core/                   # Parsers and shared utilities
│   ├── plotting/               # Bands, DOS/PDOS, and fatbands
│   ├── analysis/               # Band-gap and bilayer analysis
│   ├── structure/              # Structure parsing and 3D model support
│   ├── kpath/                  # Native high-symmetry path engine
│   ├── kmesh/                  # Irreducible k-grid reducer
│   ├── symmetry/               # Γ representations and SALCs
│   └── converters/             # Projection conversion utilities
├── docs/                       # Scientific methods and user documentation
├── examples/                   # Example inputs, outputs, and structures
├── tests/                      # Regression and scientific validation tests
├── archive/monolithic/         # Historical qep.py and legacy GUI
├── requirements.txt
├── pyproject.toml
└── setup.py
```

The modular implementation is exercised by the Python 3.9–3.12 CI matrix. The
archived monolith is preserved for older scripts but is not the source of new
features.
