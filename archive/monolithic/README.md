# Archived monolithic interface

This directory preserves QEPlotter's former all-in-one implementation:

- `qep.py`: standalone plotting and converter module;
- `gui_app.py`: legacy Streamlit interface built on `qep.py`.

These files are archived for reproducibility and migration of older scripts.
They are not the canonical QEPlotter 2.0 implementation and do not receive the
same feature coverage as the modular `qeplotter/` package and `gui_mod.py`.
The standalone plotting API aims to remain compatible where practical, but new
work should use:

```python
from qeplotter import plot_from_file
```

To inspect the historical GUI from the repository root:

```bash
python -m streamlit run archive/monolithic/gui_app.py
```

No files in this archive are imported by the main application.
