# WS₂/MoS₂ layer-projection example

This directory contains historical Quantum ESPRESSO-format band and
`projwfc.x` projection files that have been bundled with QEPlotter since the
repository's first clean commit (2025-07-01). They are useful for exercising
band, fatband, and layer-colour plotting.

The original SCF/NSCF inputs, pseudopotentials, Quantum ESPRESSO version,
convergence settings, and calculation record are not present. Consequently,
these files must be treated as **unverified example data**, not as a validated
or reproducible scientific result.

The layer assignment preserved by the original notebook is:

- lower WS₂ layer: `W1`, `S4`, `S6`;
- upper MoS₂ layer: `Mo2`, `S3`, `S5`.

Figures made from this dataset should include a visible provenance note, for
example:

```python
data_note=(
    "Bundled QE-format example data · "
    "original calculation provenance not included"
)
```

For research figures, replace these files with traceable calculation outputs
and record the method, pseudopotentials, cutoffs, k-mesh, and calculation ID in
the figure note or accompanying metadata.
