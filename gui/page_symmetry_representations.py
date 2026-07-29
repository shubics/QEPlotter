"""Interactive Γ-point symmetry and orbital-representation analyser."""
from __future__ import annotations

import hashlib

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from gui.io_helpers import save_file
from qeplotter.structure import build_3dmol_html, read_structure
from qeplotter.structure.bonds import atom_labels
from qeplotter.symmetry import GammaRepresentationAnalyzer, SymmetryAnalysisError


@st.cache_resource(show_spinner=False)
def _cached_analyzer(path: str, content_digest: str, symprec: float):
    del content_digest  # Included to invalidate the cache when file bytes change.
    return GammaRepresentationAnalyzer(read_structure(path), symprec=symprec)


def _display_irrep(row):
    if row["label_source"] == "conventional":
        return row["label"]
    return f"{row['gamma']} ({row['label']}-like)"


def _decomposition_text(result):
    terms = []
    for row in result["decomposition"]:
        prefix = f"{row['multiplicity']}×" if row["multiplicity"] > 1 else ""
        terms.append(f"{prefix}{_display_irrep(row)}")
    return " ⊕ ".join(terms)


def _selected_orbit_highlight(analyzer, orbit):
    labels = atom_labels(analyzer.atoms)
    return {
        "kind": "orbit",
        "points": [
            analyzer.atoms.positions[index].tolist() for index in orbit.indices
        ],
        "point_labels": [labels[index] for index in orbit.indices],
    }


def _render_decomposition(result):
    st.markdown(
        f"""
        <div class="sym-result">
          <span class="sym-kicker">RESULT AT Γ</span>
          <div class="sym-equation">Γ<sub>{result['basis']}</sub> =
          {_decomposition_text(result)}</div>
          <div class="sym-caption">{result['dimension']} basis functions are
          reproduced exactly by the irreducible components below.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    frame = pd.DataFrame([
        {
            "irrep": _display_irrep(row),
            "dimension": row["dimension"],
            "copies": row["multiplicity"],
            "basis states": row["states"],
            "meaning": row["basis_hint"],
        }
        for row in result["decomposition"]
    ])
    st.dataframe(frame, hide_index=True, width="stretch")


def _render_characters(analyzer, result):
    st.subheader("How the reducible representation was built")
    st.caption(
        "For every symmetry class, QEPlotter counts atoms that remain on the "
        "same site and traces the selected orbital transformation. Moved basis "
        "functions contribute zero to the trace."
    )
    st.dataframe(
        pd.DataFrame(result["character_rows"]),
        hide_index=True,
        width="stretch",
    )
    st.info(
        f"The complete {analyzer.pointgroup} group is used: "
        f"{analyzer.operation_count} operations in {analyzer.class_count} classes. "
        "No symmetry class is silently omitted."
    )
    with st.expander("Computed character table"):
        st.dataframe(
            pd.DataFrame(analyzer.character_table_rows()),
            hide_index=True,
            width="stretch",
        )
        if analyzer.has_internal_labels:
            st.caption(
                "Γ labels are deterministic QEPlotter labels. A/E/T indicates "
                "the real irrep dimension; conventional Mulliken labels are only "
                "shown when the mapping is unambiguous."
            )


def _render_salcs(analyzer, result):
    st.subheader("Symmetry-adapted linear combinations")
    st.caption(
        "Each row is an orthonormal combination of the selected atomic basis. "
        "Different bases inside a degenerate irrep are equivalent rotations of "
        "the same symmetry subspace."
    )
    if not result["salcs"]:
        if result["dimension"] <= 48:
            with st.spinner("Generating projection-operator combinations…"):
                analyzer.generate_salcs(result)
        elif not st.button(
            "Generate SALCs",
            type="primary",
            help="Large orbital bases are generated on demand to keep the "
                 "rest of the page responsive.",
        ):
            st.info(
                f"This {result['dimension']}-dimensional basis is large. "
                "Generate its SALCs when you need the explicit coefficients; "
                "the character and irrep results above are already complete."
            )
            return
        else:
            with st.spinner("Generating projection-operator combinations…"):
                analyzer.generate_salcs(result)
    choices = {
        _display_irrep(row): row["gamma"] for row in result["decomposition"]
    }
    selected = st.selectbox("Irrep to inspect", list(choices), key="salc_irrep")
    expressions = result["salcs"][choices[selected]]
    visible = expressions[:24]
    st.dataframe(
        pd.DataFrame([
            {"combination": row["salc"], "normalised coefficients": row["expression"]}
            for row in visible
        ]),
        hide_index=True,
        width="stretch",
    )
    if len(expressions) > len(visible):
        st.caption(
            f"Showing the first {len(visible)} of {len(expressions)} combinations."
        )


def _render_matching(analyzer, result):
    st.subheader("Which orbital sets can interact?")
    st.caption(
        "Two sets are symmetry-compatible when they share at least one irrep. "
        "This is a selection rule; energy and spatial overlap still determine "
        "whether the physical interaction is strong."
    )
    with st.spinner("Comparing s, p and d channels…"):
        rows = analyzer.compatibility(result)
    frame = pd.DataFrame(rows)
    allowed_only = st.toggle(
        "Show symmetry-allowed matches only", value=True, key="allowed_only"
    )
    if allowed_only:
        frame = frame[frame["symmetry allowed"] == "Yes"]
    st.dataframe(frame, hide_index=True, width="stretch", height=360)
    source_orbit = result["orbit"].label
    st.caption(
        f"Source: {source_orbit} · {result['basis_title']}. "
        "The source row itself is retained as a useful consistency check."
    )


def _render_method(analyzer):
    st.subheader("What this tool calculates")
    st.markdown(
        """
1. Standardises the uploaded structure to a primitive cell.
2. Finds the space group, point group, Wyckoff orbits and site symmetries.
3. Builds the full transformation matrices for the selected `s`, `p`, `d`
   or displacement basis.
4. Computes the point-group irreducible characters numerically from the finite
   group and decomposes the reducible representation.
5. Applies projection operators to generate normalised SALCs.
        """
    )
    st.warning(
        "This is a structure-derived Γ-point analysis. It does not assign "
        "irreps to individual QE electronic bands without their wavefunctions, "
        "and it does not include spinor/double-group representations yet."
    )
    st.caption(
        f"Primitive-cell convention · {analyzer.operation_count} point operations · "
        f"symmetry tolerance {analyzer.symprec:g} Å"
    )


def render_symmetry_representations():
    st.title("Symmetry & Orbital Representations")
    st.caption(
        "Γ-point character decomposition, symmetry-adapted combinations, and "
        "orbital compatibility from a periodic crystal structure."
    )

    uploaded = st.file_uploader(
        "Structure file",
        type=None,
        key="symrep_upload",
        help="CIF, POSCAR/CONTCAR, Quantum ESPRESSO input/output, XSF or XYZ.",
    )
    if not uploaded:
        st.info(
            "Upload a periodic structure. QEPlotter will find its primitive "
            "symmetry orbits automatically."
        )
        return

    path = save_file(uploaded, subdir="symmetry-representations")
    content_digest = hashlib.sha1(uploaded.getvalue()).hexdigest()
    symprec = st.slider(
        "Symmetry tolerance (Å)",
        min_value=1e-5,
        max_value=1e-1,
        value=1e-3,
        format="%.5f",
        key="symrep_tol",
        help="Increase slightly for structures with small numerical distortions.",
    )
    try:
        with st.spinner("Building the complete symmetry representation…"):
            analyzer = _cached_analyzer(
                path, content_digest, float(symprec)
            )
    except Exception as error:
        st.error(f"Symmetry representation could not be built: {error}")
        return

    orbit_labels = {orbit.label: orbit.id for orbit in analyzer.orbits}
    controls, context = st.columns([1.1, 1], gap="large")
    with controls:
        selected_label = st.selectbox(
            "Atomic symmetry orbit",
            list(orbit_labels),
            help="A complete Wyckoff orbit is selected so the basis remains "
                 "closed under every point-group operation.",
        )
        basis = st.selectbox(
            "Basis to transform",
            list(analyzer.supported_bases),
            format_func=lambda key: analyzer.basis_titles[key],
        )
        orbit = analyzer.orbit(orbit_labels[selected_label])
        try:
            result = analyzer.analyse(orbit.id, basis)
        except SymmetryAnalysisError as error:
            st.error(str(error))
            return
    with context:
        html = build_3dmol_html(
            analyzer.atoms,
            style="ball-stick",
            show_cell=True,
            height=300,
            highlight=_selected_orbit_highlight(analyzer, orbit),
            periodic_bonds=False,
        )
        components.html(html, height=320, scrolling=False)
        st.caption(
            "The highlighted atoms form the selected complete symmetry orbit "
            "in the standardised primitive cell."
        )

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Point group", analyzer.pointgroup)
    m2.metric("Space group", f"{analyzer.spacegroup} (#{analyzer.spacegroup_number})")
    m3.metric("Primitive atoms", len(analyzer.atoms))
    m4.metric("Representation size", result["dimension"])

    tabs = st.tabs([
        "Result",
        "Character calculation",
        "SALCs",
        "Orbital matching",
        "Method & limits",
    ])
    with tabs[0]:
        _render_decomposition(result)
        st.caption(
            "Degenerate dimensions are kept together; repeated copies are "
            "shown explicitly instead of being removed as duplicates."
        )
    with tabs[1]:
        _render_characters(analyzer, result)
    with tabs[2]:
        _render_salcs(analyzer, result)
    with tabs[3]:
        _render_matching(analyzer, result)
    with tabs[4]:
        _render_method(analyzer)
