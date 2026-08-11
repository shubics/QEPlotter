"""
Crystal Structure page: upload CIF / POSCAR / QE and view it in 3D (client-side)
together with space group, bond lengths and bond angles.

Heavy 3D rendering runs in the browser via 3Dmol.js; the server only does the
lightweight parsing / symmetry / bond analysis.
"""
import numpy as np
import streamlit as st
import streamlit.components.v1 as components
import re
from functools import partial

from gui.io_helpers import save_file
from gui.theme import remember_tab, remembered_tabs
from qeplotter.structure import (
    read_structure,
    structure_summary,
    get_spacegroup,
    analyse_bonds,
    analyse_angles,
    build_3dmol_html,
)

_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_STRUCTURE_TAB_STATE = "structure_analysis_tab"


def _display_formula(formula):
    """Render atom counts as real chemical subscripts (MoS2 -> MoS₂)."""
    return re.sub(r"\d+", lambda match: match.group(0).translate(_SUBSCRIPT), formula)


def _formula_identity(atoms):
    """Return a material-aware formula, preserving bilayer identity."""
    total = atoms.get_chemical_formula(mode="metal", empirical=False)
    reduced = atoms.get_chemical_formula(mode="metal", empirical=True)
    result = {"display": _display_formula(reduced),
              "total": _display_formula(total), "layers": None}
    try:
        from qeplotter.analysis.bilayer import detect_bilayer
        detected = detect_bilayer(
            np.asarray(atoms.cell[:]), list(atoms.get_chemical_symbols()),
            np.asarray(atoms.get_scaled_positions()))
        if detected["is_bilayer"]:
            lower = atoms[detected["lower"]].get_chemical_formula(
                mode="metal", empirical=True)
            upper = atoms[detected["upper"]].get_chemical_formula(
                mode="metal", empirical=True)
            lower, upper = _display_formula(lower), _display_formula(upper)
            result.update(display=f"{lower} / {upper}", layers=(lower, upper))
    except Exception:
        pass
    return result


def _selected_index(key, frame):
    state = st.session_state.get(key)
    if not state:
        return None
    selection = state.get("selection", {}) if hasattr(state, "get") else state.selection
    rows = selection.get("rows", []) if hasattr(selection, "get") else selection.rows
    if not rows or rows[0] >= len(frame):
        return None
    return int(rows[0])


def _activate_selection(kind):
    st.session_state.structure_active_selection = kind
    remember_tab(
        _STRUCTURE_TAB_STATE, "Bonds" if kind == "bond" else "Angles"
    )


def _clear_selection():
    st.session_state.structure_active_selection = None


def _viewer_highlight(bonds, bond_geometry, angles, angle_geometry):
    kind = st.session_state.get("structure_active_selection")
    frame = bonds if kind == "bond" else angles
    geometry = bond_geometry if kind == "bond" else angle_geometry
    selected_index = _selected_index(f"structure_{kind}_table", frame)
    if selected_index is None or selected_index >= len(geometry):
        return None, None
    row = frame.iloc[selected_index]
    record = geometry[selected_index]
    if kind == "bond":
        labels = [row["atom_1"], row["atom_2"]]
        description = f"Selected bond: {labels[0]} — {labels[1]} ({row['length (Å)']:.3f} Å)"
        value_label = f"{row['length (Å)']:.3f} Å"
    else:
        labels = [row["neighbor_1"], row["vertex"], row["neighbor_2"]]
        description = (f"Selected angle: {labels[0]} — {labels[1]} — {labels[2]} "
                       f"({row['angle (°)']:.2f}°)")
        value_label = f"{row['angle (°)']:.2f}°"
    return {"kind": kind, "indices": record["indices"],
            "points": record["points"],
            "point_elements": record["point_elements"],
            "point_labels": record["point_labels"],
            "value_label": value_label}, description


def _render_summary(atoms):
    s = structure_summary(atoms)
    identity = _formula_identity(atoms)
    c1, c2, c3 = st.columns(3)
    c1.metric("Material", identity["display"])
    c2.metric("Atoms", s["natoms"])
    c3.metric("Volume", f"{s['volume']:.2f} Å³")
    if identity["layers"]:
        st.caption(f"Bilayer recognised automatically • total stoichiometry: {identity['total']}")

    values = [
        ("a", f"{s['a']:.3f} Å"), ("b", f"{s['b']:.3f} Å"),
        ("c", f"{s['c']:.3f} Å"), ("α", f"{s['alpha']:.2f}°"),
        ("β", f"{s['beta']:.2f}°"), ("γ", f"{s['gamma']:.2f}°"),
    ]
    cells = "".join(
        f'<div><span style="color:#A8B0BA;font-size:12px">{name}</span>'
        f'<strong style="display:block;color:#E7EAED">{value}</strong></div>'
        for name, value in values)
    st.markdown(
        '<div style="margin-top:.75rem;padding:.8rem 1rem;border:1px solid #343B44;'
        'border-radius:4px;background:#181C21">'
        '<div style="font-size:12px;color:#A8B0BA;margin-bottom:.55rem">UNIT CELL</div>'
        f'<div style="display:grid;grid-template-columns:repeat(3,1fr);gap:.7rem">{cells}</div>'
        '</div>', unsafe_allow_html=True)
    st.caption(f"Periodic boundaries: {s['pbc']}")


def _render_symmetry(atoms):
    symprec = st.slider("Symmetry tolerance (Å)", 1e-4, 1e-1, 1e-3,
                        format="%.4f", key="sym_tol",
                        help="Larger values are more forgiving of small distortions.",
                        on_change=remember_tab,
                        args=(_STRUCTURE_TAB_STATE, "Symmetry"))
    sg = get_spacegroup(atoms, symprec=symprec)
    if "error" in sg:
        st.warning(f"Could not determine symmetry: {sg['error']}")
        return
    c1, c2 = st.columns(2)
    c1.metric("Space group", f"{sg['international']}")
    c2.metric("Number", f"#{sg['number']}")
    c1, c2 = st.columns(2)
    c1.metric("Crystal system", sg["crystal_system"])
    c2.metric("Point group", sg["pointgroup"])
    st.caption(f"Hall: {sg['hall']}  •  symmetry operations: {sg['n_symmetry_ops']}")


def _render_bonds(df):
    if df.empty:
        st.info("No bonds found at this tolerance. Try increasing it.")
        return
    st.caption(f"{len(df)} bonds detected.")
    # Per element-pair summary
    summary = (df.groupby("elements")["length (Å)"]
                 .agg(["count", "min", "mean", "max"]).round(3).reset_index())
    st.markdown("**Per bond type**")
    st.dataframe(summary, width="stretch", hide_index=True)
    st.markdown("**All bonds**")
    st.caption("Select a row to draw that exact bond in yellow. A suffix such as "
               "`[+1,+0,+0]` identifies the neighbouring periodic cell.")
    st.dataframe(df, width="stretch", hide_index=True, height=280,
                 key="structure_bond_table",
                 on_select=partial(_activate_selection, "bond"),
                 selection_mode="single-row")


def _render_angles(df):
    if df.empty:
        st.info("No angles found at this tolerance.")
        return
    st.caption(f"{len(df)} PBC-aware angles detected. Equivalent/repeated angle "
               "values from different atom combinations are included.")
    st.caption("Select a row: yellow arms show the two bonds, pink marks the vertex, "
               "and the cyan arc shows the measured angle.")
    st.dataframe(df, width="stretch", hide_index=True, height=320,
                 key="structure_angle_table",
                 on_select=partial(_activate_selection, "angle"),
                 selection_mode="single-row")


def _render_stacking(atoms):
    st.caption(
        "Automatic bilayer check with ordered R-type (AA / AB / BA) and "
        "H-type (AA′ / AB′ / A′B) registry analysis."
    )
    try:
        from qeplotter.analysis.bilayer import analyse_stacking
        cell = np.asarray(atoms.cell[:])
        species = list(atoms.get_chemical_symbols())
        frac = np.asarray(atoms.get_scaled_positions())

        result = analyse_stacking(cell, species, frac)
        if not result["is_bilayer"]:
            st.info("This structure does not look like a bilayer. Stacking analysis was skipped.")
            st.caption(result["reason"])
            return

        c1, c2, c3 = st.columns(3)
        c1.metric("Registry", result["label"])
        c2.metric("Orientation", result["family"])
        c3.metric("Interlayer Δz", f"{result.get('spacing', 0.0):.3f} Å")

        st.success(result["bilayer_type"])
        st.caption(
            f"Lower → upper: {result['lower_formula']} / {result['upper_formula']}  •  "
            f"facing interface: {result['interface']}"
        )
        if result["label"] == "General registry":
            shift = result.get("shift")
            if shift is not None:
                st.info(
                    f"{result['description']}  Fractional core shift: "
                    f"({shift[0]:.3f}, {shift[1]:.3f})."
                )
            else:
                st.info(result["description"])
        else:
            st.caption(
                f"{result['description']}  •  confidence: {result['confidence']}"
            )

        with st.expander("Stacking label convention"):
            st.markdown(
                """
                **R-type (parallel):** AA, AB, BA

                **H-type (antiparallel):** AA′, AB′, A′B

                Labels are ordered from the lower layer to the upper layer, so
                AB and BA remain distinct for heterobilayers and Janus systems.
                The facing atomic species are reported separately. Twisted,
                incommensurate, or unsupported lattices are shown as
                **General registry** instead of being forced into AA.
                """
            )
    except Exception as e:
        st.warning(f"Stacking analysis not available for this structure: {e}")


def render_structure():
    st.title("Crystal Structure Explorer")
    st.caption("Upload once, inspect the structure, and get symmetry, bonds, all angles "
               "and bilayer stacking without manual setup.")

    f_struc = st.file_uploader(
        "Structure file",
        type=None,  # POSCAR/CONTCAR have no extension, so accept anything
        key="struct_upload",
        help="Supported: CIF, POSCAR/CONTCAR, Quantum ESPRESSO (.in/.out), XSF, XYZ.",
    )

    if not f_struc:
        st.info("Upload a CIF, POSCAR or Quantum ESPRESSO file to begin.")
        return

    path = save_file(f_struc)
    try:
        atoms = read_structure(path)
    except Exception as e:
        st.error(f"Could not read structure: {e}")
        return

    if atoms.cell.volume <= 1e-6:
        st.warning("This file has no (or zero-volume) unit cell. "
                   "Symmetry and bond analysis may be unreliable.")

    identity = _formula_identity(atoms)
    st.success(f"Loaded {identity['display']} • {len(atoms)} atoms")

    # Widget values and table selections are already in session_state when a
    # rerun begins, allowing the earlier 3D component to react to later tables.
    tol = float(st.session_state.get("bond_tol", 1.15))
    periodic_bonds = bool(st.session_state.get("periodic_bonds", False))
    live_bonds, bond_geometry = analyse_bonds(
        atoms, tol=tol, include_periodic=periodic_bonds)
    live_angles, angle_geometry = analyse_angles(
        atoms, tol=tol, include_periodic=periodic_bonds)
    highlight, highlight_description = _viewer_highlight(
        live_bonds, bond_geometry, live_angles, angle_geometry)
    col_view, col_info = st.columns([1.45, 1], gap="large")

    # ---- 3D viewer (client-side) ----
    with col_view:
        st.subheader("3D preview")
        cc1, cc2, cc3 = st.columns([1, 1, 1.2])
        style = cc1.selectbox("View style", ["ball-stick", "spacefill", "stick", "wireframe"])
        show_cell = cc2.checkbox("Show unit cell", value=True)
        periodic_bonds = cc3.toggle(
            "Periodic bonds", value=False, key="periodic_bonds",
            on_change=_clear_selection,
            help="Off: connect only atoms visible inside the displayed cell. "
                 "On: include bonds to periodic images outside its boundaries.")
        sc1, sc2, sc3 = st.columns(3)
        nx = sc1.number_input("Repeat a", 1, 4, 1)
        ny = sc2.number_input("Repeat b", 1, 4, 1)
        nz = sc3.number_input("Repeat c", 1, 4, 1)

        html = build_3dmol_html(
            atoms, style=style, show_cell=show_cell,
            supercell=(nx, ny, nz), height=480, highlight=highlight,
            bond_tol=tol, periodic_bonds=periodic_bonds,
        )
        components.html(html, height=500, scrolling=False)
        if highlight_description:
            st.info(highlight_description)
        boundary_mode = "including periodic images" if periodic_bonds else "visible atoms only"
        st.caption(f"Drag to rotate • scroll to zoom • click an atom to pin its identity • {boundary_mode}")

    # ---- Analysis panel ----
    with col_info:
        st.subheader("Structure analysis")
        tol = st.slider("Bond tolerance ×(covalent radii)", 0.8, 1.6, 1.15, 0.05,
                        key="bond_tol",
                        on_change=_clear_selection,
                        help="Scales covalent-radius cutoffs used to perceive bonds. "
                             "Lower it if unrelated atoms get bonded. Results update as you drag.")
        m1, m2 = st.columns(2)
        m1.metric("Detected bonds", len(live_bonds))
        m2.metric("All angles", len(live_angles))
        st.caption("The counters and the Bonds/Angles tables update immediately with the tolerance.")
        tabs = remembered_tabs(
            ["Summary", "Symmetry", "Bonds", "Angles", "Stacking"],
            _STRUCTURE_TAB_STATE,
        )
        with tabs[0]:
            _render_summary(atoms)
        with tabs[1]:
            _render_symmetry(atoms)
        with tabs[2]:
            _render_bonds(live_bonds)
        with tabs[3]:
            _render_angles(live_angles)
        with tabs[4]:
            _render_stacking(atoms)
