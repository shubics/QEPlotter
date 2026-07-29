"""Native high-symmetry K-path recommendation and Brillouin-zone page."""
from io import StringIO

import pandas as pd
import streamlit as st

from gui.io_helpers import save_file
from qeplotter.kpath import (
    build_bz_figure, format_qe_kpoints, parse_path_expression, primary_path,
    recommend_kpath, with_path,
)
from qeplotter.structure import read_structure


def _display_label(label):
    return "Γ" if label == "GAMMA" else label


def _path_text(path):
    groups, current = [], []
    for start, end in path:
        if not current:
            current = [start, end]
        elif current[-1] == start:
            current.append(end)
        else:
            groups.append(current)
            current = [start, end]
    if current:
        groups.append(current)
    return "  |  ".join("–".join(_display_label(label) for label in group)
                         for group in groups)


def _point_table(result):
    used = []
    for start, end in result["path"]:
        for label in (start, end):
            if label not in used:
                used.append(label)
    rows = []
    for label in used:
        x, y, z = result["point_coords"][label]
        rows.append({"Point": _display_label(label), "k₁": x, "k₂": y, "k₃": z})
    return pd.DataFrame(rows)


def _explicit_csv(result):
    rows = []
    for item in result["explicit"]:
        x, y, z = item["frac"]
        rows.append({"segment": item["segment"] + 1,
                     "from": _display_label(item["start"]),
                     "to": _display_label(item["end"]),
                     "label": _display_label(item["label"]) if item["label"] else "",
                     "k1": x, "k2": y, "k3": z})
    buffer = StringIO()
    pd.DataFrame(rows).to_csv(buffer, index=False)
    return buffer.getvalue()


def render_kpath():
    st.title("K-path & Brillouin Zone")
    st.caption("Native QEPlotter engine: spglib primitive-cell extraction, "
               "Setyawan–Curtarolo recipes for all 14 Bravais lattices, and a "
               "separate Wigner–Seitz geometry engine. No SeeK-path dependency.")

    uploaded = st.file_uploader(
        "Structure file", type=None, key="kpath_structure_upload",
        help="CIF, POSCAR/CONTCAR, Quantum ESPRESSO input/output, XSF or XYZ.")
    if not uploaded:
        st.info("Upload a periodic structure to generate its recommended band path.")
        return

    path = save_file(uploaded, subdir="kpath")
    try:
        atoms = read_structure(path)
    except Exception as error:
        st.error(f"Could not read the structure: {error}")
        return
    if atoms.cell.volume <= 1e-8 or not atoms.pbc.any():
        st.error("A non-zero periodic unit cell is required for reciprocal-space analysis.")
        return

    settings, _ = st.columns([1.4, 2])
    with settings:
        symprec = st.number_input("Symmetry tolerance (Å)", min_value=1e-6,
                                  max_value=1e-1, value=1e-3, format="%.5f")
        spacing = st.number_input("Target k-point spacing (Å⁻¹)", min_value=0.005,
                                  max_value=0.5, value=0.05, step=0.005,
                                  help="Smaller spacing creates more explicit k-points.")

    try:
        result = recommend_kpath(atoms, symprec=symprec,
                                 reference_distance=spacing)
    except Exception as error:
        st.error(f"K-path generation failed: {error}")
        return

    st.markdown("#### Path selection")
    variant = st.segmented_control(
        "Choose a path variant",
        ["Full recommended", "Primary branch", "Custom path"],
        default="Full recommended", key="kpath_variant",
        help="Primary branch removes auxiliary cross-links. Custom accepts branches such as Γ-X-M-Γ | R-X.")
    if variant == "Primary branch":
        result = with_path(result, primary_path(result["path"]), spacing)
    elif variant == "Custom path":
        expression = st.text_input(
            "Custom path", value=_path_text(result["path"]),
            key=f"kpath_custom_{uploaded.name}",
            placeholder="Γ-X-M-Γ | R-X",
            help="Use -, – or → between points and | between disconnected branches.")
        try:
            custom_path = parse_path_expression(
                expression, result["point_coords"].keys())
            result = with_path(result, custom_path, spacing)
        except ValueError as error:
            st.error(str(error))
            return

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Bravais case", f"{result['bravais_lattice_extended']} · {result['recipe_variant']}")
    m2.metric("Space group", f"{result['spacegroup_international']} (#{result['spacegroup_number']})")
    m3.metric("Path segments", len(result["path"]))
    m4.metric("Explicit points", len(result["explicit"]))
    st.success(f"Conventional path: {result['convention']}")

    st.markdown("#### Conventional recommended path")
    st.code(_path_text(result["path"]), language=None)

    col_bz, col_data = st.columns([1.55, 1], gap="large")
    with col_bz:
        st.plotly_chart(build_bz_figure(result), use_container_width=True,
                        config={"displaylogo": False, "scrollZoom": True})
        st.caption("Drag to rotate • scroll to zoom • orange: recommended path • blue: special points")
    with col_data:
        tabs = st.tabs(["Special points", "QE input", "Details"])
        with tabs[0]:
            st.dataframe(_point_table(result), hide_index=True, width="stretch",
                         column_config={"k₁": st.column_config.NumberColumn(format="%.6f"),
                                        "k₂": st.column_config.NumberColumn(format="%.6f"),
                                        "k₃": st.column_config.NumberColumn(format="%.6f")})
            st.caption("Fractional coordinates are expressed in the uploaded cell's reciprocal basis.")
        with tabs[1]:
            output_mode = st.radio(
                "QE card format", ["Band path (crystal_b)", "Explicit points (crystal)"],
                horizontal=True, help="crystal_b is recommended for pw.x bands calculations.")
            qe_text = format_qe_kpoints(
                result, explicit=output_mode.startswith("Explicit"))
            st.code(qe_text, language=None, line_numbers=True)
            st.download_button("Download K_POINTS", qe_text, "K_POINTS.dat",
                               "text/plain", width="stretch")
            st.download_button("Download explicit CSV", _explicit_csv(result),
                               "kpath.csv", "text/csv", width="stretch")
        with tabs[2]:
            parameters = result["recipe_parameters"]
            parameter_text = (", ".join(f"{name}={value:.8f}" for name, value in parameters.items())
                              if parameters else "not required for this case")
            st.markdown(
                f"**Crystal system:** {result['crystal_system'].title()}  \n"
                f"**Geometry engine:** reciprocal Wigner–Seitz / Voronoi  \n"
                f"**Band-path engine:** {result['convention']}  \n"
                f"**Extended Bravais case:** {result['bravais_lattice_extended']} ({result['recipe_variant']})  \n"
                f"**Metric-dependent parameters:** {parameter_text}  \n"
                f"**Inversion symmetry:** {'Yes' if result['has_inversion'] else 'No'}  \n"
                f"**BZ vertices:** {len(result['bz']['vertices'])}  \n"
                f"**BZ edges:** {len(result['bz']['edges'])}  \n"
                f"**Coordinate basis:** {result['coordinate_basis']}  \n"
                f"**Reciprocal convention:** `{result['reciprocal_convention']}`  \n"
                f"**Basis mapping residual:** {result['basis_mapping_error']:.3e} Å²")
            st.caption("The first BZ is calculated as the Wigner–Seitz cell of the "
                       "primitive reciprocal lattice. BZ vertices are geometry only; "
                       "they are never used as high-symmetry path points.")
