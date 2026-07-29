"""Irreducible uniform k-grid workspace."""

from io import StringIO
import re

import numpy as np
import pandas as pd
import streamlit as st

from gui.io_helpers import save_file
from qeplotter.kmesh import (
    build_kmesh_figure,
    format_qe_automatic,
    format_qe_ir_kpoints,
    full_grid_points,
    irreducible_kmesh,
    orbit_members,
)
from qeplotter.structure import read_structure, structure_summary


def _to_csv(frame):
    buffer = StringIO()
    frame.to_csv(buffer, index=False)
    return buffer.getvalue()


def _ibz_frame(result, limit=None):
    rows = []
    points = result["points"] if limit is None else result["points"][:limit]
    for point in points:
        k1, k2, k3 = point["frac"]
        rows.append({
            "IBZ #": point["index"],
            "k₁": k1,
            "k₂": k2,
            "k₃": k3,
            "multiplicity": point["multiplicity"],
            "weight": point["normalized_weight"],
        })
    return pd.DataFrame(rows)


def _members_frame(result, ir_index):
    rows = []
    for member in orbit_members(result, ir_index):
        k1, k2, k3 = member["frac"]
        a1, a2, a3 = member["address"]
        rows.append({
            "full-grid #": member["grid_index"] + 1,
            "address₁": a1,
            "address₂": a2,
            "address₃": a3,
            "k₁": k1,
            "k₂": k2,
            "k₃": k3,
        })
    return pd.DataFrame(rows)


def _full_mapping_frame(result):
    rows = []
    for point in full_grid_points(result):
        k1, k2, k3 = point["frac"]
        a1, a2, a3 = point["address"]
        rows.append({
            "full-grid #": point["grid_index"] + 1,
            "address₁": a1,
            "address₂": a2,
            "address₃": a3,
            "k₁": k1,
            "k₂": k2,
            "k₃": k3,
            "IBZ #": point["irreducible_index"],
        })
    return pd.DataFrame(rows)


def _magnetic_input_hint(uploaded):
    """Find QE flags that require an explicit time-reversal decision."""
    try:
        text = uploaded.getvalue().decode("utf-8", errors="ignore").lower()
    except Exception:
        return []
    hints = []
    if re.search(r"\bnspin\s*=\s*2\b", text):
        hints.append("nspin=2")
    if re.search(r"\bnoncolin\s*=\s*\.?true\.?", text):
        hints.append("noncolin=.true.")
    magnetizations = re.findall(
        r"\bstarting_magnetization\s*\(\s*\d+\s*\)\s*=\s*"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[de][+-]?\d+)?)",
        text,
    )
    if any(abs(float(value.replace("d", "e"))) > 1e-14
           for value in magnetizations):
        hints.append("non-zero starting_magnetization")
    return hints


def _suggest_mesh(atoms, spacing, slab_mode):
    reciprocal_lengths = np.linalg.norm(
        2.0 * np.pi * np.linalg.inv(np.asarray(atoms.cell[:])).T, axis=1)
    mesh = np.maximum(1, np.ceil(reciprocal_lengths / spacing).astype(int))
    if slab_mode:
        mesh[2] = 1
    return tuple(int(value) for value in mesh)


def render_kmesh():
    st.title("K-grid & Irreducible Brillouin Zone")
    st.caption(
        "Generate a complete uniform k-grid, reduce it with the uploaded "
        "crystal's actual symmetry, inspect every orbit, and export weighted "
        "Quantum ESPRESSO input. This engine is separate from band paths."
    )

    uploaded = st.file_uploader(
        "Structure file",
        type=None,
        key="kmesh_structure_upload",
        help="CIF, POSCAR/CONTCAR, Quantum ESPRESSO input/output, XSF or XYZ.",
    )
    if not uploaded:
        st.info(
            "Upload a periodic structure, then choose the k-grid and shift. "
            "The result updates when a setting changes."
        )
        return

    path = save_file(uploaded, subdir="kmesh")
    try:
        atoms = read_structure(path)
    except Exception as error:
        st.error(f"Could not read the structure: {error}")
        return
    if atoms.cell.volume <= 1e-8:
        st.error("A finite, non-zero unit cell is required.")
        return

    summary = structure_summary(atoms)
    st.markdown("#### Grid definition")
    definition_mode = st.segmented_control(
        "Grid input",
        ["Manual mesh", "Target spacing"],
        default="Manual mesh",
        key="kmesh_definition_mode",
        help="Target spacing estimates each nk from the reciprocal-vector length.",
    )

    lengths = atoms.cell.lengths()
    likely_slab = bool(
        lengths[2] > 2.2 * max(min(lengths[0], lengths[1]), 1e-8))
    slab_mode = False
    if definition_mode == "Target spacing":
        left, right = st.columns([1, 1.7])
        with left:
            target_spacing = st.number_input(
                "Target spacing (Å⁻¹)",
                min_value=0.005,
                max_value=1.0,
                value=0.15,
                step=0.005,
                format="%.3f",
                key="kmesh_target_spacing",
            )
        with right:
            slab_mode = st.checkbox(
                "Treat as a 2D slab (nk₃ = 1)",
                value=likely_slab,
                key="kmesh_slab_mode",
                help="Use only when the third cell direction is vacuum/non-periodic.",
            )
        suggested = _suggest_mesh(atoms, target_spacing, slab_mode)
        st.caption(
            f"Suggested mesh from reciprocal-vector lengths: "
            f"`{suggested[0]} × {suggested[1]} × {suggested[2]}`"
        )
        mesh = suggested
    else:
        default_n3 = 1 if likely_slab else 6
        c1, c2, c3 = st.columns(3)
        nk1 = c1.number_input(
            "nk₁", min_value=1, max_value=200, value=6, step=1,
            key="kmesh_nk1")
        nk2 = c2.number_input(
            "nk₂", min_value=1, max_value=200, value=6, step=1,
            key="kmesh_nk2")
        nk3 = c3.number_input(
            "nk₃", min_value=1, max_value=200, value=default_n3, step=1,
            key="kmesh_nk3")
        mesh = (int(nk1), int(nk2), int(nk3))

    shift_col, symmetry_col = st.columns([1.35, 1])
    with shift_col:
        st.markdown("##### Half-grid shift")
        s1, s2, s3 = st.columns(3)
        shift = (
            int(s1.checkbox("shift k₁", key="kmesh_shift1")),
            int(s2.checkbox("shift k₂", key="kmesh_shift2")),
            int(s3.checkbox("shift k₃", key="kmesh_shift3")),
        )
        st.caption(
            "`0` includes the unshifted mesh origin; `1` moves that axis by "
            "half of one grid step."
        )
    with symmetry_col:
        symprec = st.number_input(
            "Symmetry tolerance (Å)",
            min_value=1e-7,
            max_value=1e-1,
            value=1e-3,
            format="%.6f",
            key="kmesh_symprec",
        )
        time_reversal = st.checkbox(
            "Include time-reversal symmetry",
            value=True,
            key="kmesh_time_reversal",
            help="Makes k and −k equivalent. Disable for calculations whose "
                 "magnetic order breaks time reversal.",
        )

    magnetic_hints = _magnetic_input_hint(uploaded)
    if magnetic_hints and time_reversal:
        st.warning(
            "This QE input contains magnetic/non-collinear indicators "
            f"({', '.join(magnetic_hints)}). Confirm that the physical state "
            "preserves time reversal, or disable it above."
        )
    elif time_reversal:
        st.caption(
            "Time reversal is appropriate for non-magnetic calculations. "
            "SOC by itself does not necessarily break it; magnetic order can."
        )

    try:
        result = irreducible_kmesh(
            atoms,
            mesh=mesh,
            shift=shift,
            symprec=symprec,
            time_reversal=time_reversal,
        )
    except Exception as error:
        st.error(f"Irreducible k-grid generation failed: {error}")
        return

    reduction = 100.0 * (
        1.0 - result["irreducible_count"] / result["total_grid_points"])
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Full grid", f"{result['total_grid_points']:,}")
    m2.metric("Irreducible points", f"{result['irreducible_count']:,}")
    m3.metric("Reduction", f"{reduction:.1f}%")
    m4.metric(
        f"Space group #{result['spacegroup_number']}",
        result["spacegroup_international"],
    )

    if result["full_crystal_symmetry_preserved"]:
        st.success(
            "The selected grid is compatible with every detected point "
            "rotation of the crystal."
        )
    else:
        st.warning(
            f"This mesh/shift is not invariant under the full crystal point "
            f"group. {result['dropped_spatial_rotations']} of "
            f"{result['detected_unique_point_rotations']} rotations were "
            "excluded because they do not map the selected grid onto itself. "
            "The reported reduction is exact for the grid-preserving subgroup."
        )

    show_full_grid = st.checkbox(
        "Show the complete grid behind the IBZ representatives",
        value=result["total_grid_points"] <= 2000,
        key="kmesh_show_full",
        help="Large meshes are visually sampled to keep the 3D view responsive.",
    )
    plot_col, table_col = st.columns([1.45, 1], gap="large")
    with plot_col:
        st.plotly_chart(
            build_kmesh_figure(result, show_full_grid=show_full_grid),
            use_container_width=True,
            config={"displaylogo": False, "scrollZoom": True},
        )
        st.caption(
            "Drag to rotate • scroll to zoom • marker size and colour show "
            "multiplicity • coordinates are folded into the first BZ"
        )
    with table_col:
        st.markdown("#### Irreducible points")
        table_limit = 10_000
        ibz_frame = _ibz_frame(result, limit=table_limit)
        st.dataframe(
            ibz_frame,
            hide_index=True,
            width="stretch",
            height=525,
            column_config={
                "k₁": st.column_config.NumberColumn(format="%.8f"),
                "k₂": st.column_config.NumberColumn(format="%.8f"),
                "k₃": st.column_config.NumberColumn(format="%.8f"),
                "weight": st.column_config.NumberColumn(format="%.10f"),
            },
        )
        if result["irreducible_count"] > table_limit:
            st.caption(
                f"The preview shows the first {table_limit:,} of "
                f"{result['irreducible_count']:,} irreducible points."
            )
        if result["irreducible_count"] <= 200_000:
            download_frame = (
                ibz_frame if result["irreducible_count"] <= table_limit
                else _ibz_frame(result)
            )
            st.download_button(
                "Download IBZ table",
                _to_csv(download_frame),
                "irreducible_kpoints.csv",
                "text/csv",
                width="stretch",
            )
        else:
            st.info(
                "Browser CSV export is limited to 200,000 irreducible points "
                "to keep the interactive session responsive."
            )

    st.markdown("#### Export and verification")
    tab_qe, tab_orbit, tab_mapping, tab_method = st.tabs(
        ["QE input", "Orbit inspector", "Full mapping", "Method & checks"])

    with tab_qe:
        output_mode = st.radio(
            "Output",
            ["Recommended automatic grid", "Explicit normalized IBZ",
             "Explicit multiplicity IBZ"],
            horizontal=True,
            key="kmesh_qe_output",
        )
        if output_mode == "Recommended automatic grid":
            qe_text = format_qe_automatic(result)
            file_name = "K_POINTS_automatic.dat"
            st.info(
                "For ordinary QE SCF/NSCF calculations this is the safest "
                "input: QE performs its own symmetry handling. The explicit "
                "cards are useful for inspection and reproducible custom grids."
            )
        else:
            weight_mode = (
                "normalized" if "normalized" in output_mode.lower()
                else "multiplicity")
            qe_text = format_qe_ir_kpoints(result, weight_mode=weight_mode)
            file_name = "K_POINTS_irreducible.dat"
        if not result["time_reversal"]:
            st.warning(
                "To reproduce this no-time-reversal grid with pw.x, set "
                "`noinv=.true.` in `&SYSTEM`. The downloaded card includes "
                "the same reminder as QE comment lines."
            )
        st.code(qe_text, language=None, line_numbers=True)
        st.download_button(
            "Download QE K_POINTS",
            qe_text,
            file_name,
            "text/plain",
            width="stretch",
        )

    with tab_orbit:
        selected = st.selectbox(
            "Irreducible representative",
            options=range(result["irreducible_count"]),
            format_func=lambda index: (
                f"IBZ #{index + 1} · multiplicity "
                f"{result['points'][index]['multiplicity']}"
            ),
            key="kmesh_orbit_selection",
        )
        member_frame = _members_frame(result, selected)
        st.dataframe(
            member_frame,
            hide_index=True,
            width="stretch",
            column_config={
                "k₁": st.column_config.NumberColumn(format="%.8f"),
                "k₂": st.column_config.NumberColumn(format="%.8f"),
                "k₃": st.column_config.NumberColumn(format="%.8f"),
            },
        )
        st.caption(
            "These are all points in the complete grid represented by the "
            "selected IBZ point."
        )

    with tab_mapping:
        if result["total_grid_points"] <= 200_000:
            prepare_mapping = st.checkbox(
                "Prepare the complete mapping table",
                value=result["total_grid_points"] <= 10_000,
                key="kmesh_prepare_mapping",
                help="Large CSV tables are generated only on request.",
            )
            if prepare_mapping:
                mapping_frame = _full_mapping_frame(result)
                st.dataframe(
                    mapping_frame.head(10_000),
                    hide_index=True,
                    width="stretch",
                    column_config={
                        "k₁": st.column_config.NumberColumn(format="%.8f"),
                        "k₂": st.column_config.NumberColumn(format="%.8f"),
                        "k₃": st.column_config.NumberColumn(format="%.8f"),
                    },
                )
                if len(mapping_frame) > 10_000:
                    st.caption(
                        "The preview shows the first 10,000 rows; the download "
                        "contains the complete mapping."
                    )
                st.download_button(
                    "Download full grid → IBZ mapping",
                    _to_csv(mapping_frame),
                    "full_grid_to_ibz.csv",
                    "text/csv",
                    width="stretch",
                )
            else:
                st.caption(
                    "Enable the option above when you need every full-grid "
                    "address and its representative."
                )
        else:
            st.info(
                "The complete mapping is omitted from the browser for this "
                "very large grid. Reduce the mesh below 200,000 points to "
                "inspect or download every row."
            )

    with tab_method:
        weight_sum = sum(
            point["normalized_weight"] for point in result["points"])
        multiplicity_sum = sum(
            point["multiplicity"] for point in result["points"])
        st.markdown(
            f"**Structure:** {summary['reduced_formula']} · "
            f"{summary['natoms']} atoms  \n"
            f"**Crystal system:** {result['crystal_system'].title()}  \n"
            f"**Engine:** {result['engine']}  \n"
            f"**Symmetry source:** {result['symmetry_source']}  \n"
            f"**Detected space-group operations:** "
            f"{result['detected_spacegroup_operations']}  \n"
            f"**Unique spatial point rotations:** "
            f"{result['detected_unique_point_rotations']}  \n"
            f"**Grid-compatible spatial rotations:** "
            f"{result['compatible_spatial_rotations']}  \n"
            f"**Final equivalence transforms:** "
            f"{result['equivalence_transforms']}  \n"
            f"**Time reversal:** "
            f"{'Included' if result['time_reversal'] else 'Not included'}  \n"
            f"**Magnetic symmetry:** structural spatial symmetry only; "
            f"magnetic moments are not inferred from a geometry-only file  \n"
            f"**Coordinate basis:** {result['coordinate_basis']}  \n"
            f"**Reciprocal convention:** "
            f"`{result['reciprocal_convention']}`"
        )
        c1, c2, c3 = st.columns(3)
        c1.metric("Σ multiplicity", f"{multiplicity_sum:,}")
        c2.metric("Full-grid size", f"{result['total_grid_points']:,}")
        c3.metric("Σ normalized weight", f"{weight_sum:.12f}")
        st.success(
            "Internal checks passed: every full-grid point maps to exactly "
            "one representative, multiplicities cover the full grid, and "
            "normalized weights sum to one."
        )
