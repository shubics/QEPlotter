"""Interactive visualisation for uniform and irreducible k-point meshes."""

from itertools import product

import numpy as np
import plotly.graph_objects as go


def _fold_to_first_bz(fractional, reciprocal):
    """Fold reciprocal coordinates into the Wigner--Seitz first BZ."""
    fractional = np.atleast_2d(np.asarray(fractional, dtype=float))
    reciprocal = np.asarray(reciprocal, dtype=float)
    shifts = np.asarray(list(product(range(-2, 3), repeat=3)), dtype=float)
    output = np.empty((len(fractional), 3), dtype=float)
    chunk_size = 4000
    for start in range(0, len(fractional), chunk_size):
        chunk = fractional[start:start + chunk_size]
        candidates = (chunk[:, None, :] - shifts[None, :, :]) @ reciprocal
        norms = np.einsum("ijk,ijk->ij", candidates, candidates)
        selected = np.argmin(norms, axis=1)
        output[start:start + len(chunk)] = candidates[
            np.arange(len(chunk)), selected]
    return output


def build_kmesh_figure(result, show_full_grid=False, max_full_points=5000):
    """Build a first-BZ plot with weighted IBZ representatives."""
    vertices = np.asarray(result["bz"]["vertices"], dtype=float)
    edge_x, edge_y, edge_z = [], [], []
    for first, second in result["bz"]["edges"]:
        edge_x.extend([vertices[first, 0], vertices[second, 0], None])
        edge_y.extend([vertices[first, 1], vertices[second, 1], None])
        edge_z.extend([vertices[first, 2], vertices[second, 2], None])

    figure = go.Figure()
    figure.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="#68727D", width=4), name="First BZ",
        hoverinfo="skip",
    ))

    total = result["total_grid_points"]
    if show_full_grid:
        indices = np.asarray(result["grid_indices"])
        if len(indices) > max_full_points:
            selection = np.linspace(
                0, len(indices) - 1, max_full_points, dtype=int)
            indices = indices[selection]
        fractional = (
            indices + 0.5 * np.asarray(result["shift"])
        ) / np.asarray(result["mesh"])
        fractional = np.mod(fractional + 0.5, 1.0) - 0.5
        cartesian = _fold_to_first_bz(
            fractional, result["reciprocal_input"])
        figure.add_trace(go.Scatter3d(
            x=cartesian[:, 0], y=cartesian[:, 1], z=cartesian[:, 2],
            mode="markers", name="Full grid",
            marker=dict(size=2.5, color="#87919C", opacity=0.25),
            hoverinfo="skip",
        ))

    fractional = np.asarray([point["frac"] for point in result["points"]])
    cartesian = _fold_to_first_bz(
        fractional, result["reciprocal_input"])
    multiplicities = np.asarray(
        [point["multiplicity"] for point in result["points"]], dtype=float)
    max_multiplicity = max(float(multiplicities.max(initial=1.0)), 1.0)
    sizes = 5.5 + 7.5 * np.sqrt(multiplicities / max_multiplicity)
    customdata = np.column_stack([
        np.arange(1, len(fractional) + 1),
        fractional,
        multiplicities,
        multiplicities / total,
    ])
    figure.add_trace(go.Scatter3d(
        x=cartesian[:, 0], y=cartesian[:, 1], z=cartesian[:, 2],
        mode="markers", name="Irreducible points",
        marker=dict(
            size=sizes, color=multiplicities, colorscale="Cividis",
            colorbar=dict(title="Multiplicity", thickness=12),
            line=dict(color="#E7EAED", width=0.6), opacity=0.95,
        ),
        customdata=customdata,
        hovertemplate=(
            "IBZ #%{customdata[0]:.0f}"
            "<br>k = (%{customdata[1]:.6f}, %{customdata[2]:.6f}, "
            "%{customdata[3]:.6f})"
            "<br>multiplicity = %{customdata[4]:.0f}"
            "<br>weight = %{customdata[5]:.8f}<extra></extra>"
        ),
    ))

    figure.update_layout(
        height=620,
        margin=dict(l=0, r=0, t=35, b=0),
        paper_bgcolor="#15191E",
        plot_bgcolor="#15191E",
        font=dict(color="#E7EAED"),
        legend=dict(orientation="h", y=1.03),
        scene=dict(
            aspectmode="data",
            bgcolor="#15191E",
            xaxis=dict(title="kₓ (Å⁻¹)", showbackground=False),
            yaxis=dict(title="kᵧ (Å⁻¹)", showbackground=False),
            zaxis=dict(title="k_z (Å⁻¹)", showbackground=False),
        ),
    )
    return figure
