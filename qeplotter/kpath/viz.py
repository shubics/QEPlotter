"""Client-rendered Plotly Brillouin-zone visualisation."""
import numpy as np
import plotly.graph_objects as go


def _label(label):
    return "Γ" if label == "GAMMA" else label


def build_bz_figure(result):
    vertices = np.asarray(result["bz"]["vertices"])
    edge_x, edge_y, edge_z = [], [], []
    for first, second in result["bz"]["edges"]:
        for axis, target in enumerate((edge_x, edge_y, edge_z)):
            target.extend([vertices[first, axis], vertices[second, axis], None])

    figure = go.Figure()
    figure.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="#68727D", width=4), name="First BZ", hoverinfo="skip"))

    for start, end in result["path"]:
        points = np.vstack([result["point_coords_cart"][start],
                            result["point_coords_cart"][end]])
        figure.add_trace(go.Scatter3d(
            x=points[:, 0], y=points[:, 1], z=points[:, 2], mode="lines",
            line=dict(color="#C2A15A", width=7), showlegend=False,
            hovertemplate=f"{_label(start)} → {_label(end)}<extra></extra>"))

    labels = list(result["point_coords_cart"])
    points = np.vstack([result["point_coords_cart"][label] for label in labels])
    figure.add_trace(go.Scatter3d(
        x=points[:, 0], y=points[:, 1], z=points[:, 2], mode="markers+text",
        marker=dict(size=5, color="#7F9DB9"),
        text=[_label(label) for label in labels], textposition="top center",
        customdata=np.asarray([result["point_coords"][label] for label in labels]),
        hovertemplate="%{text}<br>(%{customdata[0]:.4f}, %{customdata[1]:.4f}, %{customdata[2]:.4f})<extra></extra>",
        name="High-symmetry points"))

    figure.update_layout(
        height=600, margin=dict(l=0, r=0, t=35, b=0),
        paper_bgcolor="#15191E", plot_bgcolor="#15191E",
        font=dict(color="#E7EAED"), legend=dict(orientation="h"),
        scene=dict(aspectmode="data", bgcolor="#15191E",
                   xaxis=dict(title="kₓ (Å⁻¹)", showbackground=False),
                   yaxis=dict(title="kᵧ (Å⁻¹)", showbackground=False),
                   zaxis=dict(title="k_z (Å⁻¹)", showbackground=False)))
    return figure
