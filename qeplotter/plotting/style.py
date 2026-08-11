"""Shared validation and styling helpers for Matplotlib plots."""

from __future__ import annotations

import math

from matplotlib import colormaps
from matplotlib.colors import ListedColormap


LEGEND_LOCATIONS = (
    "best",
    "upper right",
    "upper left",
    "lower left",
    "lower right",
    "center right",
    "center left",
    "lower center",
    "upper center",
    "center",
)


def colormap_colors(cmap_name, count):
    """Return visually separated colours sampled from a Matplotlib colormap."""
    if count < 0:
        raise ValueError("count must be zero or greater")
    if count == 0:
        return []

    cmap = colormaps.get_cmap(cmap_name)
    if isinstance(cmap, ListedColormap) and cmap.N <= 20:
        denominator = max(cmap.N - 1, 1)
        return [cmap((index % cmap.N) / denominator) for index in range(count)]
    if count == 1:
        return [cmap(0.65)]
    return [cmap(0.08 + 0.84 * index / (count - 1)) for index in range(count)]


def display_text(value, default):
    """Use a non-empty custom label, otherwise retain the plot default."""
    if value is None:
        return default
    text = str(value).strip()
    return text if text else default


def figure_size(figsize, default):
    """Validate a user-provided Matplotlib ``figsize`` pair."""
    if figsize is None:
        return default
    try:
        width, height = (float(value) for value in figsize)
    except (TypeError, ValueError):
        raise ValueError("figsize must contain exactly two numeric values") from None
    if not (math.isfinite(width) and math.isfinite(height)):
        raise ValueError("figsize values must be finite")
    if width <= 0 or height <= 0:
        raise ValueError("figsize values must be greater than zero")
    return width, height


def apply_axis_style(
    ax,
    *,
    default_title,
    default_xlabel,
    default_ylabel,
    plot_title=None,
    x_label=None,
    y_label=None,
    show_title=True,
    show_grid=True,
    grid_alpha=0.3,
):
    """Apply the plot-wide text and grid options to an axis."""
    ax.set_title(display_text(plot_title, default_title) if show_title else "")
    ax.set_xlabel(display_text(x_label, default_xlabel))
    ax.set_ylabel(display_text(y_label, default_ylabel))
    if show_grid:
        ax.grid(True, ls="--", alpha=grid_alpha)
    else:
        ax.grid(False)


def apply_legend(
    ax,
    *,
    show_legend=True,
    location="best",
    title=None,
    **kwargs,
):
    """Create or remove a categorical legend with validated placement."""
    existing = ax.get_legend()
    if not show_legend:
        if existing is not None:
            existing.remove()
        return None

    if location not in LEGEND_LOCATIONS:
        allowed = ", ".join(LEGEND_LOCATIONS)
        raise ValueError(f"Unknown legend location '{location}'. Use one of: {allowed}")

    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    return ax.legend(
        handles,
        labels,
        loc=location,
        title=str(title).strip() if title and str(title).strip() else None,
        **kwargs,
    )
