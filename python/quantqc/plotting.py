"""Shared plot themes and helpers — equivalent to R misc.R formatting objects."""

from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib as mpl


def dot_plot_style(ax: plt.Axes) -> plt.Axes:
    """Apply the ``dot_plot`` ggplot theme: bw + large fonts."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=16)
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)
    ax.title.set_size(22)
    return ax


def um_plot_style(ax: plt.Axes) -> plt.Axes:
    """Apply the ``um_plot`` ggplot theme: classic + large axis titles."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=12)
    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)
    ax.title.set_size(22)
    return ax


def diverging_cmap(values, midpoint=None):
    """Return a TwoSlopeNorm centred on *midpoint* (default: median)."""
    import numpy as np

    # Drop NaN / Inf before computing statistics
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]

    # Fall back to simple Normalize when there is nothing useful
    if len(clean) == 0:
        return mpl.colors.Normalize(vmin=0, vmax=1)

    vmin, vmax = float(clean.min()), float(clean.max())

    if vmin == vmax:
        return mpl.colors.Normalize(vmin=vmin - 1, vmax=vmax + 1)

    if midpoint is None:
        midpoint = float(np.median(clean))

    # Clamp midpoint strictly between vmin and vmax
    eps = (vmax - vmin) * 1e-6
    midpoint = max(vmin + eps, min(vmax - eps, midpoint))

    return mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=midpoint, vmax=vmax)


MY_COL3 = ["purple", "black"]
