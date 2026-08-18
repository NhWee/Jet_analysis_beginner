"""Plotting helpers — publication-style figures via mplhep."""

from __future__ import annotations

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np


def use_style(experiment: str = "ATLAS") -> None:
    """Apply an experiment plot style: ATLAS, CMS, ALICE or LHCb."""
    styles = {
        "ATLAS": hep.style.ATLAS,
        "CMS": hep.style.CMS,
        "ALICE": hep.style.ALICE,
        "LHCb": hep.style.LHCb2,
    }
    if experiment not in styles:
        raise ValueError(f"unknown experiment {experiment!r}; choose from {sorted(styles)}")
    hep.style.use(styles[experiment])


def plot_eec(
    centres: np.ndarray,
    values: np.ndarray,
    label: str | None = None,
    ax: "plt.Axes | None" = None,
    **kwargs,
) -> "plt.Axes":
    """Plot an EEC distribution on log-log axes.

    The characteristic shape has three regions: a free-hadron/small-angle region,
    a scaling (perturbative) region, and a turnover near the jet radius.
    """
    ax = ax or plt.gca()
    ax.plot(centres, values, label=label, **kwargs)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R_L$")
    ax.set_ylabel(r"$\frac{1}{\sigma}\frac{d\sigma}{d\log R_L}$")
    if label:
        ax.legend()
    return ax


def save(fig: "plt.Figure", path: str, dpi: int = 300) -> None:
    """Save a figure with tight bounding box (PDF for papers, PNG for notes)."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
