"""Small kinematic helpers shared by the clustering, correlator and analysis code.

Kept free of fastjet so that :mod:`jeteec.eec` can be used on its own.
"""

from __future__ import annotations

import awkward as ak
import numpy as np
import vector

vector.register_awkward()


def wrap_phi(phi):
    """Map an azimuthal angle into [-pi, pi)."""
    return (phi + np.pi) % (2 * np.pi) - np.pi


def delta_phi(phi1, phi2):
    """``phi1 - phi2`` wrapped into [-pi, pi)."""
    return wrap_phi(phi1 - phi2)


def delta_r(a, b):
    """``sqrt(d_eta^2 + d_phi^2)`` between records ``a`` and ``b`` (pseudorapidity)."""
    return np.sqrt((a.eta - b.eta) ** 2 + delta_phi(a.phi, b.phi) ** 2)


def sum_p4(v: ak.Array, axis: int = -1) -> ak.Array:
    """Sum four-vectors along ``axis`` (default: the innermost list)."""
    return ak.zip(
        {"px": ak.sum(v.px, axis=axis), "py": ak.sum(v.py, axis=axis),
         "pz": ak.sum(v.pz, axis=axis), "E": ak.sum(v.E, axis=axis)},
        with_name="Momentum4D",
    )
