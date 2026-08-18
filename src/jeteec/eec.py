"""Energy-Energy Correlators (EEC).

The two-point EEC of a jet is the pT-weighted distribution of pairwise angular
separations between its constituents:

    EEC(R_L) = sum_{i,j}  (pT_i * pT_j / pT_jet^2)  * delta(R_L - dR_ij)

Conventions made explicit (these are the usual source of disagreement between papers):

* **Normalisation** — weights are divided by ``pT_jet**2``, where ``pT_jet`` is taken as the
  scalar sum of constituent pT (``sum_pt``) rather than the clustered jet pT. The two differ
  slightly; ``sum_pt`` makes ``sum_ij w_ij == 1`` exactly, including the ``i == j`` terms.
* **Pair counting** — all ordered pairs ``(i, j)`` including ``i == j`` are summed, so the
  total weight is exactly 1. Self-pairs sit at ``dR = 0`` and fall outside any log-spaced
  binning, so they are harmlessly dropped by the histogram. Set ``include_self=False`` to
  exclude them explicitly.
* **Angle** — ``dR_ij = sqrt(dEta^2 + dPhi^2)`` with dPhi wrapped to (-pi, pi].

Reference: Energy correlation functions for jet substructure,
Larkoski, Moult, Neill, JHEP 06 (2013) 108, arXiv:1305.0007.
"""

from __future__ import annotations

import awkward as ak
import numpy as np


def _delta_r(constituents: ak.Array) -> tuple[ak.Array, ak.Array]:
    """Return (dR, weight-numerator) for all constituent pairs within each jet.

    weight-numerator is ``pT_i * pT_j``; normalisation is applied by the caller.
    """
    pairs = ak.argcombinations(constituents, 2, replacement=True, axis=-1)
    i, j = ak.unzip(pairs)

    ci = constituents[i]
    cj = constituents[j]

    deta = ci.eta - cj.eta
    dphi = (ci.phi - cj.phi + np.pi) % (2 * np.pi) - np.pi
    dr = np.sqrt(deta**2 + dphi**2)

    # argcombinations gives unordered pairs; off-diagonal terms count twice
    w = ci.pt * cj.pt
    w = ak.where(i == j, w, 2 * w)
    return dr, w


def eec(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the jet-normalised two-point EEC, averaged over all jets.

    Parameters
    ----------
    constituents : ak.Array
        Jagged array of ``Momentum4D`` constituents, one sublist per jet
        (as returned by :func:`jeteec.jets.cluster`, flattened over events).
    bins : np.ndarray, optional
        Bin edges in dR. Defaults to 50 log-spaced bins over [1e-3, 1].
    include_self : bool
        Include ``i == j`` self-pairs. They land at dR = 0 and are dropped by
        log-spaced bins either way; keeping them makes the weights sum to 1.

    Returns
    -------
    centres : np.ndarray
        Geometric bin centres.
    values : np.ndarray
        ``dSigma/dlog(R_L)`` per jet — i.e. summed weights per bin, divided by bin
        width in log space and by the number of jets.
    """
    if bins is None:
        bins = np.logspace(-3, 0, 51)

    dr, w = _delta_r(constituents)

    # Normalise per jet: sum of constituent pT, squared
    sum_pt = ak.sum(constituents.pt, axis=-1)
    w = w / sum_pt**2

    if not include_self:
        mask = dr > 0
        dr, w = dr[mask], w[mask]

    n_jets = len(constituents)
    if n_jets == 0:
        raise ValueError("no jets supplied")

    dr_flat = ak.to_numpy(ak.flatten(dr, axis=None))
    w_flat = ak.to_numpy(ak.flatten(w, axis=None))

    counts, edges = np.histogram(dr_flat, bins=bins, weights=w_flat)

    dlog = np.diff(np.log(edges))
    values = counts / dlog / n_jets
    centres = np.sqrt(edges[:-1] * edges[1:])
    return centres, values


def normalisation_check(constituents: ak.Array) -> float:
    """Return the total EEC weight summed over all pairs, averaged per jet.

    Should be ~1.0 by construction. A useful unit test on any new sample: if it is
    not 1, the constituent collection or the pT normalisation is wrong.
    """
    _, w = _delta_r(constituents)
    sum_pt = ak.sum(constituents.pt, axis=-1)
    total = ak.sum(w, axis=-1) / sum_pt**2
    return float(ak.mean(total))
