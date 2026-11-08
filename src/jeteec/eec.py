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
* **Angle** — ``dR_ij = sqrt(dEta^2 + dPhi^2)`` with dPhi wrapped to [-pi, pi).
* **Uncertainty** — pairs inside one jet are correlated, jets are (to a good approximation)
  not. The statistical error of a bin is therefore the spread of the *per-jet* histograms
  divided by sqrt(n_jets), not sqrt of the summed pair weights. See :func:`eec_per_jet`.

References:
  Komiske, Moult, Thaler, Zhu, "Analyzing N-point energy correlators inside
  jets with CMS Open Data", Phys. Rev. Lett. 130 (2023) 051901,
  arXiv:2201.07800 -- the two-point EEC in jets, as computed here.
  Larkoski, Salam, Thaler, "Energy correlation functions for jet
  substructure", JHEP 06 (2013) 108, arXiv:1305.0007 -- the related e2/C2
  discriminants, not the R_L distribution.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

from .kinematics import delta_r


def _bins(bins: np.ndarray | None) -> np.ndarray:
    """Default binning: 50 log-spaced bins over [1e-3, 1]."""
    return np.logspace(-3, 0, 51) if bins is None else np.asarray(bins, dtype=float)


def _pairs(constituents: ak.Array) -> tuple[ak.Array, ak.Array, ak.Array]:
    """Return (dR, weight-numerator, is_self) for all constituent pairs within each jet.

    weight-numerator is ``pT_i * pT_j``; normalisation is applied by the caller.
    """
    pairs = ak.argcombinations(constituents, 2, replacement=True, axis=-1)
    i, j = ak.unzip(pairs)

    ci = constituents[i]
    cj = constituents[j]
    dr = delta_r(ci, cj)

    # argcombinations gives unordered pairs; off-diagonal terms count twice
    is_self = i == j
    w = ci.pt * cj.pt
    w = ak.where(is_self, w, 2 * w)
    return dr, w, is_self


def _normalised_pairs(constituents, include_self):
    if len(constituents) == 0:
        raise ValueError("no jets supplied")
    dr, w, is_self = _pairs(constituents)

    # Normalise per jet: sum of constituent pT, squared
    sum_pt = ak.sum(constituents.pt, axis=-1)
    w = w / sum_pt**2

    if not include_self:
        # select on the pair index, not on dR > 0: two distinct particles can share
        # exactly the same (eta, phi), e.g. after packing kinematics to a fixed grid
        dr, w = dr[~is_self], w[~is_self]
    return dr, w


def eec(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
    return_errors: bool = False,
):
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
    return_errors : bool
        Also return the statistical uncertainty of each bin (see :func:`eec_per_jet`).

    Returns
    -------
    centres : np.ndarray
        Geometric bin centres.
    values : np.ndarray
        ``dSigma/dlog(R_L)`` per jet — i.e. summed weights per bin, divided by bin
        width in log space and by the number of jets.
    errors : np.ndarray
        Only with ``return_errors=True``: standard error of ``values`` from the
        jet-to-jet spread.
    """
    bins = _bins(bins)
    dr, w = _normalised_pairs(constituents, include_self)

    n_jets = len(constituents)
    dr_flat = ak.to_numpy(ak.flatten(dr, axis=None))
    w_flat = ak.to_numpy(ak.flatten(w, axis=None))

    counts, edges = np.histogram(dr_flat, bins=bins, weights=w_flat)

    dlog = np.diff(np.log(edges))
    values = counts / dlog / n_jets
    centres = np.sqrt(edges[:-1] * edges[1:])
    if not return_errors:
        return centres, values

    per_jet = _per_jet_histograms(dr, w, n_jets, edges)
    errors = per_jet.std(axis=0, ddof=1) / np.sqrt(n_jets) if n_jets > 1 \
        else np.full_like(values, np.nan)
    return centres, values, errors


def eec_per_jet(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(edges, per_jet)``, where row ``k`` is jet ``k``'s dSigma/dlog(R_L).

    ``per_jet.mean(axis=0)`` reproduces :func:`eec` (up to floating-point rounding).
    Keeping the rows is what makes honest uncertainties possible: their spread gives the
    statistical error of each bin, and re-weighting them gives a bootstrap.
    """
    bins = _bins(bins)
    dr, w = _normalised_pairs(constituents, include_self)
    return bins, _per_jet_histograms(dr, w, len(constituents), bins)


def _per_jet_histograms(dr, w, n_jets, edges) -> np.ndarray:
    n_bins = len(edges) - 1
    pairs_per_jet = ak.to_numpy(ak.num(dr, axis=-1))
    jet = np.repeat(np.arange(n_jets), pairs_per_jet)
    x = ak.to_numpy(ak.flatten(dr, axis=None))
    wt = ak.to_numpy(ak.flatten(w, axis=None))

    # same bin convention as np.histogram: [lo, hi), last bin closed on the right
    b = np.searchsorted(edges, x, side="right") - 1
    b[x == edges[-1]] = n_bins - 1
    ok = (b >= 0) & (b < n_bins)
    hist = np.bincount(jet[ok] * n_bins + b[ok], weights=wt[ok],
                       minlength=n_jets * n_bins)
    return hist.reshape(n_jets, n_bins) / np.diff(np.log(edges))


def normalisation_check(constituents: ak.Array) -> float:
    """Return the total EEC weight summed over all pairs, averaged per jet.

    Should be ~1.0 by construction. A useful unit test on any new sample: if it is
    not 1, the constituent collection or the pT normalisation is wrong.
    """
    _, w, _ = _pairs(constituents)
    sum_pt = ak.sum(constituents.pt, axis=-1)
    total = ak.sum(w, axis=-1) / sum_pt**2
    return float(ak.mean(total))
