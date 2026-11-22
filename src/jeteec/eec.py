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
* **Angle** — ``dR_ij = sqrt(dEta^2 + dPhi^2)`` with dPhi wrapped to [-pi, pi) (``kinematics.delta_r``).
* **Alternative normalisation** — ``norm="jet_pt"`` divides by the clustered jet pT squared
  instead, as CMS does (arXiv:2402.13864). The weights then no longer sum to exactly 1, so
  state which convention was used whenever a number is compared with the literature.
* **Uncertainties** — pairs inside one jet are correlated, so ``sqrt(sum w^2)`` per bin is
  wrong. The error returned by :func:`eec_with_errors` is the standard error of the mean of
  the *per-jet* histograms, which treats jets (not pairs) as the independent unit. Jets from
  the same event are still mildly correlated; for a final result use an event-level bootstrap.

References:
  Komiske, Moult, Thaler, Zhu, "Analyzing N-point energy correlators inside
  jets with CMS Open Data", Phys. Rev. Lett. 130 (2023) 051901,
  arXiv:2201.07800 -- the two-point EEC in jets, as computed here.
  Larkoski, Salam, Thaler, "Energy correlation functions for jet
  substructure", JHEP 06 (2013) 108, arXiv:1305.0007 -- the related e2/C2
  discriminants, not the R_L distribution.
"""

from __future__ import annotations

from dataclasses import dataclass

import awkward as ak
import numpy as np

from .kinematics import delta_r


def _delta_r(constituents: ak.Array) -> tuple[ak.Array, ak.Array, ak.Array]:
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
    """(dR, weight) of all pairs with ``norm="sum_pt"``; self-pairs dropped by *index*."""
    if len(constituents) == 0:
        raise ValueError("no jets supplied")
    dr, w, is_self = _weights(constituents, "sum_pt", None)
    if not include_self:
        dr, w = dr[~is_self], w[~is_self]
    return dr, w


@dataclass(frozen=True)
class EECResult:
    """A binned EEC with its statistical uncertainty.

    ``values`` and ``errors`` are ``dSigma/dlog(R_L)`` per jet; ``errors`` is the standard
    error of the mean over jets. ``n_jets`` is the number of jets averaged over and ``n_eff``
    the effective number when jets are weighted.
    """

    centres: np.ndarray
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    n_jets: int
    n_eff: float = float("nan")  # (sum w)^2 / sum w^2; equals n_jets when unweighted


def _iter_chunks(n: int, chunk_size: int | None):
    step = n if not chunk_size else int(chunk_size)
    for lo in range(0, n, max(step, 1)):
        yield lo, min(lo + step, n)


def _weights(constituents, norm, jet_pt):
    """Pair separations and normalised pair weights for a block of jets."""
    dr, w, is_self = _delta_r(constituents)
    if norm == "sum_pt":
        denom = ak.sum(constituents.pt, axis=-1)
    elif norm == "jet_pt":
        if jet_pt is None:
            raise ValueError('norm="jet_pt" needs jet_pt (one value per jet)')
        denom = jet_pt
    else:
        raise ValueError(f'unknown norm {norm!r}; choose "sum_pt" or "jet_pt"')
    return dr, w / denom**2, is_self


def eec_with_errors(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
    norm: str = "sum_pt",
    jet_pt=None,
    chunk_size: int | None = None,
    jet_weights=None,
) -> EECResult:
    """Compute the jet-averaged two-point EEC together with its statistical error.

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
    norm : {"sum_pt", "jet_pt"}
        Per-jet normalisation; see the module docstring.
    jet_pt : array-like, optional
        Clustered jet pT, one per jet. Required for ``norm="jet_pt"``.
    chunk_size : int, optional
        Process this many jets at a time. All pairs of a block are held in memory at
        once (O(sum n_i^2)), so set this for large samples. ``None`` = one block.
    jet_weights : array-like, optional
        One weight per jet, e.g. the event weight of a pT-hat-biased generator sample.
        The EEC is then the weighted mean of the per-jet histograms, and the error is the
        corresponding weighted standard error. ``None`` = unweighted (exactly as before).
    """
    if bins is None:
        bins = np.logspace(-3, 0, 51)
    edges = np.asarray(bins, dtype=float)
    nb = len(edges) - 1

    n_jets = len(constituents)
    if n_jets == 0:
        raise ValueError("no jets supplied")
    if jet_pt is not None:
        jet_pt = np.asarray(ak.to_numpy(jet_pt), dtype=float)
        if len(jet_pt) != n_jets:
            raise ValueError("jet_pt must have one entry per jet")
    if jet_weights is not None:
        jet_weights = np.asarray(ak.to_numpy(jet_weights), dtype=float)
        if len(jet_weights) != n_jets:
            raise ValueError("jet_weights must have one entry per jet")
        sw = sw2 = 0.0
        s_wh = np.zeros(nb)    # sum w h
        s_w2h = np.zeros(nb)   # sum w^2 h
        s_w2h2 = np.zeros(nb)  # sum w^2 h^2

    total = np.zeros(nb)     # sum over jets of per-jet bin content
    total_sq = np.zeros(nb)  # sum over jets of (per-jet bin content)^2

    for lo, hi in _iter_chunks(n_jets, chunk_size):
        block = constituents[lo:hi]
        dr, w, is_self = _weights(block, norm, None if jet_pt is None else jet_pt[lo:hi])
        if not include_self:
            # select on the pair index, not on dR > 0: two distinct particles can share
            # exactly the same (eta, phi), e.g. after packing kinematics to a fixed grid
            dr, w = dr[~is_self], w[~is_self]

        # Jet-summed content via np.histogram, so a single block reproduces the
        # original implementation exactly.
        dr_flat = ak.to_numpy(ak.flatten(dr, axis=None))
        w_flat = ak.to_numpy(ak.flatten(w, axis=None))
        counts, _ = np.histogram(dr_flat, bins=edges, weights=w_flat)
        total += counts

        # Per-jet content for the variance: index pairs by (jet, bin).
        jet_idx = np.repeat(np.arange(hi - lo), ak.to_numpy(ak.num(dr, axis=-1)))
        b = np.searchsorted(edges, dr_flat, side="right") - 1
        b[dr_flat == edges[-1]] = nb - 1  # np.histogram: last edge is inclusive
        ok = (b >= 0) & (b < nb)
        per_jet = np.bincount(jet_idx[ok] * nb + b[ok], weights=w_flat[ok],
                              minlength=(hi - lo) * nb).reshape(hi - lo, nb)
        total_sq += (per_jet**2).sum(axis=0)
        if jet_weights is not None:
            wj = jet_weights[lo:hi]
            sw += wj.sum()
            sw2 += (wj**2).sum()
            s_wh += wj @ per_jet
            s_w2h += (wj**2) @ per_jet
            s_w2h2 += (wj**2) @ per_jet**2

    dlog = np.diff(np.log(edges))
    centres = np.sqrt(edges[:-1] * edges[1:])
    if jet_weights is not None:
        if sw <= 0:
            raise ValueError("sum of jet_weights must be positive")
        m = s_wh / sw
        # weighted SEM: sqrt(sum w_j^2 (h_j - m)^2) / sum w_j, with the small-sample
        # factor n_eff/(n_eff-1) so that unit weights reproduce the unweighted error
        n_eff = sw**2 / sw2
        var_m = np.clip(s_w2h2 - 2 * m * s_w2h + m**2 * sw2, 0, None) / sw**2
        var_m = var_m * n_eff / (n_eff - 1) if n_eff > 1 else np.full(nb, np.nan)
        return EECResult(centres=centres, edges=edges, values=m / dlog,
                         errors=np.sqrt(var_m) / dlog, n_jets=n_jets,
                         n_eff=sw**2 / sw2)
    mean = total / n_jets
    if n_jets > 1:
        var = np.clip(total_sq / n_jets - mean**2, 0, None) * n_jets / (n_jets - 1)
        sem = np.sqrt(var / n_jets)
    else:
        sem = np.full(nb, np.nan)

    return EECResult(
        centres=np.sqrt(edges[:-1] * edges[1:]),
        edges=edges,
        values=total / dlog / n_jets,
        errors=sem / dlog,
        n_jets=n_jets,
        n_eff=float(n_jets),
    )


def eec(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
    norm: str = "sum_pt",
    jet_pt=None,
    chunk_size: int | None = None,
    jet_weights=None,
    return_errors: bool = False,
):
    """Jet-averaged two-point EEC as ``(centres, values)``, or ``(centres, values, errors)``
    with ``return_errors=True``.

    ``values`` is ``dSigma/dlog(R_L)`` per jet: summed weights per bin, divided by the bin
    width in log space and by the number of jets. Use :func:`eec_with_errors` when the
    statistical uncertainty is needed; arguments are the same.
    """
    r = eec_with_errors(constituents, bins=bins, include_self=include_self,
                        norm=norm, jet_pt=jet_pt, chunk_size=chunk_size,
                        jet_weights=jet_weights)
    if return_errors:
        return r.centres, r.values, r.errors
    return r.centres, r.values


def eec_per_jet(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    include_self: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(edges, per_jet)``, where row ``k`` is jet ``k``'s dSigma/dlog(R_L).

    ``per_jet.mean(axis=0)`` reproduces :func:`eec` (up to floating-point rounding). Same as
    :func:`per_jet_histograms` divided by the log bin widths. ``include_self`` only matters
    for a binning that reaches dR = 0.
    """
    if len(constituents) == 0:
        raise ValueError("no jets supplied")
    edges = np.logspace(-3, 0, 51) if bins is None else np.asarray(bins, dtype=float)
    if include_self:
        edges, h = per_jet_histograms(constituents, edges)
        return edges, h / np.diff(np.log(edges))
    n, nb = len(constituents), len(edges) - 1
    dr, w, is_self = _weights(constituents, "sum_pt", None)
    dr, w = dr[~is_self], w[~is_self]
    x = ak.to_numpy(ak.flatten(dr, axis=None))
    wt = ak.to_numpy(ak.flatten(w, axis=None))
    jet = np.repeat(np.arange(n), ak.to_numpy(ak.num(dr, axis=-1)))
    b = np.searchsorted(edges, x, side="right") - 1
    b[x == edges[-1]] = nb - 1
    ok = (b >= 0) & (b < nb)
    h = np.bincount(jet[ok] * nb + b[ok], weights=wt[ok], minlength=n * nb).reshape(n, nb)
    return edges, h / np.diff(np.log(edges))


def eec_in_pt_bins(
    constituents: ak.Array,
    jet_pt,
    pt_edges,
    bins: np.ndarray | None = None,
    **kwargs,
) -> list[tuple[float, float, EECResult | None]]:
    """EEC in jet-pT slices: ``[(pt_lo, pt_hi, EECResult or None), ...]``.

    The hadronisation turnover sits at ``R_L ~ Lambda_QCD / pT``, so it only becomes
    visible once jets are separated by pT. Slices with no jets give ``None``.
    Extra keyword arguments go to :func:`eec_with_errors`.
    """
    jet_pt = np.asarray(ak.to_numpy(jet_pt), dtype=float)
    if len(jet_pt) != len(constituents):
        raise ValueError("jet_pt must have one entry per jet")
    out = []
    for lo, hi in zip(pt_edges[:-1], pt_edges[1:]):
        sel = (jet_pt >= lo) & (jet_pt < hi)
        if not sel.any():
            out.append((float(lo), float(hi), None))
            continue
        kw = dict(kwargs)
        if kw.get("norm") == "jet_pt":
            kw["jet_pt"] = jet_pt[sel]
        if kw.get("jet_weights") is not None:
            kw["jet_weights"] = np.asarray(ak.to_numpy(kw["jet_weights"]), dtype=float)[sel]
        out.append((float(lo), float(hi),
                    eec_with_errors(constituents[sel], bins=bins, **kw)))
    return out


def normalisation_check(constituents: ak.Array, chunk_size: int | None = None) -> float:
    """Return the total EEC weight summed over all pairs, averaged per jet.

    Should be 1.0 by construction (``norm="sum_pt"``). A useful unit test on any new
    sample: if it is not 1, the constituent collection or the pT normalisation is wrong.
    """
    n = len(constituents)
    if n == 0:
        raise ValueError("no jets supplied")
    acc = 0.0
    for lo, hi in _iter_chunks(n, chunk_size):
        _, w, _ = _weights(constituents[lo:hi], "sum_pt", None)
        acc += float(ak.sum(ak.sum(w, axis=-1)))
    return acc / n


def per_jet_histograms(
    constituents: ak.Array,
    bins: np.ndarray | None = None,
    chunk_size: int | None = None,
    dtype=np.float64,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-jet EEC bin contents ``(edges, H)`` with ``H`` of shape ``(n_jets, n_bins)``.

    ``H[j, b]`` is the summed pair weight (``norm="sum_pt"``, self-pairs included) of jet
    ``j`` in bin ``b``, *not* divided by the bin width. Keeping the per-jet histograms lets
    several samples be combined, reweighted or bootstrapped afterwards without recomputing
    the pairs. :func:`eec_from_per_jet` turns them into the same :class:`EECResult` as
    :func:`eec_with_errors`. Memory: ``n_jets * n_bins`` numbers (use ``dtype=np.float32``
    for large samples).
    """
    if bins is None:
        bins = np.logspace(-3, 0, 51)
    edges = np.asarray(bins, dtype=float)
    nb = len(edges) - 1
    n_jets = len(constituents)
    out = np.zeros((n_jets, nb), dtype=dtype)
    for lo, hi in _iter_chunks(n_jets, chunk_size):
        dr, w, _ = _weights(constituents[lo:hi], "sum_pt", None)
        dr_flat = ak.to_numpy(ak.flatten(dr, axis=None))
        w_flat = ak.to_numpy(ak.flatten(w, axis=None))
        jet_idx = np.repeat(np.arange(hi - lo), ak.to_numpy(ak.num(dr, axis=-1)))
        b = np.searchsorted(edges, dr_flat, side="right") - 1
        b[dr_flat == edges[-1]] = nb - 1
        ok = (b >= 0) & (b < nb)
        out[lo:hi] = np.bincount(jet_idx[ok] * nb + b[ok], weights=w_flat[ok],
                                 minlength=(hi - lo) * nb).reshape(hi - lo, nb)
    return edges, out


def eec_from_per_jet(edges, per_jet, jet_weights=None) -> EECResult:
    """Weighted jet average of per-jet histograms, with the jet-level SEM.

    Same estimator and error as :func:`eec_with_errors` (weighted SEM with the
    ``n_eff/(n_eff-1)`` factor; unweighted when ``jet_weights`` is None).
    """
    edges = np.asarray(edges, dtype=float)
    h = np.asarray(per_jet, dtype=float)
    n = h.shape[0]
    if n == 0:
        raise ValueError("no jets supplied")
    w = np.ones(n) if jet_weights is None else np.asarray(jet_weights, dtype=float)
    if len(w) != n:
        raise ValueError("jet_weights must have one entry per jet")
    sw, sw2 = w.sum(), (w**2).sum()
    if sw <= 0:
        raise ValueError("sum of jet_weights must be positive")
    m = (w @ h) / sw
    n_eff = sw**2 / sw2
    var_m = (w**2) @ (h - m) ** 2 / sw**2
    var_m = var_m * n_eff / (n_eff - 1) if n_eff > 1 else np.full(len(m), np.nan)
    dlog = np.diff(np.log(edges))
    return EECResult(centres=np.sqrt(edges[:-1] * edges[1:]), edges=edges,
                     values=m / dlog, errors=np.sqrt(var_m) / dlog, n_jets=n,
                     n_eff=float(n_eff))


def shape(values: np.ndarray, edges: np.ndarray, rl_range=None) -> np.ndarray:
    """EEC normalised to unit area in ``d log R_L`` over the bins (or over ``rl_range``).

    Use this, not the raw values, to compare samples with different multiplicity. The raw
    EEC sums to 1 *including* the self-pairs at R_L = 0, which lie outside the bins and carry
    weight sum z_i^2; a sample with fewer constituents therefore has less weight in the bins
    and its raw curve is offset. ``rl_range = (lo, hi)`` normalises over bins whose centre is
    inside that range only (the common range of two samples); bins outside are still scaled.
    """
    values = np.asarray(values, dtype=float)
    edges = np.asarray(edges, dtype=float)
    dlog = np.diff(np.log(edges))
    m = np.ones(len(values), dtype=bool)
    if rl_range is not None:
        c = np.sqrt(edges[:-1] * edges[1:])
        m = (c > rl_range[0]) & (c < rl_range[1])
    area = np.sum(values[m] * dlog[m])
    if area <= 0:
        raise ValueError("EEC has no weight in the normalisation range")
    return values / area
