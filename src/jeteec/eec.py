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
    dr, w = _delta_r(constituents)
    if norm == "sum_pt":
        denom = ak.sum(constituents.pt, axis=-1)
    elif norm == "jet_pt":
        if jet_pt is None:
            raise ValueError('norm="jet_pt" needs jet_pt (one value per jet)')
        denom = jet_pt
    else:
        raise ValueError(f'unknown norm {norm!r}; choose "sum_pt" or "jet_pt"')
    return dr, w / denom**2


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
        dr, w = _weights(block, norm, None if jet_pt is None else jet_pt[lo:hi])
        if not include_self:
            mask = dr > 0
            dr, w = dr[mask], w[mask]

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
) -> tuple[np.ndarray, np.ndarray]:
    """Jet-averaged two-point EEC as ``(centres, values)``.

    ``values`` is ``dSigma/dlog(R_L)`` per jet: summed weights per bin, divided by the bin
    width in log space and by the number of jets. Use :func:`eec_with_errors` when the
    statistical uncertainty is needed; arguments are the same.
    """
    r = eec_with_errors(constituents, bins=bins, include_self=include_self,
                        norm=norm, jet_pt=jet_pt, chunk_size=chunk_size,
                        jet_weights=jet_weights)
    return r.centres, r.values


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
        _, w = _weights(constituents[lo:hi], "sum_pt", None)
        acc += float(ak.sum(ak.sum(w, axis=-1)))
    return acc / n
