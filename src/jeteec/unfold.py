"""Unfolding of binned, *weighted* distributions (D'Agostini iterative Bayesian).

The EEC is a pair-weighted observable: every jet contributes many pairs, each with a weight
``z_i z_j``. A reco pair and a gen pair are not the same object, so a response "per event"
does not exist. The forward model used here (see docs/mod_unfolding_results.md) is linear in
the *weighted* truth histogram ``t``:

    reco_i  =  sum_j A_ij t_j  +  fake_i,         A_ij = M_ij / T_j

* ``T_j``  -- summed gen-pair weight in gen bin ``j`` (all gen pairs, matched or not);
* ``M_ij`` -- summed *reco* weight of reco pairs whose two particles are both matched to the
  two particles of a gen pair, in reco bin ``i`` and gen bin ``j``;
* ``fake_i`` -- reco weight of reco pairs with at least one unmatched particle (or whose
  matched gen pair falls outside the gen binning).

``eff_j = sum_i A_ij`` is then a weighted "efficiency" that absorbs unmatched gen pairs, pairs
migrating out of range *and* the mean reco/gen weight ratio. D'Agostini's iteration only needs
this linear expectation, not a probability per event, so columns need not sum to <= 1.

Identity property (tested): if the data are exactly ``A t + fake`` and the prior is ``t``,
every iteration returns ``t``. With a diagonal ``A`` the first iteration returns the
efficiency-corrected data whatever the prior.

Reference for the algorithm: G. D'Agostini, Nucl. Instrum. Meth. A 362 (1995) 487; used for the
EEC by CMS, arXiv:2402.13864 (RooUnfold, stopping when the p-value between consecutive
iterations exceeds 0.05).
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def response_from_counts(M, T):
    """``A_ij = M_ij / T_j`` (0 where ``T_j == 0``)."""
    M = np.asarray(M, dtype=float)
    T = np.asarray(T, dtype=float)
    A = np.zeros_like(M)
    ok = T > 0
    A[:, ok] = M[:, ok] / T[ok]
    return A


def fold(A, t, fake=None):
    """Expected reco histogram ``A @ t (+ fake)``."""
    r = np.asarray(A, float) @ np.asarray(t, float)
    return r if fake is None else r + np.asarray(fake, float)


def dagostini(data, A, prior, n_iter, fake=None, fake_fraction=None, return_all=False):
    """D'Agostini iterative unfolding of a (weighted) histogram.

    Parameters
    ----------
    data : reco histogram (1-D; flatten multi-dimensional binnings first).
    A : response, shape ``(n_reco, n_gen)``, ``reco = A @ truth`` for the matched part.
    prior : starting truth (iteration 0), shape ``(n_gen,)``.
    n_iter : number of iterations (>= 1).
    fake : absolute fake contribution to subtract, or
    fake_fraction : fraction of each reco bin that is fake (``data * (1 - f)`` is unfolded).
    return_all : also return the list of all iterations ``[t_1, ..., t_n]``.

    Gen bins with zero efficiency stay at 0 (no information).
    """
    d = np.asarray(data, dtype=float).copy()
    A = np.asarray(A, dtype=float)
    if fake is not None and fake_fraction is not None:
        raise ValueError("give fake or fake_fraction, not both")
    if fake is not None:
        d = d - np.asarray(fake, float)
    if fake_fraction is not None:
        d = d * (1.0 - np.asarray(fake_fraction, float))
    if n_iter < 1:
        raise ValueError("n_iter must be >= 1")
    eff = A.sum(axis=0)
    t = np.asarray(prior, dtype=float).copy()
    out = []
    has = eff > 0
    for _ in range(n_iter):
        nu = A @ t
        ratio = np.zeros_like(d)
        okr = nu > 0
        ratio[okr] = d[okr] / nu[okr]
        new = np.zeros_like(t)
        new[has] = t[has] * (A[:, has].T @ ratio)[...] / eff[has]
        t = new
        out.append(t.copy())
    return (t, out) if return_all else t


def bin_by_bin(data, reco_sim, truth_sim):
    """Bin-by-bin correction ``data * truth_sim / reco_sim`` (same binning on both sides)."""
    data, reco_sim, truth_sim = (np.asarray(a, float) for a in (data, reco_sim, truth_sim))
    c = np.zeros_like(data)
    ok = reco_sim > 0
    c[ok] = truth_sim[ok] / reco_sim[ok]
    return data * c, c


def iteration_pvalue(t_new, t_old, sigma):
    """p-value of chi2 = sum((t_new - t_old)/sigma)^2 with ndf = number of bins used.

    Stopping rule of arXiv:2402.13864: stop once this exceeds 0.05 (the change from one more
    iteration is statistically insignificant). Bins with ``sigma <= 0`` are skipped.
    """
    t_new, t_old, sigma = (np.asarray(a, float) for a in (t_new, t_old, sigma))
    m = sigma > 0
    chi2 = float(np.sum(((t_new[m] - t_old[m]) / sigma[m]) ** 2))
    return float(stats.chi2.sf(chi2, int(m.sum()))), chi2
