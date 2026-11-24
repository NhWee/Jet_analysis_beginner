"""Reco-gen particle matching and pair-level response histograms for the EEC.

Matching (the rule of CMS, arXiv:2402.13864): a gen particle is matched to a reco particle of
the same charge within ``dR < max_dr`` (0.05 there); among several candidates the one with
the closest pT, ``|pT_reco - pT_gen| / pT_gen``, is chosen. The match is one-to-one: if two
gen particles pick the same reco particle, the one with the smaller pT difference keeps it and
the other tries its next-best free candidate in a further round (repeated until no free gen
particle has a free candidate left).

Pairs: for every jet, all unordered pairs ``i < j`` with weight ``2 z_i z_j``,
``z = pT / sum pT`` over the particle list of that level (the ``norm="sum_pt"`` convention of
:mod:`jeteec.eec`; self-pairs sit at R_L = 0, outside any log binning, and are omitted).
A gen pair is *matched* when both of its particles are matched; its reco image is the pair of
their reco partners (a distinct reco pair, because matching is one-to-one).

Inputs are awkward arrays with fields ``pt, eta, phi, charge``, one sublist per jet.
"""

from __future__ import annotations

import awkward as ak
import numpy as np


def _flat(a, field):
    return ak.to_numpy(ak.flatten(a[field], axis=None)).astype(float)


def _dphi(a, b):
    return (a - b + np.pi) % (2 * np.pi) - np.pi


def match_particles(gen: ak.Array, reco: ak.Array, max_dr: float = 0.05):
    """One-to-one gen -> reco matching.

    Returns ``(gen_to_reco, reco_matched)`` as flat numpy arrays over all particles of all
    jets: ``gen_to_reco[k]`` is the *flat* index of the reco partner of flat gen particle
    ``k`` (or -1), ``reco_matched`` is a boolean mask over flat reco particles.
    """
    n_jets = len(gen)
    if len(reco) != n_jets:
        raise ValueError("gen and reco must have one sublist per jet each")
    ng = ak.to_numpy(ak.num(gen))
    nr = ak.to_numpy(ak.num(reco))
    g_off = np.concatenate([[0], np.cumsum(ng)])
    r_off = np.concatenate([[0], np.cumsum(nr)])
    gen_to_reco = np.full(int(ng.sum()), -1, dtype=np.int64)
    reco_matched = np.zeros(int(nr.sum()), dtype=bool)
    if gen_to_reco.size == 0 or reco_matched.size == 0:
        return gen_to_reco, reco_matched

    c = ak.cartesian({"g": gen, "r": reco}, axis=1, nested=True)
    ij = ak.argcartesian({"g": gen, "r": reco}, axis=1, nested=True)
    dr2 = (c.g.eta - c.r.eta) ** 2 + _dphi(c.g.phi, c.r.phi) ** 2
    ok = (dr2 < max_dr**2) & (c.r.charge == c.g.charge)
    metric = abs(c.r.pt - c.g.pt) / c.g.pt
    # flat list of allowed candidates (gen flat index, reco flat index, metric)
    jet = ak.broadcast_arrays(ak.local_index(ij, axis=0), ij["g"])[0]
    okf = ak.to_numpy(ak.flatten(ok, axis=None))
    jf = ak.to_numpy(ak.flatten(jet, axis=None))[okf]
    cg = ak.to_numpy(ak.flatten(ij["g"], axis=None))[okf] + g_off[jf]
    cr = ak.to_numpy(ak.flatten(ij["r"], axis=None))[okf] + r_off[jf]
    cm = ak.to_numpy(ak.flatten(metric, axis=None))[okf].astype(float)
    # rounds of mutual best choice: every free gen proposes its best free reco, every reco
    # keeps the proposal with the smallest pT difference; losers try again next round
    while cg.size:
        o = np.lexsort((cm, cg))
        first = np.ones(len(o), bool)
        first[1:] = cg[o][1:] != cg[o][:-1]
        prop = o[first]                                   # best candidate of each free gen
        o2 = np.lexsort((cm[prop], cr[prop]))
        win = np.ones(len(o2), bool)
        win[1:] = cr[prop][o2][1:] != cr[prop][o2][:-1]
        acc = prop[o2[win]]
        gen_to_reco[cg[acc]] = cr[acc]
        reco_matched[cr[acc]] = True
        free = (gen_to_reco[cg] < 0) & ~reco_matched[cr]
        cg, cr, cm = cg[free], cr[free], cm[free]
    return gen_to_reco, reco_matched


def _pairs(parts: ak.Array):
    """Flat (jet, i_flat, j_flat, R, w) of all i<j pairs; i, j are flat particle indices."""
    n = ak.to_numpy(ak.num(parts))
    off = np.concatenate([[0], np.cumsum(n)])
    pt, eta, phi = (_flat(parts, f) for f in ("pt", "eta", "phi"))
    jet_of = np.repeat(np.arange(len(n)), n)
    sum_pt = np.bincount(jet_of, weights=pt, minlength=len(n))
    sum_pt = np.where(n > 0, sum_pt, 1.0)
    z = pt / sum_pt[jet_of] if pt.size else pt
    ij = ak.argcombinations(parts.pt, 2, axis=1)
    li, lj = ak.unzip(ij)
    npair = ak.to_numpy(ak.num(li))
    jet = np.repeat(np.arange(len(n)), npair)
    i = ak.to_numpy(ak.flatten(li)).astype(np.int64) + off[jet]
    j = ak.to_numpy(ak.flatten(lj)).astype(np.int64) + off[jet]
    R = np.hypot(eta[i] - eta[j], _dphi(phi[i], phi[j]))
    w = 2.0 * z[i] * z[j]
    return jet, i, j, R, w, (pt, eta, phi, z)


def rl_bin(R, edges):
    """Bin index of R_L in ``edges`` (-1 outside), upper edge inclusive like eec.py."""
    nb = len(edges) - 1
    b = np.searchsorted(edges, R, side="right") - 1
    b[R == edges[-1]] = nb - 1
    b[(b < 0) | (b >= nb)] = -1
    return b


def per_jet_hist(parts: ak.Array, edges):
    """Per-jet EEC histogram (pairs i<j, weight 2 z_i z_j), shape ``(n_jets, n_bins)``.

    Equal to :func:`jeteec.eec.per_jet_histograms` (self-pairs fall outside the bins).
    """
    nb = len(edges) - 1
    jet, _, _, R, w, _ = _pairs(parts)
    b = rl_bin(R, edges)
    ok = b >= 0
    return np.bincount(jet[ok] * nb + b[ok], weights=w[ok],
                       minlength=len(parts) * nb).reshape(len(parts), nb)


def pair_response(gen, reco, gen_jet_bin, reco_jet_bin, jet_weight, group, n_groups, edges,
                  n_gen_jet_bins, n_reco_jet_bins, max_dr=0.05):
    """Accumulate the pair-level response of a set of jets.

    Bins are flattened as ``jet_bin * n_RL + RL_bin``. Returns a dict with

    * ``T`` (n_groups, n_gen)          -- all gen pairs, gen weight
    * ``M`` (n_groups, n_reco, n_gen)  -- matched pairs, *reco* weight, (reco bin, gen bin)
    * ``R`` (n_groups, n_reco)         -- all reco pairs, reco weight
    * ``Mg`` (n_groups, n_gen)         -- matched pairs, gen weight (for the matched fraction)

    Every jet's contribution is multiplied by ``jet_weight`` and goes to sub-sample
    ``group`` (0 .. n_groups-1). A jet with ``reco_jet_bin < 0`` (fails the reco selection)
    contributes to ``T`` only. Sub-sample so the response can be split for closure tests and
    resampled for its statistical uncertainty.
    """
    nR = len(edges) - 1
    nG, nRc = n_gen_jet_bins * nR, n_reco_jet_bins * nR
    gen_to_reco, _ = match_particles(gen, reco, max_dr)
    gj, gi, gk, gR, gw, _ = _pairs(gen)
    rj, ri, rk, rR, rw, (_, reta, rphi, rz) = _pairs(reco)
    jw = np.asarray(jet_weight, float)
    grp = np.asarray(group, np.int64)
    gjb = np.asarray(gen_jet_bin, np.int64)
    rjb = np.asarray(reco_jet_bin, np.int64)

    gb = rl_bin(gR, edges)
    gbin = np.where(gb >= 0, gjb[gj] * nR + gb, -1)
    ok = gbin >= 0
    T = np.bincount(grp[gj[ok]] * nG + gbin[ok], weights=(gw * jw[gj])[ok],
                    minlength=n_groups * nG).reshape(n_groups, nG)

    rb = rl_bin(rR, edges)
    # jets with reco_jet_bin < 0 fail the reco selection: no reco pairs (jet-level misses)
    rbin = np.where((rb >= 0) & (rjb[rj] >= 0), rjb[rj] * nR + rb, -1)
    ok = rbin >= 0
    Rh = np.bincount(grp[rj[ok]] * nRc + rbin[ok], weights=(rw * jw[rj])[ok],
                     minlength=n_groups * nRc).reshape(n_groups, nRc)

    # matched gen pairs -> their reco image
    mi, mk = gen_to_reco[gi], gen_to_reco[gk]
    m = (mi >= 0) & (mk >= 0) & (gbin >= 0)
    mi, mk, jets_m = mi[m], mk[m], gj[m]
    R_reco = np.hypot(reta[mi] - reta[mk], _dphi(rphi[mi], rphi[mk]))
    w_reco = 2.0 * rz[mi] * rz[mk]
    b = rl_bin(R_reco, edges)
    okr = (b >= 0) & (rjb[jets_m] >= 0)
    rbin_m = rjb[jets_m[okr]] * nR + b[okr]
    idx = (grp[jets_m[okr]] * nRc + rbin_m) * nG + gbin[m][okr]
    M = np.bincount(idx, weights=(w_reco * jw[jets_m])[okr],
                    minlength=n_groups * nRc * nG).reshape(n_groups, nRc, nG)
    Mg = np.bincount(grp[gj[m]] * nG + gbin[m], weights=(gw * jw[gj])[m],
                     minlength=n_groups * nG).reshape(n_groups, nG)
    return {"T": T, "M": M, "R": Rh, "Mg": Mg}
