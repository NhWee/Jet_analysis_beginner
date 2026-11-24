"""Hand-checkable tests of jeteec.matching (particle matching and pair-level response)."""

import awkward as ak
import numpy as np
import pytest

from jeteec import eec as eec_mod
from jeteec import matching


def _parts(rows):
    """rows: per jet, list of (pt, eta, phi, charge)."""
    return ak.Array([[dict(pt=a, eta=b, phi=c, charge=q) for a, b, c, q in jet] for jet in rows])


GEN = _parts([[(10.0, 0.0, 0.0, 1), (30.0, 0.1, 0.0, -1)]])
RECO = _parts([[(11.0, 0.001, 0.0, 1), (29.0, 0.1, 0.002, -1), (5.0, -0.2, 0.0, 0)]])
EDGES = np.array([0.05, 0.15, 0.25, 0.35])


def test_match_simple():
    g2r, rm = matching.match_particles(GEN, RECO)
    assert g2r.tolist() == [0, 1]
    assert rm.tolist() == [True, True, False]


def test_match_requires_same_charge_and_dr():
    gen = _parts([[(10.0, 0.0, 0.0, 1), (10.0, 0.5, 0.0, 1)]])
    reco = _parts([[(10.0, 0.0, 0.0, -1), (10.0, 0.56, 0.0, 1)]])  # wrong charge; dR 0.06
    g2r, rm = matching.match_particles(gen, reco)
    assert g2r.tolist() == [-1, -1] and not rm.any()


def test_match_is_one_to_one_closest_pt_wins():
    # both gens are within 0.05 of the single reco (pT 11.5): metrics 0.15 (A) and 0.042 (D)
    gen = _parts([[(10.0, 0.0, 0.0, 1), (12.0, 0.02, 0.0, 1)], [(5.0, 1.0, 1.0, 0)]])
    reco = _parts([[(11.5, 0.01, 0.0, 1)], [(5.0, 1.0, 1.0, 0)]])
    g2r, rm = matching.match_particles(gen, reco)
    assert g2r.tolist() == [-1, 0, 1]          # flat indices; jet 1's reco is flat index 1
    assert rm.tolist() == [True, True]


def test_pair_response_by_hand():
    out = matching.pair_response(GEN, RECO, gen_jet_bin=[0], reco_jet_bin=[0], jet_weight=[2.0],
                                 group=[1], n_groups=2, edges=EDGES, n_gen_jet_bins=1,
                                 n_reco_jet_bins=1)
    # gen pair: R = 0.1 (bin 0), w = 2 * 0.25 * 0.75 = 0.375; jet weight 2
    assert out["T"][1].tolist() == pytest.approx([0.75, 0.0, 0.0])
    assert out["T"][0].sum() == 0
    # reco z = 11/45, 29/45, 5/45: pairs A'B' (R 0.099, bin 0), A'C' (0.201, bin 1),
    # B'C' (0.300, bin 2)
    w_ab, w_ac, w_bc = 2 * 11 * 29 / 45**2, 2 * 11 * 5 / 45**2, 2 * 29 * 5 / 45**2
    assert out["R"][1] == pytest.approx(2 * np.array([w_ab, w_ac, w_bc]))
    # matched: the gen pair's image A'B', reco weight, (reco bin 0, gen bin 0)
    M = np.zeros((3, 3))
    M[0, 0] = 2 * w_ab
    assert out["M"][1] == pytest.approx(M)
    assert out["Mg"][1] == pytest.approx([0.75, 0, 0])
    # fakes = R - sum_j M = the two pairs with the unmatched photon
    assert (out["R"][1] - out["M"][1].sum(axis=1)) == pytest.approx(2 * np.array([0, w_ac, w_bc]))


def test_pair_response_jet_bins_flatten():
    out = matching.pair_response(GEN, RECO, gen_jet_bin=[1], reco_jet_bin=[0], jet_weight=[1.0],
                                 group=[0], n_groups=1, edges=EDGES, n_gen_jet_bins=2,
                                 n_reco_jet_bins=1)
    assert out["T"].shape == (1, 6) and out["T"][0, 3] == pytest.approx(0.375)
    assert out["M"][0, 0, 3] > 0 and out["M"][0, 0, :3].sum() == 0


def test_per_jet_hist_equals_eec_module():
    rng = np.random.default_rng(3)
    rows = [[(rng.exponential(5) + 0.5, rng.normal(0, 0.15), rng.normal(0, 0.15), 1)
             for _ in range(rng.integers(0, 12))] for _ in range(40)]
    parts = _parts(rows)
    edges = np.logspace(-3, 0, 31)
    import vector
    vector.register_awkward()
    four = ak.zip({"pt": parts.pt, "eta": parts.eta, "phi": parts.phi,
                   "mass": parts.pt * 0}, with_name="Momentum4D")
    keep = ak.to_numpy(ak.num(parts)) > 0
    _, h_ref = eec_mod.per_jet_histograms(four[keep], edges)
    h = matching.per_jet_hist(parts, edges)
    assert np.allclose(h[keep], h_ref) and np.all(h[~keep] == 0)


def test_loser_takes_second_choice():
    # gen A (10) and D (12) both prefer reco X (11.5); D wins (0.042 < 0.15), A falls back to
    # reco Y (pT 9, metric 0.1), which D never wanted
    gen = _parts([[(10.0, 0.0, 0.0, 1), (12.0, 0.02, 0.0, 1)]])
    reco = _parts([[(11.5, 0.01, 0.0, 1), (9.0, 0.0, 0.03, 1)]])
    g2r, rm = matching.match_particles(gen, reco)
    assert g2r.tolist() == [1, 0] and rm.all()


def test_jet_failing_reco_selection_is_a_miss_only():
    out = matching.pair_response(GEN, RECO, gen_jet_bin=[0], reco_jet_bin=[-1], jet_weight=[1.0],
                                 group=[0], n_groups=1, edges=EDGES, n_gen_jet_bins=1,
                                 n_reco_jet_bins=1)
    assert out["T"][0, 0] == pytest.approx(0.375)
    assert out["R"].sum() == 0 and out["M"].sum() == 0
