"""Hand-calculable checks of the helpers used by analyze_mod_sim.py and of the per-jet EEC
histograms added to jeteec.eec."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest
import vector

pytest.importorskip("energyflow")

from analyze_mod_sim import (boot_means, n_eff, region_mean_ratio,  # noqa: E402
                             slice_weight_scale)
from jeteec import eec  # noqa: E402

vector.register_awkward()


def _jets(rows):
    return ak.zip({k: ak.Array([[p[i] for p in jet] for jet in rows])
                   for i, k in enumerate(("pt", "eta", "phi", "mass"))}, with_name="Momentum4D")


# jet 0: two particles pT 30 and 10, dR = 0.1  -> off-diagonal weight 2*30*10/40^2 = 0.375
# jet 1: three particles; pairs at dR 0.05 (pT 10,10) and 0.5, 0.55
J = _jets([[(30.0, 0.0, 0.0, 0.0), (10.0, 0.1, 0.0, 0.0)],
           [(10.0, 0.0, 0.0, 0.0), (10.0, 0.0, 0.05, 0.0), (20.0, 0.5, 0.0, 0.0)]])
EDGES = np.array([0.01, 0.07, 0.2, 1.0])


def test_per_jet_histograms_by_hand():
    edges, h = eec.per_jet_histograms(J, EDGES)
    np.testing.assert_allclose(h[0], [0.0, 0.375, 0.0])
    # jet 1: sum pT 40. (10,10) at 0.05: 2*100/1600 = 0.125; (10,20) at 0.5 and at
    # sqrt(0.5^2+0.05^2) = 0.502: 2*200/1600 each = 0.25 -> 0.5 in the last bin
    np.testing.assert_allclose(h[1], [0.125, 0.0, 0.5])


def test_per_jet_matches_eec_with_errors():
    w = np.array([2.0, 1.0])
    edges, h = eec.per_jet_histograms(J, EDGES, chunk_size=1)
    a = eec.eec_from_per_jet(edges, h, w)
    b = eec.eec_with_errors(J, EDGES, jet_weights=w)
    np.testing.assert_allclose(a.values, b.values)
    np.testing.assert_allclose(a.errors, b.errors)
    assert a.n_eff == pytest.approx(b.n_eff)
    a = eec.eec_from_per_jet(edges, h)
    b = eec.eec_with_errors(J, EDGES)
    np.testing.assert_allclose(a.values, b.values)
    np.testing.assert_allclose(a.errors, b.errors)


def test_shape_unit_area_and_common_range():
    edges = np.array([1e-3, 1e-2, 1e-1, 1.0])  # equal log widths ln(10)
    v = np.array([1.0, 2.0, 3.0])
    s = eec.shape(v, edges)
    assert np.sum(s * np.diff(np.log(edges))) == pytest.approx(1.0)
    # a sample with half the in-range weight (more self-pair weight) has the same shape
    np.testing.assert_allclose(eec.shape(0.5 * v, edges), s)
    s2 = eec.shape(v, edges, rl_range=(0.01, 1.0))  # last two bins only: area 5 ln10
    np.testing.assert_allclose(s2, v / (5 * np.log(10)))


def test_slice_weight_scale():
    # slice total 0.6 nb, two of its files read carrying 0.1 and 0.2 nb -> x2
    assert slice_weight_scale(0.6, [0.1, 0.2]) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        slice_weight_scale(1.0, [0.0])


def test_n_eff():
    assert n_eff([1, 1, 1, 1]) == pytest.approx(4.0)
    assert n_eff([3, 1]) == pytest.approx(16 / 10)


def test_region_mean_ratio():
    c = np.array([0.005, 0.05, 0.3])
    assert region_mean_ratio([1, 2, 3], [2, 2, 0], c, 0.0, 0.1) == pytest.approx(0.75)
    # division by zero is skipped, not propagated
    assert region_mean_ratio([1, 2, 3], [2, 2, 0], c, 0.0, 1.0) == pytest.approx(0.75)


def test_boot_means_shared_weights():
    rng = np.random.default_rng(1)
    h = np.tile([[1.0], [3.0]], (10, 1))  # 20 jets, so no replica has all weights zero
    a, b = boot_means([h, 2 * h], np.ones(20), 50, rng, chunk=7)
    np.testing.assert_allclose(b, 2 * a)  # same replica weights for both arrays
    assert np.all((a >= 1.0) & (a <= 3.0))
