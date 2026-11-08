"""The correlator: normalisation, analytic cases, angle wrapping, uncertainties."""

import awkward as ak
import numpy as np
import pytest

from jeteec import eec as eec_mod


def make_jets(jets):
    """``jets``: list of jets, each a list of (pt, eta, phi) constituents."""
    arr = ak.Array([[{"pt": pt, "eta": eta, "phi": phi, "mass": 0.0}
                     for pt, eta, phi in jet] for jet in jets])
    return ak.zip({k: arr[k] for k in arr.fields}, with_name="Momentum4D")


@pytest.fixture(scope="module")
def random_jets():
    rng = np.random.default_rng(1)
    jets = []
    for _ in range(200):
        n = rng.integers(2, 30)
        jets.append(list(zip(rng.exponential(10.0, n), rng.normal(0, 0.1, n),
                             rng.normal(0, 0.1, n))))
    return make_jets(jets)


def test_weights_sum_to_one(random_jets):
    assert eec_mod.normalisation_check(random_jets) == pytest.approx(1.0, abs=1e-12)


def test_two_particle_jet_matches_hand_calculation():
    # z = 0.75 and 0.25 at dR = 0.1: the off-diagonal weight is 2 * 0.75 * 0.25 = 0.375
    jets = make_jets([[(75.0, 0.0, 0.0), (25.0, 0.06, 0.08)]])
    bins = np.array([0.05, 0.2, 0.5])
    _, values = eec_mod.eec(jets, bins=bins)
    assert values[0] == pytest.approx(0.375 / np.log(0.2 / 0.05))
    assert values[1] == 0.0


def test_delta_phi_wraps_around_pi():
    jets = make_jets([[(50.0, 0.0, np.pi - 0.01), (50.0, 0.0, -np.pi + 0.01)]])
    _, values = eec_mod.eec(jets, bins=np.array([0.01, 0.03, 1.0]))
    assert values[0] > 0 and values[1] == 0.0


def test_exclude_self_pairs_by_index_not_by_angle():
    # two distinct particles at the same (eta, phi), as packing to a grid can produce:
    # their pair is a real pair and must survive include_self=False
    jets = make_jets([[(40.0, 0.1, 0.2), (40.0, 0.1, 0.2), (20.0, 0.3, 0.2)]])
    _, w = eec_mod._normalised_pairs(jets, include_self=False)
    z = np.array([0.4, 0.4, 0.2])
    assert float(ak.sum(w)) == pytest.approx(1.0 - np.sum(z**2))


def test_per_jet_rows_reproduce_the_average(random_jets):
    centres, values, errors = eec_mod.eec(random_jets, return_errors=True)
    edges, per_jet = eec_mod.eec_per_jet(random_jets)
    assert per_jet.shape == (len(random_jets), len(centres))
    np.testing.assert_allclose(per_jet.mean(axis=0), values, rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(errors, per_jet.std(axis=0, ddof=1) / np.sqrt(len(per_jet)))
    np.testing.assert_allclose(np.sqrt(edges[:-1] * edges[1:]), centres)


def test_empty_input_is_an_error():
    with pytest.raises(ValueError, match="no jets"):
        eec_mod.eec(make_jets([[(1.0, 0.0, 0.0)]])[:0])
