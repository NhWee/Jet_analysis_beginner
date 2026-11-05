"""Invariants of the EEC implementation. If one of these fails, no physics number from
the pipeline can be trusted."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest
import vector

from jeteec import eec

vector.register_awkward()


def jets_from_lists(*jets):
    """Build a jet-level constituent array from [(pt, eta, phi), ...] per jet."""
    return ak.zip(
        {
            "pt": [[c[0] for c in j] for j in jets],
            "eta": [[c[1] for c in j] for j in jets],
            "phi": [[c[2] for c in j] for j in jets],
            "mass": [[0.0 for _ in j] for j in jets],
        },
        with_name="Momentum4D",
    )


# ---- normalisation ---------------------------------------------------------

def test_normalisation_is_exactly_one(toy_jets):
    _, const = toy_jets
    assert eec.normalisation_check(const) == pytest.approx(1.0, abs=1e-12)


def test_normalisation_chunked_matches(toy_jets):
    _, const = toy_jets
    assert eec.normalisation_check(const, chunk_size=7) == pytest.approx(
        eec.normalisation_check(const), abs=1e-12)


def test_pair_weights_sum_to_one():
    # Every ordered pair, self-pairs included, sums to exactly 1 per jet.
    const = jets_from_lists([(10, 0, 0), (5, 0.1, 0.2), (1, -0.2, 0.05)],
                            [(2, 1.0, 1.0), (2, 1.3, 0.9)])
    _, w = eec._weights(const, "sum_pt", None)
    np.testing.assert_allclose(ak.to_numpy(ak.sum(w, axis=-1)), [1.0, 1.0])


# ---- hand-computed case ----------------------------------------------------

def test_two_particle_jet_by_hand():
    # pT 3 and 1, separated by dR = 0.1. Off-diagonal weight = 2*3*1/4^2 = 0.375.
    const = jets_from_lists([(3, 0.0, 0.0), (1, 0.1, 0.0)])
    edges = np.array([0.05, 0.2])
    c, v = eec.eec(const, bins=edges)
    assert v[0] == pytest.approx(0.375 / np.log(0.2 / 0.05))
    assert c[0] == pytest.approx(np.sqrt(0.05 * 0.2))


def test_include_self_only_changes_dr_zero():
    const = jets_from_lists([(3, 0.0, 0.0), (1, 0.1, 0.0)])
    edges = np.logspace(-3, 0, 31)
    _, a = eec.eec(const, bins=edges, include_self=True)
    _, b = eec.eec(const, bins=edges, include_self=False)
    np.testing.assert_allclose(a, b)


def test_delta_phi_wraps_across_pi():
    # phi = pi - 0.01 and -pi + 0.01 are 0.02 apart, not 2*pi - 0.02.
    const = jets_from_lists([(1, 0.0, np.pi - 0.01), (1, 0.0, -np.pi + 0.01)])
    edges = np.array([0.01, 0.03, 7.0])
    _, v = eec.eec(const, bins=edges)
    assert v[0] > 0 and v[1] == 0


# ---- chunking, errors, conventions -----------------------------------------

def test_chunked_equals_unchunked(toy_jets):
    _, const = toy_jets
    a = eec.eec_with_errors(const)
    b = eec.eec_with_errors(const, chunk_size=13)
    # equal up to floating-point summation order
    np.testing.assert_allclose(a.values, b.values, rtol=1e-9, atol=1e-13)
    np.testing.assert_allclose(a.errors, b.errors, rtol=1e-9, atol=1e-13)


def test_eec_wrapper_matches_result(toy_jets):
    _, const = toy_jets
    c, v = eec.eec(const)
    r = eec.eec_with_errors(const)
    np.testing.assert_array_equal(c, r.centres)
    np.testing.assert_array_equal(v, r.values)


def test_errors_are_sem_over_jets():
    # Two identical jets: zero spread, so zero error.
    j = [(3, 0.0, 0.0), (1, 0.1, 0.0)]
    r = eec.eec_with_errors(jets_from_lists(j, j), bins=np.array([0.05, 0.2]))
    assert r.errors[0] == pytest.approx(0.0, abs=1e-15)
    # One jet with the pair in the bin, one without: per-jet contents x and 0.
    far = [(3, 0.0, 0.0), (1, 0.5, 0.0)]
    r = eec.eec_with_errors(jets_from_lists(j, far), bins=np.array([0.05, 0.2]))
    x = 0.375 / np.log(4)
    assert r.values[0] == pytest.approx(x / 2)
    assert r.errors[0] == pytest.approx(x / 2)  # sample std x/sqrt(2), / sqrt(2)


def test_errors_shrink_with_more_jets(toy_jets):
    _, const = toy_jets
    half = eec.eec_with_errors(const[: len(const) // 4])
    full = eec.eec_with_errors(const)
    ok = (half.values > 0) & (full.values > 0)
    assert np.median(full.errors[ok] / half.errors[ok]) < 0.75


def test_jet_pt_norm(toy_jets):
    jets, const = toy_jets
    a = eec.eec(const, norm="jet_pt", jet_pt=jets.pt)[1]
    b = eec.eec(const)[1]
    assert not np.allclose(a, b)
    # same shape to within the sum_pt / jet_pt difference (a few percent)
    assert np.sum(a) / np.sum(b) == pytest.approx(1.0, rel=0.1)


def test_jet_pt_norm_requires_jet_pt(toy_jets):
    _, const = toy_jets
    with pytest.raises(ValueError, match="jet_pt"):
        eec.eec(const, norm="jet_pt")


def test_unknown_norm_rejected(toy_jets):
    _, const = toy_jets
    with pytest.raises(ValueError, match="norm"):
        eec.eec(const, norm="bogus")


def test_empty_input_rejected():
    with pytest.raises(ValueError, match="no jets"):
        eec.eec(jets_from_lists()[:0])


# ---- pT slices ---------------------------------------------------------------

def test_pt_bins_partition_jets(toy_jets):
    jets, const = toy_jets
    edges = [100, 150, 250, 2000]
    res = eec.eec_in_pt_bins(const, jets.pt, edges)
    assert [(lo, hi) for lo, hi, _ in res] == [(100, 150), (150, 250), (250, 2000)]
    n = sum(r.n_jets for _, _, r in res if r is not None)
    assert n == int(np.sum((ak.to_numpy(jets.pt) >= 100) & (ak.to_numpy(jets.pt) < 2000)))


def test_empty_pt_bin_is_none(toy_jets):
    jets, const = toy_jets
    res = eec.eec_in_pt_bins(const, jets.pt, [5000, 6000])
    assert res[0][2] is None
