"""Clustering wrapper, and the guard against the fastjet extra-field trap."""

from __future__ import annotations

import awkward as ak
import numpy as np
import pytest

from jeteec import jets as jets_mod


def test_four_strips_extra_fields(toy_particles):
    assert "charge" in toy_particles.fields
    assert sorted(jets_mod.four(toy_particles).fields) == ["eta", "mass", "phi", "pt"]


def test_cluster_ignores_extra_fields(toy_particles):
    """Passing records with ``charge`` must give the same jets as stripped records.

    Regression test for the silent-wrong-jets trap (README: up to 168 GeV pT error)."""
    j1, c1 = jets_mod.cluster(toy_particles)
    j2, c2 = jets_mod.cluster(jets_mod.four(toy_particles))
    np.testing.assert_allclose(ak.to_numpy(ak.flatten(j1.pt)), ak.to_numpy(ak.flatten(j2.pt)))
    assert ak.all(ak.num(c1, axis=2) == ak.num(c2, axis=2))


def test_constituents_sum_to_jet(toy_jets):
    jets, const = toy_jets
    s = ak.sum(const, axis=-1)
    np.testing.assert_allclose(ak.to_numpy(s.pt), ak.to_numpy(jets.pt), rtol=1e-6)
    np.testing.assert_allclose(ak.to_numpy(s.mass), ak.to_numpy(jets.mass), rtol=1e-4, atol=1e-3)


def test_cuts_applied(toy_jets):
    jets, _ = toy_jets
    cfg = jets_mod.JetConfig()
    assert ak.all(jets.pt >= cfg.min_pt)
    assert ak.all(abs(jets.eta) < cfg.max_abs_eta)


def test_unknown_algorithm():
    with pytest.raises(ValueError, match="unknown algorithm"):
        jets_mod.JetConfig(algorithm="siscone").definition()


def test_as_particles_needs_mass_or_energy():
    with pytest.raises(ValueError):
        jets_mod.as_particles([[1.0]], [[0.0]], [[0.0]])
