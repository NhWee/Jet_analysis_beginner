"""Clustering: closure, ordering, field propagation and the fastjet layout trap."""

import awkward as ak
import fastjet
import numpy as np
import pytest

from jeteec import jets as jets_mod
from jeteec.kinematics import sum_p4

CFG = jets_mod.JetConfig(algorithm="antikt", R=0.4, min_pt=100.0, max_abs_eta=2.5)


@pytest.fixture(scope="module")
def clustered(toy_particles):
    return jets_mod.cluster(toy_particles, CFG)


def test_jet_is_sum_of_its_constituents(clustered):
    jets, const = clustered
    total = sum_p4(const)
    assert ak.max(abs(total.pt - jets.pt)) < 1e-9
    assert ak.max(abs(total.mass - jets.mass)) < 1e-6


def test_jets_are_pt_ordered(clustered):
    # fastjet itself returns them unordered; jets[:, 0] must be the leading jet
    jets, _ = clustered
    assert ak.all(jets.pt[:, :-1] >= jets.pt[:, 1:])


def test_cuts_are_applied(clustered):
    jets, _ = clustered
    assert ak.all(ak.flatten(jets.pt) > CFG.min_pt)
    assert ak.all(abs(ak.flatten(jets.eta)) < CFG.max_abs_eta)


def test_extra_fields_reach_the_constituents(toy_particles, clustered):
    """charge survives clustering, attached to the right particle — no dR re-matching."""
    _, const = clustered
    assert "charge" in const.fields
    for event, jets_in_event in zip(toy_particles, const):
        truth = {(p["pt"], p["charge"]) for p in ak.to_list(event)}
        for jet in ak.to_list(jets_in_event):
            assert {(c["pt"], c["charge"]) for c in jet} <= truth


def _subjet_closure(constituents):
    sub = jets_mod.subjets(constituents, algorithm="kt", R=0.2)
    total, parent = sum_p4(sub), sum_p4(constituents)
    return (float(ak.max(abs(total.pt - parent.pt))),
            float(ak.max(abs(total.mass - parent.mass))))


@pytest.mark.parametrize("reshuffle", ["acceptance cut", "reversed"])
def test_reclustering_fastjet_constituents_closes(toy_particles, reshuffle):
    """Regression test for the fastjet <= 3.5.1.4 layout trap.

    Once the constituents fastjet returns are selected (an acceptance cut, as the
    original code did) or reordered, their sublists no longer sit back to back in memory.
    fastjet <= 3.5.1.4 misread such input: reclustering it missed the parent jet by
    hundreds of GeV. Verified: this test fails on 3.5.1.4 if ``clusterable`` is a no-op.
    """
    seq = fastjet.ClusterSequence(jets_mod.clusterable(toy_particles), CFG.definition())
    jets = seq.inclusive_jets(min_pt=CFG.min_pt)
    const = seq.constituents(min_pt=CFG.min_pt)
    if reshuffle == "acceptance cut":
        raw = ak.flatten(const[abs(jets.eta) < CFG.max_abs_eta], axis=1)
    else:
        raw = ak.flatten(const, axis=1)[::-1]
    assert isinstance(raw.layout, ak.contents.ListArray)       # the dangerous layout
    d_pt, d_mass = _subjet_closure(raw)
    assert d_pt < 1e-9 and d_mass < 1e-6


def test_reclustering_cluster_output_closes(clustered):
    _, const = clustered
    cflat = ak.flatten(const, axis=1)
    d_pt, d_mass = _subjet_closure(cflat[ak.num(cflat) >= 2])
    assert d_pt < 1e-9 and d_mass < 1e-6


def test_clusterable_packs_the_layout(clustered):
    _, const = clustered
    gathered = ak.flatten(const, axis=1)[::2]          # strided: not contiguous
    packed = jets_mod.clusterable(gathered)
    assert isinstance(packed.layout, ak.contents.ListOffsetArray)
    assert isinstance(packed.layout.content, ak.contents.RecordArray)
    assert packed.layout.offsets[0] == 0
    assert packed.fields == gathered.fields


def test_four_strips_to_momentum_fields(toy_particles):
    assert jets_mod.four(toy_particles).fields == ["pt", "eta", "phi", "mass"]


def test_unknown_algorithm_is_rejected():
    with pytest.raises(ValueError, match="unknown algorithm"):
        jets_mod.JetConfig(algorithm="cone").definition()


def test_as_particles_needs_mass_or_energy():
    pt = ak.Array([[1.0]])
    with pytest.raises(ValueError):
        jets_mod.as_particles(pt, pt, pt)
    assert np.isclose(jets_mod.as_particles(pt, pt * 0, pt * 0, mass=pt * 0).pt[0, 0], 1.0)
