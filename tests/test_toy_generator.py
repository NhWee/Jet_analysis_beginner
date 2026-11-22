"""The toy generator: reproducible, physical ranges, truth that matches the particles."""

import awkward as ak
import numpy as np

from toy_generator import ToyConfig, generate, ue_density


def test_same_seed_gives_identical_events():
    a = generate(50, seed=7)
    b = generate(50, seed=7)
    for field in a.fields:
        assert np.array_equal(ak.to_numpy(ak.flatten(a[field])),
                              ak.to_numpy(ak.flatten(b[field])))
    c = generate(50, seed=8)
    assert not np.array_equal(ak.to_numpy(ak.flatten(a.pt))[:20],
                              ak.to_numpy(ak.flatten(c.pt))[:20])


def test_phi_is_wrapped(toy_particles):
    phi = ak.to_numpy(ak.flatten(toy_particles.phi))
    assert phi.min() >= -np.pi and phi.max() < np.pi


def test_truth_does_not_change_the_particles():
    plain = generate(40, seed=11)
    particles, truth = generate(40, seed=11, return_truth=True)
    assert ak.to_list(plain) == ak.to_list(particles)
    assert ak.all((ak.num(truth) == 2) | (ak.num(truth) == 3))


def test_jet_fragments_carry_the_axis_pt():
    # the Dirichlet momentum fractions of each jet sum to one, so all particles above
    # the soft UE scale add up to the generated axes' pT, event by event
    cfg = ToyConfig()
    particles, truth = generate(30, seed=12, cfg=cfg, return_truth=True)
    total = ak.to_numpy(ak.sum(particles.pt, axis=1))
    axes = ak.to_numpy(ak.sum(truth.pt, axis=1))
    ue_max = cfg.ue_particles[1] * cfg.ue_pt_scale * 3     # generous bound on the UE pT
    assert np.all(total >= axes - 1e-9)
    assert np.all(total - axes < ue_max)


def test_ue_density_formula():
    cfg = ToyConfig()
    n_mean = (cfg.ue_particles[0] + cfg.ue_particles[1] - 1) / 2
    area = 2 * cfg.ue_eta_max * 2 * np.pi
    assert ue_density(cfg) == np.float64(n_mean * cfg.ue_pt_scale / area)
    assert ue_density(cfg, pt_min=0.3) < ue_density(cfg)
