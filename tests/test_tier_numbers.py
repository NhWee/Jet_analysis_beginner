"""The documented numbers, and the closure tests behind them.

The README promises that re-running ``src/tier_analyses.py --n-events 2000 --seed 12345``
reproduces every number in docs/tier_analyses_results.md. This file holds the code to it.

* ``test_documented_numbers`` pins the values at the precision the docs print them. If it
  fails, either a bug crept in or the docs and figures must be regenerated with the code
  (the pinned environment, requirements-lock.txt, is what produced them).
* The other tests are closure tests: they compare measurements with what the toy
  generator put in, and should hold for any seed.
"""

import awkward as ak
import dataclasses
import numpy as np
import pytest

import tier_analyses as ta
from toy_generator import ToyConfig, generate, ue_density

N_EVENTS, SEED = 2000, 12345


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    cfg = ToyConfig()
    particles, truth = generate(N_EVENTS, seed=SEED, cfg=cfg, return_truth=True)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ta, "FIGDIR", tmp_path_factory.mktemp("figures"))
        yield {
            "nanoaod": ta.run_nanoaod(particles, truth),
            "pfnano": ta.run_pfnano(particles),
            "miniaod": ta.run_miniaod(particles, cfg),
            "aod": ta.run_aod(particles),
        }


def test_documented_numbers(results):
    nano, pf, mini, aod = (results[k] for k in ("nanoaod", "pfnano", "miniaod", "aod"))
    assert nano["n_jets"] == 3237
    assert nano["n_fit"] == pytest.approx(5.04, abs=0.005)
    assert nano["n_err"] == pytest.approx(0.07, abs=0.005)
    assert nano["n_true"] == pytest.approx(5.08, abs=0.005)
    assert nano["chi2"] == pytest.approx(21.8, abs=0.05) and nano["ndf"] == 14
    assert nano["r32"] == pytest.approx(0.151, abs=0.0005)
    assert nano["r32_err"] == pytest.approx(0.008, abs=0.0005)
    assert nano["r32_true"] == pytest.approx(0.158, abs=0.0005)

    assert pf["n_jets"] == 3237
    assert pf["n_const"] == pytest.approx(14.8, abs=0.05)
    assert pf["girth"] == pytest.approx(0.0837, abs=0.00005)
    assert pf["ptd"] == pytest.approx(0.4374, abs=0.00005)

    assert mini["gain"] == pytest.approx(3.01, abs=0.005)
    assert mini["gain_err"] == pytest.approx(0.10, abs=0.005)
    assert mini["ue_expected"] == pytest.approx(2.12, abs=0.005)
    assert mini["m_raw"] == pytest.approx(20.58, abs=0.005)
    assert mini["m_trim"] == pytest.approx(15.82, abs=0.005)
    assert mini["empty_bins"] == [0, 0, 4]
    assert mini["packing_dev_large"] < 0.02

    assert aod["n_jets"] == {"nominal": 3237, "track eff -5%": 2945,
                             "track eff -10%": 2679, "angular smear 2 mrad": 3239}
    s, s_err, l, l_err = aod["shifts"]["track eff -10%"]
    assert (round(100 * s, 1), round(100 * s_err, 1)) == (1.4, 2.1)
    assert (round(100 * l, 1), round(100 * l_err, 1)) == (-1.1, 1.2)


def test_spectrum_index_closes_against_generator(results):
    nano = results["nanoaod"]
    assert abs(nano["n_fit"] - 5.0) < 3 * nano["n_err"]
    assert abs(nano["n_fit"] - nano["n_true"]) < 2 * nano["n_err"]


def test_r32_closes_against_generator(results):
    nano = results["nanoaod"]
    assert abs(nano["r32"] - nano["r32_true"]) < 3 * nano["r32_err"]


def test_eec_is_normalised(results):
    assert results["pfnano"]["norm"] == pytest.approx(1.0, abs=1e-12)


def test_tracking_shifts_are_within_statistics(results):
    # a toy with charge assigned independently of kinematics has no real tracking
    # effect on the EEC shape to find at this sample size
    for s, s_err, l, l_err in results["aod"]["shifts"].values():
        assert abs(s) < 3 * s_err and abs(l) < 3 * l_err


def test_radius_scan_gain_is_the_generated_underlying_event():
    """Switching the UE off (same random stream) must remove ue_density x area of gain."""
    cfg_on = ToyConfig()
    cfg_off = dataclasses.replace(cfg_on, ue_pt_scale=1e-9)
    gains = []
    for cfg in (cfg_on, cfg_off):
        particles = generate(N_EVENTS, seed=SEED, cfg=cfg)
        gains.append(ta.matched_gain(particles, 0.4, 1.0, max_abs_eta=2.0).mean())
    expected = ue_density(cfg_on) * (ta.cos_weighted_area(1.0) - ta.cos_weighted_area(0.4))
    assert gains[0] - gains[1] == pytest.approx(expected, rel=0.08)


def test_cos_weighted_area_small_radius_limit():
    # 2 pi R J1(R) -> pi R^2 (1 - R^2 / 8) for small R
    for r in (0.05, 0.1):
        assert ta.cos_weighted_area(r) == pytest.approx(np.pi * r**2 * (1 - r**2 / 8),
                                                        rel=1e-6)


def test_minaod_packing_keeps_every_candidate(toy_particles):
    packed = ta.compress(toy_particles)
    assert ak.all(ak.num(packed) == ak.num(toy_particles))
    deta = ak.to_numpy(ak.flatten(abs(packed.eta - toy_particles.eta)))
    assert deta.max() <= ta.ETA_LSB / 2 + 1e-12
