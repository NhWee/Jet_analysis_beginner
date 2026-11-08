"""Toy event generator — a stand-in for Pythia that needs no compiled generator.

Produces dijet (sometimes trijet) events with enough structure that tier-level analyses
give recognisable shapes:

* jet pT drawn from a falling power law, dN/dpT ~ pT^-n
* dijet rapidity difference drawn flat in chi = exp|y1-y2|, which is what t-channel
  (Rutherford-like) scattering gives
* constituent multiplicity growing logarithmically with jet pT
* fragmentation from a Dirichlet, giving a steeply falling z spectrum
* angular spread ~ 1/theta inside the jet, plus a soft underlying event

Because the generator knows the answer, it can hand it back: ``generate(...,
return_truth=True)`` also returns the jet axes it fragmented, and :func:`ue_density` gives
the underlying-event pT density it put in. Analyses use these for closure tests.

**This is not physics.** Shapes are qualitatively right and the observables respond in the
right direction, but no matrix element, no parton shower, no hadronisation model is involved.
Use it to exercise and validate analysis code, never to draw a conclusion.
"""

from __future__ import annotations

from dataclasses import dataclass

import awkward as ak
import numpy as np
import vector

vector.register_awkward()

PION_MASS = 0.13957


def _wrap_phi(phi):
    """Map into [-pi, pi), as real data stores it."""
    return (phi + np.pi) % (2 * np.pi) - np.pi


@dataclass(frozen=True)
class ToyConfig:
    pt_min: float = 100.0          # GeV, minimum hard-scatter pT
    pt_max: float = 900.0          # GeV
    spectrum_index: float = 5.0    # n in dN/dpT ~ pT^-n
    chi_max: float = 16.0          # dijet chi range (flat in chi)
    third_jet_prob: float = 0.18   # ~ alpha_s suppression of an extra jet
    ue_particles: tuple = (25, 70)  # underlying-event multiplicity range
    ue_pt_scale: float = 0.75      # GeV, exponential slope of UE particles
    ue_eta_max: float = 3.0        # UE is flat in |eta| < ue_eta_max and in phi
    charged_fraction: float = 0.62  # fraction of constituents that are charged


def ue_density(cfg: ToyConfig | None = None, pt_min: float = 0.0) -> float:
    """Mean underlying-event pT per unit eta-phi area that :func:`generate` produces.

    This is the truth a radius scan or jet-area subtraction should recover. ``pt_min``
    applies a particle-level threshold: for an exponential pT spectrum with slope s,
    the pT carried above t is (t + s) exp(-t/s) per particle.
    """
    cfg = cfg or ToyConfig()
    lo, hi = cfg.ue_particles
    n_mean = (lo + hi - 1) / 2          # rng.integers(lo, hi) is uniform on lo..hi-1
    s = cfg.ue_pt_scale
    pt_per_particle = (pt_min + s) * np.exp(-pt_min / s)
    area = 2 * cfg.ue_eta_max * 2 * np.pi
    return n_mean * pt_per_particle / area


def _sample_power_law(rng, n, lo, hi, index):
    """Sample from dN/dpT ~ pT^-index on [lo, hi] by inverse CDF."""
    u = rng.uniform(0, 1, n)
    a = 1.0 - index
    return (lo**a + u * (hi**a - lo**a)) ** (1.0 / a)


def _fragment(rng, axis_eta, axis_phi, jet_pt, radius, cfg):
    """Turn one jet axis into a list of constituents."""
    # multiplicity grows like log(pT), as in real fragmentation
    n_mean = 9.0 + 7.0 * np.log(jet_pt / 50.0)
    n_const = max(4, int(rng.poisson(max(n_mean, 4.0))))

    # momentum sharing: small concentration -> steeply falling z spectrum
    z = rng.dirichlet(np.full(n_const, 0.42))

    # angular spread ~ 1/theta out to the jet radius
    u = rng.uniform(0, 1, n_const)
    dr = 0.004 * np.exp(u * np.log(radius / 0.004))
    ang = rng.uniform(0, 2 * np.pi, n_const)

    charge = (rng.uniform(0, 1, n_const) < cfg.charged_fraction).astype(np.int8)
    charge = charge * rng.choice([-1, 1], n_const).astype(np.int8)

    return {
        "pt": jet_pt * z,
        "eta": axis_eta + dr * np.cos(ang),
        "phi": _wrap_phi(axis_phi + dr * np.sin(ang)),
        "mass": np.full(n_const, PION_MASS),
        "charge": charge,
    }


def generate(n_events: int, seed: int = 12345, cfg: ToyConfig | None = None,
             radius: float = 0.4, return_truth: bool = False):
    """Generate events as flat particle lists.

    Returns an awkward array of Momentum4D records with an extra ``charge`` field,
    one sublist per event. With ``return_truth=True`` also returns the generated jet
    axes (``pt, eta, phi`` of every fragmented parton, one sublist per event) -- what a
    perfect jet reconstruction would find. Recording them draws no random numbers, so the
    particles are identical either way.
    """
    cfg = cfg or ToyConfig()
    rng = np.random.default_rng(seed)

    hard_pt = _sample_power_law(rng, n_events, cfg.pt_min, cfg.pt_max,
                                cfg.spectrum_index)
    # flat in chi = exp|dy|  ->  |dy| = ln(chi)
    chi = rng.uniform(1.0, cfg.chi_max, n_events)
    dy = np.log(chi)
    y_boost = rng.normal(0.0, 0.7, n_events)
    phi0 = rng.uniform(-np.pi, np.pi, n_events)

    events, axes = [], []
    for i in range(n_events):
        parts = {"pt": [], "eta": [], "phi": [], "mass": [], "charge": []}
        truth = {"pt": [], "eta": [], "phi": []}

        def add(d):
            for k in parts:
                parts[k].extend(np.asarray(d[k]))

        def axis(pt, eta, phi):
            truth["pt"].append(pt)
            truth["eta"].append(eta)
            truth["phi"].append(_wrap_phi(phi))

        y1 = y_boost[i] + dy[i] / 2
        y2 = y_boost[i] - dy[i] / 2
        # slight pT imbalance, as real dijets have
        pt1 = hard_pt[i]
        pt2 = hard_pt[i] * rng.uniform(0.80, 1.0)

        add(_fragment(rng, y1, phi0[i], pt1, radius, cfg))
        axis(pt1, y1, phi0[i])
        phi2 = (phi0[i] + np.pi + rng.normal(0, 0.08) + np.pi) % (2 * np.pi) - np.pi
        add(_fragment(rng, y2, phi2, pt2, radius, cfg))
        axis(pt2, y2, phi2)

        # occasional third jet, softer and at a random azimuth
        if rng.uniform() < cfg.third_jet_prob:
            pt3 = hard_pt[i] * rng.uniform(0.45, 0.85)
            y3 = rng.normal(0.0, 1.2)
            phi3 = rng.uniform(-np.pi, np.pi)
            add(_fragment(rng, y3, phi3, pt3, radius, cfg))
            axis(pt3, y3, phi3)

        # underlying event, flat in eta-phi
        n_ue = rng.integers(*cfg.ue_particles)
        ue_charge = (rng.uniform(0, 1, n_ue) < cfg.charged_fraction).astype(np.int8)
        add({
            "pt": rng.exponential(cfg.ue_pt_scale, n_ue),
            "eta": rng.uniform(-cfg.ue_eta_max, cfg.ue_eta_max, n_ue),
            "phi": rng.uniform(-np.pi, np.pi, n_ue),
            "mass": np.full(n_ue, PION_MASS),
            "charge": ue_charge * rng.choice([-1, 1], n_ue).astype(np.int8),
        })

        events.append(parts)
        axes.append(truth)

    arr = ak.Array(events)
    particles = ak.zip(
        {"pt": arr.pt, "eta": arr.eta, "phi": arr.phi,
         "mass": arr.mass, "charge": arr.charge},
        with_name="Momentum4D",
    )
    if not return_truth:
        return particles
    t = ak.Array(axes)
    truth_axes = ak.zip({"pt": t.pt, "eta": t.eta, "phi": t.phi,
                         "mass": ak.zeros_like(t.pt)}, with_name="Momentum4D")
    return particles, truth_axes
