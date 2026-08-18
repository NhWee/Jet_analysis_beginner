"""Validate the jet + EEC pipeline on Pythia8 truth particles.

Run this BEFORE touching Open Data. With generator truth you control everything, so any
problem here is a bug in the pipeline, not in the data. Checks performed:

  1. jets are found and have sensible pT / multiplicity
  2. EEC weights sum to ~1 per jet (normalisation check)
  3. the EEC distribution has the expected shape, saved to results/figures/

Usage:  python src/validate_eec_pythia.py --n-events 2000 --seed 12345
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from jeteec import eec as eec_mod
from jeteec import jets as jets_mod
from jeteec import plot as plot_mod

REPO = pathlib.Path(__file__).resolve().parents[1]


def generate(n_events: int, seed: int, pthat_min: float, ecm: float) -> ak.Array:
    """Generate dijet events with Pythia8, returning final-state particles per event."""
    import pythia8

    pythia = pythia8.Pythia("", False)
    pythia.readString("HardQCD:all = on")
    pythia.readString(f"PhaseSpace:pTHatMin = {pthat_min}")
    pythia.readString(f"Beams:eCM = {ecm}")
    pythia.readString("Random:setSeed = on")
    pythia.readString(f"Random:seed = {seed}")
    pythia.readString("Print:quiet = on")
    pythia.init()

    events = []
    for _ in range(n_events):
        if not pythia.next():
            continue
        pt, eta, phi, m = [], [], [], []
        for p in pythia.event:
            # final-state, visible, within tracker-like acceptance
            if not p.isFinal() or not p.isVisible():
                continue
            if abs(p.eta()) > 4.0 or p.pT() < 0.2:
                continue
            pt.append(p.pT())
            eta.append(p.eta())
            phi.append(p.phi())
            m.append(p.m())
        events.append({"pt": pt, "eta": eta, "phi": phi, "mass": m})

    arr = ak.Array(events)
    return jets_mod.as_particles(arr.pt, arr.eta, arr.phi, mass=arr.mass)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-events", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--pthat-min", type=float, default=100.0, help="GeV")
    ap.add_argument("--ecm", type=float, default=13000.0, help="GeV")
    ap.add_argument("--jet-pt-min", type=float, default=100.0, help="GeV")
    ap.add_argument("--R", type=float, default=0.4)
    args = ap.parse_args()

    np.random.seed(args.seed)

    print(f"[1/3] generating {args.n_events} Pythia events (seed={args.seed}) ...")
    particles = generate(args.n_events, args.seed, args.pthat_min, args.ecm)
    print(f"      {len(particles)} events, "
          f"<multiplicity> = {ak.mean(ak.num(particles)):.1f}")

    print("[2/3] clustering jets ...")
    cfg = jets_mod.JetConfig(algorithm="antikt", R=args.R, min_pt=args.jet_pt_min)
    jets, constituents = jets_mod.cluster(particles, cfg)

    jets_flat = ak.flatten(jets)
    const_flat = ak.flatten(constituents, axis=1)
    n_jets = len(jets_flat)
    print(f"      {n_jets} jets, <pT> = {ak.mean(jets_flat.pt):.1f} GeV, "
          f"<n_const> = {ak.mean(ak.num(const_flat)):.1f}")
    if n_jets == 0:
        print("      FAIL: no jets found — check pT cuts / input units (GeV expected)")
        return 1

    print("[3/3] computing EEC ...")
    norm = eec_mod.normalisation_check(const_flat)
    print(f"      normalisation = {norm:.6f} (expect 1.0)")
    ok = abs(norm - 1.0) < 1e-6

    centres, values = eec_mod.eec(const_flat)

    outdir = REPO / "results" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    plot_mod.use_style("ATLAS")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    plot_mod.plot_eec(centres, values, label=f"Pythia8, anti-$k_T$ R={args.R}", ax=ax)
    ax.set_title(f"EEC validation, $p_T^{{jet}}>{args.jet_pt_min:.0f}$ GeV")
    out = outdir / "eec_pythia_validation.png"
    plot_mod.save(fig, str(out))
    print(f"      figure -> {out.relative_to(REPO)}")

    print("\nRESULT:", "PASS" if ok else "FAIL (normalisation off)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
