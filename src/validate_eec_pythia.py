"""Validate the jet + EEC pipeline on Pythia8 truth particles.

Run this BEFORE touching Open Data. With generator truth you control everything, so any
problem here is a bug in the pipeline, not in the data. Checks performed:

  1. jets are found and have sensible pT / multiplicity
  2. EEC weights sum to ~1 per jet (normalisation check)
  3. the EEC distribution has the expected shape, saved to results/figures/
  4. the EEC in jet-pT slices, with statistical errors, plotted against R_L and R_L*pT:
     the hadronisation turnover sits at R_L ~ Lambda_QCD/pT, so it should line up in R_L*pT
  5. (--compare-hadronisation) the same jets generated with hadronisation switched off,
     which shows directly that the small-angle turnover comes from hadronisation

Usage:  python src/validate_eec_pythia.py --n-events 2000 --seed 12345
        python src/validate_eec_pythia.py --n-events 5000 --seed 12345 --compare-hadronisation
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


def _import_pythia():
    """conda-forge ships the bindings as ``pythia8``; the PyPI wheel as ``pythia8mc``."""
    try:
        import pythia8
    except ImportError:
        import pythia8mc as pythia8
    return pythia8


def generate(n_events: int, seed: int, pthat_min: float, ecm: float,
             hadronisation: bool = True) -> ak.Array:
    """Generate dijet events with Pythia8, returning final-state particles per event.

    With ``hadronisation=False`` the final state is the parton shower output (quarks and
    gluons) -- unphysical as a detector input, but the cleanest way to see which part of
    the EEC is perturbative.
    """
    pythia8 = _import_pythia()

    pythia = pythia8.Pythia("", False)
    pythia.readString("HardQCD:all = on")
    pythia.readString(f"PhaseSpace:pTHatMin = {pthat_min}")
    pythia.readString(f"Beams:eCM = {ecm}")
    pythia.readString("Random:setSeed = on")
    pythia.readString(f"Random:seed = {seed}")
    pythia.readString("Print:quiet = on")
    if not hadronisation:
        pythia.readString("HadronLevel:all = off")
    pythia.init()

    events = []
    for _ in range(n_events):
        if not pythia.next():
            continue
        pt, eta, phi, m = [], [], [], []
        for p in pythia.event:
            # final-state, visible, within tracker-like acceptance
            # partons are "final" but not "visible" when hadronisation is off
            if not p.isFinal() or not (p.isVisible() or not hadronisation):
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


PT_EDGES = (100.0, 150.0, 200.0, 300.0, 500.0)


def jets_and_constituents(particles, cfg):
    jets, constituents = jets_mod.cluster(particles, cfg)
    return ak.flatten(jets), ak.flatten(constituents, axis=1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-events", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--pthat-min", type=float, default=100.0, help="GeV")
    ap.add_argument("--ecm", type=float, default=13000.0, help="GeV")
    ap.add_argument("--jet-pt-min", type=float, default=100.0, help="GeV")
    ap.add_argument("--R", type=float, default=0.4)
    ap.add_argument("--chunk-size", type=int, default=2000,
                    help="jets per EEC block (bounds memory)")
    ap.add_argument("--compare-hadronisation", action="store_true",
                    help="also generate with HadronLevel:all = off and overlay")
    args = ap.parse_args()

    np.random.seed(args.seed)
    n_steps = 5 if args.compare_hadronisation else 4

    print(f"[1/{n_steps}] generating {args.n_events} Pythia events (seed={args.seed}) ...")
    particles = generate(args.n_events, args.seed, args.pthat_min, args.ecm)
    print(f"      {len(particles)} events, "
          f"<multiplicity> = {ak.mean(ak.num(particles)):.1f}")

    print(f"[2/{n_steps}] clustering jets ...")
    cfg = jets_mod.JetConfig(algorithm="antikt", R=args.R, min_pt=args.jet_pt_min)
    jets_flat, const_flat = jets_and_constituents(particles, cfg)
    n_jets = len(jets_flat)
    if n_jets == 0:
        print("      FAIL: no jets found — check pT cuts / input units (GeV expected)")
        return 1
    print(f"      {n_jets} jets, <pT> = {ak.mean(jets_flat.pt):.1f} GeV, "
          f"<n_const> = {ak.mean(ak.num(const_flat)):.1f}")

    print(f"[3/{n_steps}] computing EEC ...")
    norm = eec_mod.normalisation_check(const_flat, chunk_size=args.chunk_size)
    print(f"      normalisation = {norm:.6f} (expect 1.0)")
    ok = abs(norm - 1.0) < 1e-6

    res = eec_mod.eec_with_errors(const_flat, chunk_size=args.chunk_size)
    peak = res.centres[np.argmax(res.values)]
    print(f"      EEC peak at R_L = {peak:.3f}")

    outdir = REPO / "results" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    plot_mod.use_style("ATLAS")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    plot_mod.plot_eec(res.centres, res.values, ax=ax,
                      label=f"Pythia8, anti-$k_T$ R={args.R}")
    ax.fill_between(res.centres, *plot_mod.log_band(res.values, res.errors), alpha=0.3)
    plot_mod.fit_log_ylim(ax, res.values)
    ax.set_title(f"EEC validation, $p_T^{{jet}}>{args.jet_pt_min:.0f}$ GeV")
    out = outdir / "eec_pythia_validation.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")

    print(f"[4/{n_steps}] EEC in jet-pT slices ...")
    slices = eec_mod.eec_in_pt_bins(const_flat, jets_flat.pt, PT_EDGES,
                                    chunk_size=args.chunk_size)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    plotted = []
    print("      pT slice [GeV]   jets   peak R_L   peak R_L*<pT> [GeV]   median rel. err")
    for lo, hi, r in slices:
        if r is None:
            print(f"      {lo:5.0f}-{hi:<5.0f}         0   (empty)")
            continue
        sel = (ak.to_numpy(jets_flat.pt) >= lo) & (ak.to_numpy(jets_flat.pt) < hi)
        mean_pt = float(np.mean(ak.to_numpy(jets_flat.pt)[sel]))
        pk = r.centres[np.argmax(r.values)]
        good = r.values > 0
        rel = float(np.median(r.errors[good] / r.values[good]))
        print(f"      {lo:5.0f}-{hi:<5.0f} {r.n_jets:9d}   {pk:8.3f}   {pk * mean_pt:12.1f}"
              f"            {rel:6.1%}")
        lab = f"{lo:.0f}-{hi:.0f} GeV"
        v, e = r.values[good], r.errors[good]
        lo_edge, _ = plot_mod.log_band(v, e)
        plotted.append(v)
        for ax, x in ((ax1, r.centres), (ax2, r.centres * mean_pt)):
            ax.errorbar(x[good], v, yerr=(v - lo_edge, e), fmt="o-", ms=3, lw=1, label=lab)
    for ax, xl in ((ax1, r"$R_L$"), (ax2, r"$R_L \cdot \langle p_T \rangle$ [GeV]")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.legend(title="jet $p_T$", fontsize=10)
    ax1.set_ylabel(r"$\frac{1}{\sigma}\frac{d\sigma}{d\log R_L}$")
    plot_mod.fit_log_ylim(ax1, *plotted)
    fig.tight_layout()
    out = outdir / "eec_pythia_pt_slices.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")

    if args.compare_hadronisation:
        print(f"[5/{n_steps}] same seed, hadronisation off ...")
        partons = generate(args.n_events, args.seed, args.pthat_min, args.ecm,
                           hadronisation=False)
        pj, pc = jets_and_constituents(partons, cfg)
        pres = eec_mod.eec_with_errors(pc, chunk_size=args.chunk_size)
        ppeak = pres.centres[np.argmax(pres.values)]
        print(f"      {len(pj)} parton-level jets, <n_const> = {ak.mean(ak.num(pc)):.1f}, "
              f"EEC peak at R_L = {ppeak:.3f}")
        # Compare SHAPES. Each EEC sums to 1 including self-pairs at dR = 0, which lie
        # outside the bins; parton-level jets have few hard constituents, so their
        # self-pair weight (sum z_i^2) is much larger and every visible bin is lower.
        # Normalising both to unit area in the binned range removes that offset.
        dlog = np.diff(np.log(res.edges))
        h_shape = res.values / np.sum(res.values * dlog)
        p_shape = pres.values / np.sum(pres.values * dlog)
        print(f"      binned EEC area (1 - self-pairs - out of range): "
              f"hadron {np.sum(res.values * dlog):.3f}, parton {np.sum(pres.values * dlog):.3f}")
        with np.errstate(divide="ignore", invalid="ignore"):
            shape_ratio = p_shape / h_shape
        big = (res.centres > 0.1) & (res.centres < args.R)
        small = (res.centres < 0.01) & (res.values > 0)
        print(f"      shape ratio parton/hadron, 0.1 < R_L < R : "
              f"{np.mean(shape_ratio[big]):.2f}  (perturbative region)")
        print(f"      shape ratio parton/hadron, R_L < 0.01    : "
              f"{np.mean(shape_ratio[small]):.2f}  (hadronisation region)")

        fig, (ax, axr) = plt.subplots(2, 1, figsize=(6, 7), sharex=True,
                                      gridspec_kw={"height_ratios": [3, 1]})
        for r, shp, lab in ((res, h_shape, "hadron level"),
                            (pres, p_shape, "parton level (no hadronisation)")):
            scale = shp[r.values > 0] / r.values[r.values > 0]
            k = float(np.median(scale)) if scale.size else 1.0
            plot_mod.plot_eec(r.centres, shp, ax=ax, label=lab, lw=2)
            lo_edge, hi_edge = plot_mod.log_band(r.values, r.errors)
            ax.fill_between(r.centres, lo_edge * k, hi_edge * k, alpha=0.3)
        plot_mod.fit_log_ylim(ax, h_shape, p_shape)
        ax.set_xlabel("")
        ax.set_ylabel("EEC, unit area in range")
        axr.plot(res.centres, shape_ratio, "k-")
        axr.axhline(1, color="grey", ls=":")
        axr.set_xscale("log")
        axr.set_ylim(0, 2)
        axr.set_xlabel(r"$R_L$")
        axr.set_ylabel("parton / hadron")
        out = outdir / "eec_pythia_hadronisation.png"
        plot_mod.save(fig, str(out))
        plt.close(fig)
        print(f"      figure -> {out.relative_to(REPO)}")

    print("\nRESULT:", "PASS" if ok else "FAIL (normalisation off)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
