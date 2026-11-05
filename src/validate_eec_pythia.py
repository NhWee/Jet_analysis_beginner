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

Baseline matched to the CMS 2011A MOD dataset (7 TeV, two hardest anti-kT R=0.5 jets,
pT > 375 GeV, |eta| < 1.9), with pT-hat-biased generation so the threshold is unbiased:

        python src/validate_eec_pythia.py --preset mod2011 --n-events 20000 --seed 12345 \
            --compare-hadronisation
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


def _info(pythia):
    """``pythia.info`` in the conda bindings, ``pythia.infoPython()`` in ``pythia8mc``."""
    return pythia.infoPython() if hasattr(pythia, "infoPython") else pythia.info


def generate(n_events: int, seed: int, pthat_min: float, ecm: float,
             hadronisation: bool = True, bias_pow: float = 0.0,
             bias_ref: float | None = None):
    """Generate dijet events with Pythia8.

    Returns ``(particles, weights)``: final-state particles per event and one weight per
    event, normalised to mean 1 (all 1 unless ``bias_pow > 0``).

    With ``hadronisation=False`` the final state is the parton shower output (quarks and
    gluons) -- unphysical as a detector input, but the cleanest way to see which part of
    the EEC is perturbative.

    ``bias_pow > 0`` oversamples high pT-hat by ``(pTHat/bias_ref)**bias_pow`` and returns
    compensating weights. This allows a low ``pthat_min`` (so jets just above a high pT
    threshold are not lost) without spending most events far below that threshold.
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
    if bias_pow > 0:
        pythia.readString("PhaseSpace:bias2Selection = on")
        pythia.readString(f"PhaseSpace:bias2SelectionPow = {bias_pow}")
        pythia.readString(f"PhaseSpace:bias2SelectionRef = {bias_ref or pthat_min}")
    if not pythia.init():
        # Typical cause: Pythia cannot find its XML data. The conda-forge build sets
        # PYTHIA8DATA only when the environment is *activated*.
        raise RuntimeError(
            "Pythia initialisation failed (see the PYTHIA messages above). If they mention "
            "'unmatched version numbers' or unknown keys, Pythia cannot find its data files: "
            "activate the environment (conda activate jet-eec) or set PYTHIA8DATA to "
            "<prefix>/share/Pythia8/xmldoc.")

    events, weights = [], []
    for _ in range(n_events):
        if not pythia.next():
            continue
        weights.append(_info(pythia).weight())
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
    w = np.asarray(weights, dtype=float)
    return jets_mod.as_particles(arr.pt, arr.eta, arr.phi, mass=arr.mass), w / w.mean()


PT_EDGES = (100.0, 150.0, 200.0, 300.0, 500.0)

# Settings matched to the CMS 2011A MOD HDF5 dataset (zenodo 3340205): two hardest anti-kT
# R=0.5 jets per event, pT > 375 GeV (where HLT_Jet300 is fully efficient), |eta| < 1.9.
# pTHatMin=150 with bias: a direct check found 3.3% of pT>375 jets come from pTHat<250
# and 8.5% from pTHat<300 at R=0.5, so a cut near the threshold would bias the sample.
PRESETS = {
    "mod2011": dict(ecm=7000.0, R=0.5, jet_pt_min=375.0, max_abs_eta=1.9, leading=2,
                    pthat_min=150.0, bias_pow=4.0, bias_ref=150.0,
                    pt_edges=(375.0, 450.0, 550.0, 700.0, 1000.0)),
}


def jets_and_constituents(particles, cfg, event_weights=None, leading=0):
    """Cluster, optionally keep the ``leading`` hardest jets per event before cuts.

    Returns flat (jets, constituents, jet_weights), one entry per selected jet.
    """
    if leading:
        loose = jets_mod.JetConfig(algorithm=cfg.algorithm, R=cfg.R,
                                   min_pt=min(20.0, cfg.min_pt), max_abs_eta=99.0)
        jets, const = jets_mod.cluster(particles, loose)
        jets, const = jets_mod.leading(jets, const, leading)
        keep = (jets.pt > cfg.min_pt) & (abs(jets.eta) < cfg.max_abs_eta)
        jets, const = jets[keep], const[keep]
    else:
        jets, const = jets_mod.cluster(particles, cfg)
    if event_weights is None:
        event_weights = np.ones(len(jets))
    jw = np.repeat(np.asarray(event_weights, dtype=float), ak.to_numpy(ak.num(jets)))
    return ak.flatten(jets), ak.flatten(const, axis=1), jw


def wmean(x, w):
    return float(np.average(np.asarray(ak.to_numpy(x), dtype=float), weights=w))


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
    ap.add_argument("--max-abs-eta", type=float, default=2.0)
    ap.add_argument("--leading", type=int, default=0,
                    help="keep only the N hardest jets per event, before cuts (0 = all)")
    ap.add_argument("--bias-pow", type=float, default=0.0,
                    help="pT-hat bias power for weighted generation (0 = unweighted)")
    ap.add_argument("--bias-ref", type=float, default=None, help="GeV")
    ap.add_argument("--pt-edges", type=float, nargs="+", default=list(PT_EDGES))
    ap.add_argument("--tag", default="", help="suffix for output figure names")
    ap.add_argument("--preset", choices=sorted(PRESETS),
                    help="set all of the above to a named configuration")
    args = ap.parse_args()
    if args.preset:
        for k, v in PRESETS[args.preset].items():
            setattr(args, k, v)
        args.tag = args.tag or args.preset
    suffix = f"_{args.tag}" if args.tag else ""
    print(f"      config: ecm={args.ecm:.0f} R={args.R} jet pT>{args.jet_pt_min:.0f} "
          f"|eta|<{args.max_abs_eta} leading={args.leading or 'all'} "
          f"pTHatMin={args.pthat_min:.0f} bias_pow={args.bias_pow}")

    np.random.seed(args.seed)
    n_steps = 5 if args.compare_hadronisation else 4

    print(f"[1/{n_steps}] generating {args.n_events} Pythia events (seed={args.seed}) ...")
    gen_kw = dict(bias_pow=args.bias_pow, bias_ref=args.bias_ref)
    particles, ev_w = generate(args.n_events, args.seed, args.pthat_min, args.ecm, **gen_kw)
    print(f"      {len(particles)} events, "
          f"<multiplicity> = {ak.mean(ak.num(particles)):.1f}")

    print(f"[2/{n_steps}] clustering jets ...")
    cfg = jets_mod.JetConfig(algorithm="antikt", R=args.R, min_pt=args.jet_pt_min,
                             max_abs_eta=args.max_abs_eta)
    jets_flat, const_flat, jw = jets_and_constituents(particles, cfg, ev_w, args.leading)
    n_jets = len(jets_flat)
    if n_jets == 0:
        print("      FAIL: no jets found — check pT cuts / input units (GeV expected)")
        return 1
    n_eff = jw.sum() ** 2 / (jw**2).sum()
    print(f"      {n_jets} jets (effective {n_eff:.0f}), <pT> = {wmean(jets_flat.pt, jw):.1f} GeV, "
          f"<n_const> = {wmean(ak.num(const_flat), jw):.1f}")

    print(f"[3/{n_steps}] computing EEC ...")
    norm = eec_mod.normalisation_check(const_flat, chunk_size=args.chunk_size)
    print(f"      normalisation = {norm:.6f} (expect 1.0)")
    ok = abs(norm - 1.0) < 1e-6

    res = eec_mod.eec_with_errors(const_flat, chunk_size=args.chunk_size, jet_weights=jw)
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
    out = outdir / f"eec_pythia_validation{suffix}.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")

    tables = REPO / "results" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    table = tables / f"eec_pythia{suffix}.npz"
    np.savez(table, edges=res.edges, values=res.values, errors=res.errors,
             n_jets=res.n_jets, n_eff=res.n_eff)
    print(f"      table  -> {table.relative_to(REPO)}")

    print(f"[4/{n_steps}] EEC in jet-pT slices ...")
    slices = eec_mod.eec_in_pt_bins(const_flat, jets_flat.pt, args.pt_edges,
                                    chunk_size=args.chunk_size, jet_weights=jw)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    plotted = []
    sl_lo, sl_hi, sl_pt, sl_val = [], [], [], []
    print("      pT slice [GeV]   jets  (eff)   peak R_L   peak R_L*<pT> [GeV]   median rel. err")
    for lo, hi, r in slices:
        if r is None:
            print(f"      {lo:5.0f}-{hi:<5.0f}         0   (empty)")
            continue
        sel = (ak.to_numpy(jets_flat.pt) >= lo) & (ak.to_numpy(jets_flat.pt) < hi)
        mean_pt = wmean(ak.to_numpy(jets_flat.pt)[sel], jw[sel])
        pk = r.centres[np.argmax(r.values)]
        good = r.values > 0
        rel = float(np.median(r.errors[good] / r.values[good]))
        print(f"      {lo:5.0f}-{hi:<5.0f} {r.n_jets:7d} {r.n_eff:6.0f}   {pk:8.3f}   {pk * mean_pt:12.1f}"
              f"            {rel:6.1%}")
        sl_lo.append(lo), sl_hi.append(hi), sl_pt.append(mean_pt), sl_val.append(r.values)
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
    out = outdir / f"eec_pythia_pt_slices{suffix}.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")
    # Same table again, now with the per-slice EECs (read by analyze_mod.py).
    np.savez(table, edges=res.edges, values=res.values, errors=res.errors,
             n_jets=res.n_jets, n_eff=res.n_eff, slice_lo=np.array(sl_lo),
             slice_hi=np.array(sl_hi), slice_mean_pt=np.array(sl_pt),
             slice_values=np.array(sl_val))

    if args.compare_hadronisation:
        print(f"[5/{n_steps}] same seed, hadronisation off ...")
        partons, pev_w = generate(args.n_events, args.seed, args.pthat_min, args.ecm,
                                  hadronisation=False, **gen_kw)
        pj, pc, pjw = jets_and_constituents(partons, cfg, pev_w, args.leading)
        pres = eec_mod.eec_with_errors(pc, chunk_size=args.chunk_size, jet_weights=pjw)
        ppeak = pres.centres[np.argmax(pres.values)]
        print(f"      {len(pj)} parton-level jets, <n_const> = {wmean(ak.num(pc), pjw):.1f}, "
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
        out = outdir / f"eec_pythia_hadronisation{suffix}.png"
        plot_mod.save(fig, str(out))
        plt.close(fig)
        print(f"      figure -> {out.relative_to(REPO)}")

    print("\nRESULT:", "PASS" if ok else "FAIL (normalisation off)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
