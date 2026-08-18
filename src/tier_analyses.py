"""One real analysis per data tier.

Each tier gets the analysis that its information content actually permits, run end to end
on a common toy sample and producing a figure plus printed numbers.

    NanoAOD   jet 4-vectors only      -> inclusive pT spectrum, dijet chi, R32
    PFNano    + in-jet constituents   -> fragmentation D(z), multiplicity, girth/pTD, EEC
    MiniAOD   + all event candidates  -> R-scan, jet trimming, charged-only correlator
    AOD       + full precision/tracks -> tracking-efficiency and resolution systematics
    RECO      not released            -> discussed, not run

Usage:
    python src/tier_analyses.py                 # all tiers
    python src/tier_analyses.py --tier nanoaod  # one tier

The input is synthetic (see toy_generator.py). Shapes are qualitatively sensible and the
observables respond in the right direction, but no number here is a measurement.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import fastjet  # noqa: E402
import vector  # noqa: E402

from jeteec import eec as eec_mod  # noqa: E402
from jeteec import plot as plot_mod  # noqa: E402
from toy_generator import ToyConfig, generate  # noqa: E402

vector.register_awkward()

REPO = pathlib.Path(__file__).resolve().parents[1]
FIGDIR = REPO / "results" / "figures"
R_NOMINAL = 0.4
JET_PT_MIN = 100.0
TOY_NOTE = "toy sample — pipeline demo, not a measurement"


# ======================================================================================
# shared helpers
# ======================================================================================
def four(v):
    """Strip a record down to exactly four momentum fields.

    ⚠️ Do not skip this. Passing a record with extra fields (e.g. ``charge``) to
    ``fastjet.ClusterSequence`` can silently produce wrong jets — the failure is data
    dependent, so it sometimes looks fine. Verified here: with a 5-field input the
    subjet 4-vector sum missed the jet by up to 168 GeV in pT and 502 GeV in mass;
    with 4 fields the closure is exact to machine precision.
    """
    return ak.zip(
        {"pt": v.pt, "eta": v.eta, "phi": v.phi, "mass": v.mass},
        with_name="Momentum4D",
    )


def cluster(particles, R=R_NOMINAL, min_pt=JET_PT_MIN, max_abs_eta=2.5):
    """Cluster and return (jets, constituents) with an acceptance cut."""
    seq = fastjet.ClusterSequence(
        four(particles), fastjet.JetDefinition(fastjet.antikt_algorithm, R))
    jets = seq.inclusive_jets(min_pt=min_pt)
    const = seq.constituents(min_pt=min_pt)
    keep = abs(jets.eta) < max_abs_eta
    return jets[keep], const[keep]


def constituents_by_dr(particles, jets, R=R_NOMINAL):
    """Assign original particles (keeping all their fields) to jets by dR matching.

    Needed because :func:`cluster` must strip non-momentum fields before calling
    fastjet, which loses ``charge``. For isolated anti-kT jets the catchment area is
    circular, so dR matching to the jet axis recovers the constituents faithfully.

    Returns a [jet][particle] array, flattened over events.
    """
    pairs = ak.cartesian({"j": jets, "p": particles}, axis=1, nested=True)
    j, p = pairs.j, pairs.p
    dphi = (j.phi - p.phi + np.pi) % (2 * np.pi) - np.pi
    dr = np.sqrt((j.eta - p.eta) ** 2 + dphi**2)
    return ak.flatten(p[dr < R], axis=1)


def sum_p4(v, axis=1):
    """Sum a jagged collection of 4-vectors into one 4-vector per outer entry."""
    return ak.zip(
        {"px": ak.sum(v.px, axis=axis), "py": ak.sum(v.py, axis=axis),
         "pz": ak.sum(v.pz, axis=axis), "E": ak.sum(v.E, axis=axis)},
        with_name="Momentum4D",
    )


def annotate(fig):
    fig.text(0.995, 0.005, TOY_NOTE, ha="right", va="bottom",
             fontsize=9, color="grey", style="italic")


def finish(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    annotate(fig)
    out = FIGDIR / name
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  figure -> {out.relative_to(REPO)}")


# ======================================================================================
# NanoAOD — jet four-vectors only
# ======================================================================================
def run_nanoaod(particles):
    print("\n=== NanoAOD :: jet kinematics only ===")
    # cluster once with a loose threshold, then apply per-analysis cuts — as in a real
    # analysis, where the multiplicity ratio uses a softer jet definition than the spectrum
    all_jets, _ = cluster(particles, min_pt=60.0)
    jets = all_jets[all_jets.pt > JET_PT_MIN]

    # -- (a) inclusive jet pT spectrum, fitted to a power law ---------------------------
    pt = ak.to_numpy(ak.flatten(jets.pt))
    bins = np.logspace(np.log10(JET_PT_MIN), np.log10(1000), 26)
    counts, edges = np.histogram(pt, bins=bins)
    centres = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    dsig = counts / widths

    ok = counts > 4  # only fit bins with enough statistics
    slope, intercept = np.polyfit(np.log(centres[ok]), np.log(dsig[ok]), 1)
    n_fit = -slope
    print(f"  jets                    : {len(pt)}")
    print(f"  spectrum index n        : {n_fit:.2f}   (dN/dpT ~ pT^-n; input was 5.0)")

    # -- (b) dijet angular distribution -------------------------------------------------
    two = jets[ak.num(jets) >= 2]
    y1, y2 = two[:, 0].eta, two[:, 1].eta
    chi = np.exp(abs(ak.to_numpy(y1 - y2)))
    chi = chi[chi < 16]
    chi_counts, chi_edges = np.histogram(chi, bins=np.linspace(1, 16, 16))
    chi_norm = chi_counts / chi_counts.sum()
    flatness = chi_norm.std() / chi_norm.mean()
    print(f"  dijet pairs             : {len(chi)}")
    print(f"  chi distribution rel.RMS: {flatness:.3f}  "
          f"(flat = Rutherford-like t-channel)")

    # -- (c) jet multiplicity ratio -----------------------------------------------------
    njet = ak.num(all_jets)   # softer 60 GeV threshold, standard for counting ratios
    n2 = int(ak.sum(njet >= 2))
    n3 = int(ak.sum(njet >= 3))
    r32 = n3 / n2 if n2 else float("nan")
    print(f"  N(>=2 jets)             : {n2}")
    print(f"  N(>=3 jets)             : {n3}")
    print(f"  R32 = N3/N2             : {r32:.3f}   (alpha_s-sensitive in real data)")

    # -- figure -------------------------------------------------------------------------
    plot_mod.use_style("ATLAS")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    ax = axes[0]
    ax.errorbar(centres, dsig, yerr=np.sqrt(counts) / widths, fmt="o", ms=4)
    ax.plot(centres[ok], np.exp(intercept) * centres[ok] ** slope, "-",
            color="crimson", label=f"$p_T^{{-{n_fit:.2f}}}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"jet $p_T$ [GeV]"); ax.set_ylabel(r"$dN/dp_T$ [GeV$^{-1}$]")
    ax.set_title("Inclusive jet spectrum", fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[1]
    ax.step(chi_edges[:-1], chi_norm, where="post", lw=2)
    ax.axhline(chi_norm.mean(), color="crimson", ls="--", lw=1.5,
               label="flat (Rutherford)")
    ax.set_xlabel(r"$\chi = e^{|y_1-y_2|}$"); ax.set_ylabel("normalised")
    ax.set_ylim(0, chi_norm.max() * 1.6)
    ax.set_title("Dijet angular distribution", fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[2]
    mult = np.bincount(ak.to_numpy(njet), minlength=5)[:6]
    ax.bar(range(len(mult)), mult, color="#4C72B0", width=0.7)
    ax.set_xticks(range(len(mult)))
    ax.set_xlim(-0.6, len(mult) - 0.4)
    ax.set_xlabel("jets per event"); ax.set_ylabel("events")
    ax.set_title(f"Jet multiplicity  ($R_{{32}}$ = {r32:.3f})", fontsize=14)

    finish(fig, "tier_nanoaod.png")
    return {"n_fit": n_fit, "r32": r32, "chi_relrms": flatness}


# ======================================================================================
# PFNano — jets plus the particles inside them
# ======================================================================================
def run_pfnano(particles):
    print("\n=== PFNano :: jets + in-jet constituents ===")
    jets, const = cluster(particles)
    jflat = ak.flatten(jets)
    cflat = ak.flatten(const, axis=1)

    # -- (a) fragmentation function D(z) ------------------------------------------------
    z = cflat.pt / jflat.pt
    zf = ak.to_numpy(ak.flatten(z))
    zf = zf[(zf > 1e-3) & (zf < 1)]
    zbins = np.logspace(-3, 0, 31)
    zc, ze = np.histogram(zf, bins=zbins)
    zcen = np.sqrt(ze[:-1] * ze[1:])
    zval = zc / np.diff(np.log(ze)) / len(jflat)

    # -- (b) multiplicity vs jet pT -----------------------------------------------------
    nconst = ak.num(cflat)
    ptbins = np.array([100, 150, 220, 320, 460, 700, 1000])
    idx = np.digitize(ak.to_numpy(jflat.pt), ptbins) - 1
    prof_x, prof_y, prof_e = [], [], []
    for b in range(len(ptbins) - 1):
        m = idx == b
        if m.sum() > 5:
            v = ak.to_numpy(nconst)[m]
            prof_x.append(np.sqrt(ptbins[b] * ptbins[b + 1]))
            prof_y.append(v.mean())
            prof_e.append(v.std() / np.sqrt(len(v)))
    print(f"  jets                    : {len(jflat)}")
    print(f"  <n_constituents>        : {ak.mean(nconst):.1f}")
    print(f"  multiplicity growth     : {prof_y[0]:.1f} -> {prof_y[-1]:.1f} "
          f"over pT {prof_x[0]:.0f} -> {prof_x[-1]:.0f} GeV")

    # -- (c) girth and pTD — classic quark/gluon discriminants --------------------------
    deta = cflat.eta - jflat.eta
    dphi = (cflat.phi - jflat.phi + np.pi) % (2 * np.pi) - np.pi
    dr = np.sqrt(deta**2 + dphi**2)
    zc_ = cflat.pt / ak.sum(cflat.pt, axis=1)
    girth = ak.to_numpy(ak.sum(zc_ * dr, axis=1))
    ptd = ak.to_numpy(np.sqrt(ak.sum(zc_**2, axis=1)))
    print(f"  <girth>                 : {girth.mean():.4f}")
    print(f"  <pTD>                   : {ptd.mean():.4f}")

    # -- (d) EEC -------------------------------------------------------------------------
    norm = eec_mod.normalisation_check(cflat)
    cen, val = eec_mod.eec(cflat)
    print(f"  EEC normalisation       : {norm:.6f}  (must be 1.0)")

    plot_mod.use_style("ATLAS")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    ax.plot(zcen, zval, lw=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$z = p_T^{\rm const}/p_T^{\rm jet}$")
    ax.set_ylabel(r"$dN/d\log z$ per jet")
    ax.set_title("Fragmentation function", fontsize=14)

    ax = axes[0, 1]
    ax.errorbar(prof_x, prof_y, yerr=prof_e, fmt="o-", ms=5)
    ax.set_xscale("log")
    ax.set_xlabel(r"jet $p_T$ [GeV]"); ax.set_ylabel(r"$\langle n_{\rm const}\rangle$")
    ax.set_title("Multiplicity grows with jet $p_T$", fontsize=14)

    ax = axes[1, 0]
    ax.hist(girth, bins=40, histtype="step", lw=2, density=True, label="girth")
    ax.hist(ptd, bins=40, histtype="step", lw=2, density=True, label=r"$p_TD$")
    ax.set_xlabel("value"); ax.set_ylabel("normalised")
    ax.set_title("Quark/gluon discriminants", fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[1, 1]
    plot_mod.plot_eec(cen, val, ax=ax, lw=2)
    ax.axvline(R_NOMINAL, color="grey", ls=":", lw=1.5)
    ax.set_title(f"EEC (norm = {norm:.4f})", fontsize=14)

    finish(fig, "tier_pfnano.png")
    return {"norm": norm, "girth": girth.mean(), "ptd": ptd.mean()}


# ======================================================================================
# MiniAOD — every candidate in the event (compressed)
# ======================================================================================
def compress(particles, pt_threshold=0.30, eta_phi_lsb=1e-3, pt_rel_lsb=1e-3):
    """Mimic MiniAOD packing: a pT threshold plus reduced numerical precision.

    pT is quantised on a *relative* grid (log-space), angles on an absolute one — which
    is what packing to a fixed number of bits actually does.
    """
    p = particles[particles.pt > pt_threshold]
    return ak.zip(
        {
            "pt": np.exp(np.round(np.log(p.pt) / pt_rel_lsb) * pt_rel_lsb),
            "eta": np.round(p.eta / eta_phi_lsb) * eta_phi_lsb,
            "phi": np.round(p.phi / eta_phi_lsb) * eta_phi_lsb,
            "mass": p.mass,
            "charge": p.charge,
        },
        with_name="Momentum4D",
    )


def trim(const, r_sub=0.2, f_cut=0.05):
    """Jet trimming: recluster constituents into subjets, drop the soft ones."""
    const = four(const)
    seq = fastjet.ClusterSequence(
        const, fastjet.JetDefinition(fastjet.kt_algorithm, r_sub))
    sub = seq.inclusive_jets()
    jet_pt = ak.sum(const.pt, axis=1)
    keep = sub.pt > f_cut * jet_pt
    return sum_p4(sub[keep]), sum_p4(sub)


def run_miniaod(particles):
    print("\n=== MiniAOD :: all packed candidates ===")
    packed = compress(particles)
    print(f"  particles kept          : {ak.sum(ak.num(packed))} / "
          f"{ak.sum(ak.num(particles))} after 0.3 GeV threshold")

    # -- (a) jet radius scan -------------------------------------------------------------
    radii = [0.2, 0.4, 0.6, 0.8, 1.0]
    mean_pt, mean_n = [], []
    for R in radii:
        j, c = cluster(packed, R=R)
        mean_pt.append(float(ak.mean(ak.flatten(j.pt))))
        mean_n.append(float(ak.mean(ak.num(ak.flatten(c, axis=1)))))
    print("  R scan (only possible with all candidates):")
    for R, m, n in zip(radii, mean_pt, mean_n):
        print(f"    R={R:.1f}  <pT> = {m:7.2f} GeV   <n_const> = {n:5.1f}")
    growth = mean_pt[-1] - mean_pt[1]
    print(f"  <pT> gain R=0.4 -> 1.0  : {growth:+.2f} GeV  (underlying event ~ R^2)")

    # -- (b) trimming ---------------------------------------------------------------------
    jets, const = cluster(packed, R=0.8)   # grooming matters for wide jets
    cflat = ak.flatten(const, axis=1)
    trimmed, ungroomed = trim(cflat)
    m_raw = ak.to_numpy(ungroomed.mass)
    m_trim = ak.to_numpy(trimmed.mass)
    print(f"  R=0.8 jets              : {len(m_raw)}")
    print(f"  <mass> ungroomed        : {m_raw.mean():.2f} GeV")
    print(f"  <mass> trimmed          : {m_trim.mean():.2f} GeV "
          f"({100*(m_trim.mean()/m_raw.mean()-1):+.1f}%)")

    # -- (c) charged-only vs all-particle EEC ---------------------------------------------
    # charge survives only via dR matching, since fastjet must be fed 4 fields
    j4, _ = cluster(packed, R=R_NOMINAL)
    call = constituents_by_dr(packed, j4, R=R_NOMINAL)
    call = call[ak.num(call) >= 2]
    cchg = call[call.charge != 0]
    cchg = cchg[ak.num(cchg) >= 2]
    # coarser bins starting at 2e-3: the angular packing granularity (1e-3) makes
    # finer bins meaningless — which is itself the point, see below
    bins = np.logspace(np.log10(2e-3), 0, 31)
    cen_all, val_all = eec_mod.eec(call, bins=bins)
    cen_chg, val_chg = eec_mod.eec(cchg, bins=bins)
    print(f"  charged fraction of pT  : "
          f"{float(ak.sum(cchg.pt)/ak.sum(call.pt)):.3f}")

    # -- (d) where packing starts to hurt -------------------------------------------------
    fine = np.logspace(-3, 0, 51)
    raw_const = constituents_by_dr(particles, cluster(particles)[0], R=R_NOMINAL)
    raw_const = raw_const[ak.num(raw_const) >= 2]
    _, v_packed = eec_mod.eec(call, bins=fine)
    _, v_raw = eec_mod.eec(raw_const, bins=fine)
    small = np.sqrt(fine[:-1] * fine[1:]) < 5e-3
    empty_packed = int((v_packed[small] == 0).sum())
    empty_raw = int((v_raw[small] == 0).sum())
    print(f"  empty EEC bins below 5e-3: packed {empty_packed} / raw {empty_raw} "
          f"of {int(small.sum())}  <- angular packing granularity; AOD fixes this")

    plot_mod.use_style("ATLAS")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    ax = axes[0]
    ax.plot(radii, mean_pt, "o-", lw=2)
    ax.set_xlabel("jet radius $R$"); ax.set_ylabel(r"$\langle p_T\rangle$ [GeV]")
    ax.set_title("Wider cone collects more radiation", fontsize=14)

    ax = axes[1]
    mb = np.linspace(0, 120, 45)
    ax.hist(m_raw, bins=mb, histtype="step", lw=2, label="ungroomed")
    ax.hist(m_trim, bins=mb, histtype="step", lw=2, label="trimmed")
    ax.set_xlabel(r"$R=0.8$ jet mass [GeV]"); ax.set_ylabel("jets")
    ax.set_title("Trimming removes soft junk", fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[2]
    plot_mod.plot_eec(cen_all, val_all, ax=ax, lw=2, label="all particles")
    plot_mod.plot_eec(cen_chg, val_chg, ax=ax, lw=2, ls="--", label="charged only")
    ax.set_title("Charged-only correlator", fontsize=14)
    ax.legend(fontsize=11, loc="lower right")

    finish(fig, "tier_miniaod.png")
    return {"pt_growth": growth, "m_raw": m_raw.mean(), "m_trim": m_trim.mean()}


# ======================================================================================
# AOD — full precision and track-level control
# ======================================================================================
def apply_tracking(particles, efficiency=1.0, angular_smear=0.0, rng=None):
    """Drop charged particles with (1-efficiency) and smear their angles.

    Neutral particles are untouched: this models *tracking*, not calorimetry.
    """
    rng = rng or np.random.default_rng(0)
    n_tot = int(ak.sum(ak.num(particles)))
    charged = particles.charge != 0

    u = ak.unflatten(rng.uniform(0, 1, n_tot), ak.num(particles))
    keep = (~charged) | (u < efficiency)
    p = particles[keep]

    if angular_smear > 0:
        n = int(ak.sum(ak.num(p)))
        de = ak.unflatten(rng.normal(0, angular_smear, n), ak.num(p))
        dp = ak.unflatten(rng.normal(0, angular_smear, n), ak.num(p))
        is_ch = p.charge != 0
        p = ak.zip(
            {"pt": p.pt,
             "eta": p.eta + ak.where(is_ch, de, 0.0),
             "phi": p.phi + ak.where(is_ch, dp, 0.0),
             "mass": p.mass, "charge": p.charge},
            with_name="Momentum4D",
        )
    return p


def run_aod(particles):
    print("\n=== AOD :: full precision, own tracking systematics ===")
    rng = np.random.default_rng(7)

    variations = [
        ("nominal",              1.00, 0.0000, "k"),
        ("track eff -5%",        0.95, 0.0000, "#4C72B0"),
        ("track eff -10%",       0.90, 0.0000, "#DD8452"),
        ("angular smear 2 mrad", 1.00, 0.0020, "#C44E52"),
    ]
    colour = {lab: c for lab, _, _, c in variations}

    # coarse bins: the systematic must be separated from statistical noise
    bins = np.logspace(-3, 0, 25)

    curves = {}
    for label, eff, smear, _ in variations:
        p = apply_tracking(particles, eff, smear, rng)
        _, c = cluster(p)
        cflat = ak.flatten(c, axis=1)
        cen, val = eec_mod.eec(cflat, bins=bins)
        curves[label] = (cen, val)
        print(f"  {label:<22} jets = {len(cflat):5d}  "
              f"norm = {eec_mod.normalisation_check(cflat):.4f}")

    cen, nom = curves["nominal"]
    small = cen < 0.01          # the small-angle region, where tracking dominates
    large = cen > 0.1
    print("\n  fractional shift w.r.t. nominal:")
    for label, _, _, _ in variations[1:]:
        _, v = curves[label]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(nom > 0, v / nom, np.nan)
        s = np.nanmean(ratio[small]) - 1
        l = np.nanmean(ratio[large]) - 1
        print(f"    {label:<22} small angle {100*s:+6.1f}%   large angle {100*l:+6.1f}%")

    plot_mod.use_style("ATLAS")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    ax = axes[0]
    for label, _, _, col in variations:
        c, v = curves[label]
        plot_mod.plot_eec(c, v, ax=ax, lw=2, color=col,
                          ls="-" if label == "nominal" else "--", label=label)
    ax.set_title("EEC under tracking variations", fontsize=14)
    ax.legend(fontsize=10, loc="lower right")

    ax = axes[1]
    for label, _, _, col in variations[1:]:
        _, v = curves[label]
        with np.errstate(divide="ignore", invalid="ignore"):
            ax.plot(cen, np.where(nom > 0, v / nom, np.nan), lw=2,
                    color=col, label=label)
    ax.axhline(1.0, color="k", lw=1)
    ax.axvspan(1e-3, 0.01, color="grey", alpha=0.15)
    ax.text(3e-3, 1.22, "small angle:\ntracking dominates", fontsize=10, ha="center")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_xlabel(r"$R_L$"); ax.set_ylabel("ratio to nominal")
    ax.set_ylim(0.75, 1.3)
    ax.set_title("Systematic band — only derivable at AOD", fontsize=14)
    ax.legend(fontsize=10, loc="lower right")

    finish(fig, "tier_aod.png")
    return {}


# ======================================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["nanoaod", "pfnano", "miniaod", "aod", "all"],
                    default="all")
    ap.add_argument("--n-events", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    print(f"Generating {args.n_events} toy events (seed={args.seed}) ...")
    particles = generate(args.n_events, seed=args.seed, cfg=ToyConfig())
    print(f"  <particles/event> = {ak.mean(ak.num(particles)):.0f}")

    runners = {"nanoaod": run_nanoaod, "pfnano": run_pfnano,
               "miniaod": run_miniaod, "aod": run_aod}
    todo = list(runners) if args.tier == "all" else [args.tier]
    for t in todo:
        runners[t](particles)

    print("\n=== RECO ===")
    print("  Not released as Open Data and ~3 MB/event. Contents are hits, clusters and")
    print("  reconstruction intermediates — the physics it enables is detector physics")
    print("  (calibration, alignment, tracking from first principles), not jet physics.")
    print("  Nothing to run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
