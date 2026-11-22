"""One real analysis per data tier.

Each tier gets the analysis that its information content actually permits, run end to end
on a common toy sample and producing a figure plus printed numbers. Where the toy
generator knows the answer, it is printed next to the measurement (a closure test), and
every headline number carries its statistical uncertainty.

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
import math
import pathlib
import sys

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import vector  # noqa: E402

from jeteec import eec as eec_mod  # noqa: E402
from jeteec import jets as jets_mod  # noqa: E402
from jeteec import plot as plot_mod  # noqa: E402
from jeteec import stats  # noqa: E402
from jeteec.kinematics import delta_r, sum_p4, wrap_phi  # noqa: E402
from toy_generator import ToyConfig, generate, ue_density  # noqa: E402

vector.register_awkward()

REPO = pathlib.Path(__file__).resolve().parents[1]
FIGDIR = REPO / "results" / "figures"
R_NOMINAL = 0.4
JET_PT_MIN = 100.0
COUNT_PT_MIN = 60.0        # softer threshold for jet counting (R32), as in real analyses
MAX_ABS_ETA = 2.5
PT_FIT_MAX = 1000.0
CHI_MAX = 16.0
TOY_NOTE = "toy sample — pipeline demo, not a measurement"

# MiniAOD packing as implemented in CMSSW (DataFormats/PatCandidates/src/
# PackedCandidate.cc): eta and phi of packed PF candidates are stored as int16 over
# +-6.0 and +-3.2, pT and mass as 16-bit floats. Every PF candidate is kept.
ETA_LSB = 6.0 / 32767
PHI_LSB = 3.2 / 32767


# ======================================================================================
# shared helpers
# ======================================================================================
def cluster(particles, R=R_NOMINAL, min_pt=JET_PT_MIN, max_abs_eta=MAX_ABS_ETA):
    """anti-kT jets, pT-ordered, plus constituents that keep every particle field.

    Goes through :func:`jeteec.jets.cluster`, which guards against the fastjet layout
    trap (see the README) — so no field stripping and no dR re-matching is needed here.
    """
    cfg = jets_mod.JetConfig(algorithm="antikt", R=R, min_pt=min_pt,
                             max_abs_eta=max_abs_eta)
    return jets_mod.cluster(particles, cfg)


def annotate(fig):
    fig.text(0.995, 0.005, TOY_NOTE, ha="right", va="bottom",
             fontsize=9, color="grey", style="italic")


def finish(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    annotate(fig)
    out = FIGDIR / name
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    try:
        shown = out.relative_to(REPO)
    except ValueError:
        shown = out
    print(f"  figure -> {shown}")


# ======================================================================================
# NanoAOD — jet four-vectors only
# ======================================================================================
def run_nanoaod(particles, truth):
    print("\n=== NanoAOD :: jet kinematics only ===")
    # cluster once with a loose threshold, then apply per-analysis cuts — as in a real
    # analysis, where the multiplicity ratio uses a softer jet definition than the spectrum
    all_jets, _ = cluster(particles, min_pt=COUNT_PT_MIN)
    jets = all_jets[all_jets.pt > JET_PT_MIN]
    # the generator's own jet axes, with the same acceptance: what perfect
    # reconstruction would give
    t_acc = truth[abs(truth.eta) < MAX_ABS_ETA]

    # -- (a) inclusive jet pT spectrum: unbinned maximum-likelihood power law -----------
    pt = ak.to_numpy(ak.flatten(jets.pt))
    n_fit, n_err = stats.power_law_index(pt, JET_PT_MIN, PT_FIT_MAX)
    t_pt = ak.to_numpy(ak.flatten(t_acc.pt))
    n_true, n_true_err = stats.power_law_index(t_pt, JET_PT_MIN, PT_FIT_MAX)
    print(f"  jets                    : {len(pt)}")
    print(f"  spectrum index n        : {n_fit:.2f} ± {n_err:.2f}   "
          f"(dN/dpT ~ pT^-n, unbinned ML fit)")
    print(f"    generator jet axes    : {n_true:.2f} ± {n_true_err:.2f}   "
          f"(hard-scatter input 5.0)")

    # -- (b) dijet angular distribution -------------------------------------------------
    # jets are pT-ordered, so [:, 0] and [:, 1] are the leading pair; chi is defined
    # with rapidity, not pseudorapidity
    two = jets[ak.num(jets) >= 2]
    dy = abs(ak.to_numpy(two[:, 0].rapidity - two[:, 1].rapidity))
    chi = np.exp(dy)
    chi = chi[chi < CHI_MAX]
    chi_counts, chi_edges = np.histogram(chi, bins=np.linspace(1, CHI_MAX, 16))
    flat = chi_counts.sum() / len(chi_counts)
    chi2 = float(((chi_counts - flat) ** 2 / flat).sum())
    ndf = len(chi_counts) - 1
    print(f"  dijet pairs             : {len(chi)}")
    print(f"  flat-chi test           : chi2/ndf = {chi2:.1f}/{ndf}  "
          f"(flat = Rutherford-like t-channel)")

    # -- (c) jet multiplicity ratio -----------------------------------------------------
    njet = ak.num(all_jets)   # softer 60 GeV threshold, standard for counting ratios
    n2 = int(ak.sum(njet >= 2))
    n3 = int(ak.sum(njet >= 3))
    r32, r32_err = stats.binomial_ratio(n3, n2)
    t_n = ak.num(t_acc[t_acc.pt > COUNT_PT_MIN])
    r32_true, r32_true_err = stats.binomial_ratio(int(ak.sum(t_n >= 3)),
                                                  int(ak.sum(t_n >= 2)))
    print(f"  N(>=2 jets)             : {n2}")
    print(f"  N(>=3 jets)             : {n3}")
    print(f"  R32 = N3/N2             : {r32:.3f} ± {r32_err:.3f}   "
          f"(alpha_s-sensitive in real data)")
    print(f"    generator jet axes    : {r32_true:.3f} ± {r32_true_err:.3f}")

    # -- figure -------------------------------------------------------------------------
    plot_mod.use_style("ATLAS")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)

    ax = axes[0]
    bins = np.logspace(np.log10(JET_PT_MIN), np.log10(PT_FIT_MAX), 26)
    counts, edges = np.histogram(pt, bins=bins)
    centres = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    filled = counts > 0
    ax.errorbar(centres[filled], (counts / widths)[filled],
                yerr=(np.sqrt(counts) / widths)[filled], fmt="o", ms=4)
    in_fit = (pt >= JET_PT_MIN) & (pt < PT_FIT_MAX)
    norm = in_fit.sum() / stats.power_law_integral(n_fit, JET_PT_MIN, PT_FIT_MAX)
    xs = np.logspace(np.log10(JET_PT_MIN), np.log10(pt.max()), 100)
    ax.plot(xs, norm * xs ** -n_fit, "-", color="crimson",
            label=f"ML fit: $n = {n_fit:.2f} \\pm {n_err:.2f}$")
    plot_mod.log_pt_axis(ax); ax.set_yscale("log")
    ax.set_xlabel(r"jet $p_T$ [GeV]"); ax.set_ylabel(r"$dN/dp_T$ [GeV$^{-1}$]")
    ax.set_title("Inclusive jet spectrum", fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[1]
    chi_norm = chi_counts / chi_counts.sum()
    chi_cen = 0.5 * (chi_edges[:-1] + chi_edges[1:])
    ax.errorbar(chi_cen, chi_norm, yerr=np.sqrt(chi_counts) / chi_counts.sum(),
                xerr=0.5 * np.diff(chi_edges), fmt="o", ms=4)
    ax.axhline(1 / len(chi_counts), color="crimson", ls="--", lw=1.5,
               label="flat (Rutherford)")
    ax.set_xlabel(r"$\chi = e^{|y_1-y_2|}$"); ax.set_ylabel("normalised")
    ax.set_xlim(1, CHI_MAX)
    ax.set_ylim(0, chi_norm.max() * 1.6)
    ax.set_title(f"Dijet angular distribution ($\\chi^2$/ndf = {chi2:.0f}/{ndf})",
                 fontsize=14)
    ax.legend(fontsize=12)

    ax = axes[2]
    mult = np.bincount(ak.to_numpy(njet), minlength=5)[:6]
    ax.bar(range(len(mult)), mult, color="#4C72B0", width=0.7)
    ax.set_xticks(range(len(mult)))
    ax.set_xlim(-0.6, len(mult) - 0.4)
    ax.set_xlabel("jets per event"); ax.set_ylabel("events")
    ax.set_title(f"Jet multiplicity  ($R_{{32}}$ = {r32:.3f} ± {r32_err:.3f})",
                 fontsize=14)

    finish(fig, "tier_nanoaod.png")
    return {"n_jets": len(pt), "n_fit": n_fit, "n_err": n_err, "n_true": n_true,
            "chi2": chi2, "ndf": ndf, "r32": r32, "r32_err": r32_err,
            "r32_true": r32_true}


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
    dr = delta_r(cflat, jflat)
    zc_ = cflat.pt / ak.sum(cflat.pt, axis=1)
    girth = ak.to_numpy(ak.sum(zc_ * dr, axis=1))
    ptd = ak.to_numpy(np.sqrt(ak.sum(zc_**2, axis=1)))
    print(f"  <girth>                 : {girth.mean():.4f}")
    print(f"  <pTD>                   : {ptd.mean():.4f}")

    # -- (d) EEC -------------------------------------------------------------------------
    norm = eec_mod.normalisation_check(cflat)
    cen, val, err = eec_mod.eec(cflat, return_errors=True)
    peak = int(np.argmax(val))
    print(f"  EEC normalisation       : {norm:.6f}  (must be 1.0)")
    print(f"  EEC stat. error at peak : {100 * err[peak] / val[peak]:.1f}%  "
          f"(R_L = {cen[peak]:.3f}, jet-to-jet spread)")

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
    plot_mod.log_pt_axis(ax)
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
    ax.fill_between(cen, np.clip(val - err, 1e-12, None), val + err,
                    alpha=0.3, lw=0, label="stat. uncertainty")
    ax.axvline(R_NOMINAL, color="grey", ls=":", lw=1.5)
    ax.set_title(f"EEC (norm = {norm:.4f})", fontsize=14)
    ax.legend(fontsize=11, loc="lower center")

    finish(fig, "tier_pfnano.png")
    return {"n_jets": len(jflat), "n_const": float(ak.mean(nconst)), "norm": norm,
            "girth": float(girth.mean()), "ptd": float(ptd.mean())}


# ======================================================================================
# MiniAOD — every candidate in the event (compressed)
# ======================================================================================
def _float16(x):
    return ak.values_astype(ak.values_astype(x, np.float16), np.float64)


def compress(particles, eta_lsb=ETA_LSB, phi_lsb=PHI_LSB):
    """Mimic MiniAOD packing of PF candidates.

    Every candidate is kept; eta and phi go onto fixed grids (int16 in CMSSW), pT and
    mass to 16-bit floats. The defaults follow CMSSW. Pass coarser ``eta_lsb``/``phi_lsb``
    to see what a lossier format would do.
    """
    return ak.zip(
        {
            "pt": _float16(particles.pt),
            "eta": np.round(particles.eta / eta_lsb) * eta_lsb,
            "phi": np.round(wrap_phi(particles.phi) / phi_lsb) * phi_lsb,
            "mass": _float16(particles.mass),
            "charge": particles.charge,
        },
        with_name="Momentum4D",
    )


def trim(const, r_sub=0.2, f_cut=0.05):
    """Jet trimming: recluster constituents into subjets, drop the soft ones."""
    sub = jets_mod.subjets(const, algorithm="kt", R=r_sub)
    jet_pt = ak.sum(const.pt, axis=1)
    keep = sub.pt > f_cut * jet_pt
    return sum_p4(sub[keep]), sum_p4(sub)


def cos_weighted_area(R):
    """integral of cos(d_phi) over a disc of radius R in (eta, phi): 2 pi R J1(R).

    A particle at d_phi from the jet axis adds pT * cos(d_phi) to the jet pT, so this —
    not pi R^2 — is the area that converts an underlying-event pT density into jet pT.
    """
    j1 = sum((-1) ** m * (R / 2) ** (2 * m + 1) / (math.factorial(m) * math.factorial(m + 1))
             for m in range(25))
    return 2 * np.pi * R * j1


def matched_gain(particles, r_small, r_large, max_abs_eta):
    """Per-jet pT gain when the same jet is reclustered with a larger radius.

    Each selected R=r_small jet is matched to the nearest r_large jet (dR < r_small/2).
    Matching removes the bias of comparing mean pT of two different jet populations
    (a larger R also lifts new jets over the pT threshold).
    """
    small, _ = cluster(particles, R=r_small)
    large, _ = cluster(particles, R=r_large, min_pt=JET_PT_MIN / 5)
    pairs = ak.cartesian({"s": small, "l": large}, axis=1, nested=True)
    d = delta_r(pairs.s, pairs.l)
    best = ak.argmin(d, axis=2, keepdims=True)
    d_best = ak.to_numpy(ak.fill_none(ak.flatten(ak.firsts(d[best], axis=2)), np.inf))
    gain = ak.flatten(ak.firsts((pairs.l.pt - pairs.s.pt)[best], axis=2))
    gain = ak.to_numpy(ak.fill_none(gain, np.nan))
    eta = ak.to_numpy(ak.flatten(small.eta))
    ok = (d_best < r_small / 2) & (abs(eta) < max_abs_eta)
    return gain[ok]


def run_miniaod(particles, cfg):
    print("\n=== MiniAOD :: all packed candidates ===")
    packed = compress(particles)
    print(f"  candidates kept         : {ak.sum(ak.num(packed))} / "
          f"{ak.sum(ak.num(particles))}  (MiniAOD keeps every PF candidate)")

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
    # |eta| < 2: the R = 1.0 cone then lies inside the generated underlying event
    gain = matched_gain(packed, R_NOMINAL, 1.0, max_abs_eta=cfg.ue_eta_max - 1.0)
    gain_err = gain.std(ddof=1) / np.sqrt(len(gain))
    ue_expected = ue_density(cfg) * (cos_weighted_area(1.0) - cos_weighted_area(R_NOMINAL))
    merged = float(np.mean(gain > 20.0))
    print(f"  per-jet gain R=0.4->1.0 : {gain.mean():+.2f} ± {gain_err:.2f} GeV  "
          f"(matched jets, |eta| < {cfg.ue_eta_max - 1.0:.0f})")
    print(f"    generated UE predicts : {ue_expected:+.2f} GeV   "
          f"median gain {np.median(gain):+.2f} GeV")
    print(f"    jets gaining > 20 GeV : {100 * merged:.1f}%  "
          f"(absorb a neighbouring jet at R = 1.0)")

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
    # fastjet carries `charge` through to the constituents, so no dR re-matching
    _, c4 = cluster(packed, R=R_NOMINAL)
    call = ak.flatten(c4, axis=1)
    call = call[ak.num(call) >= 2]
    cchg = call[call.charge != 0]
    cchg = cchg[ak.num(cchg) >= 2]
    cen_all, val_all = eec_mod.eec(call)
    cen_chg, val_chg = eec_mod.eec(cchg)
    print(f"  charged fraction of pT  : "
          f"{float(ak.sum(cchg.pt)/ak.sum(call.pt)):.3f}")

    # -- (d) does packing hurt the EEC? ---------------------------------------------------
    fine = np.logspace(-3, 0, 51)
    small = np.sqrt(fine[:-1] * fine[1:]) < 5e-3

    def small_angle_eec(sample):
        _, c = cluster(sample, R=R_NOMINAL)
        c = ak.flatten(c, axis=1)
        return eec_mod.eec(c[ak.num(c) >= 2], bins=fine)[1]

    v_raw = small_angle_eec(particles)
    v_packed = small_angle_eec(packed)
    v_coarse = small_angle_eec(compress(particles, eta_lsb=1e-3, phi_lsb=1e-3))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(v_packed / v_raw - 1)
    dev_large = float(np.nanmax(rel[~small]))
    dev_small = float(np.nanmax(rel[small]))
    empty = [int((v[small] == 0).sum()) for v in (v_raw, v_packed, v_coarse)]
    print(f"  EEC, packed vs unpacked : <= {100 * dev_large:.1f}% per bin above R_L = 5e-3, "
          f"up to {100 * dev_small:.0f}% below")
    print(f"    (the int16 eta step, {ETA_LSB:.1e}, is as wide as a bin at R_L ~ 1e-3: "
          f"pairs hop between bins)")
    print(f"  empty EEC bins < 5e-3   : unpacked {empty[0]}, MiniAOD {empty[1]}, "
          f"hypothetical 1e-3 lattice {empty[2]}  (of {int(small.sum())})")

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
    return {"gain": float(gain.mean()), "gain_err": float(gain_err),
            "ue_expected": float(ue_expected), "m_raw": float(m_raw.mean()),
            "m_trim": float(m_trim.mean()), "packing_dev_large": dev_large,
            "packing_dev_small": dev_small, "empty_bins": empty}


# ======================================================================================
# AOD — full precision and track-level control
# ======================================================================================
def apply_tracking(particles, efficiency=1.0, angular_smear=0.0, u=None, rng=None):
    """Drop charged particles with probability (1-efficiency) and smear their angles.

    Neutral particles are untouched: this models *tracking*, not calorimetry. ``u`` is
    one uniform number per particle; passing the same ``u`` to every variation makes
    them nested (the -10% variation loses the -5% tracks plus more), so variations differ
    by the efficiency change rather than by independent random noise.
    """
    rng = rng or np.random.default_rng(0)
    charged = particles.charge != 0
    if u is None:
        n_tot = int(ak.sum(ak.num(particles)))
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


def _bootstrap_ratio(nominal, varied, n_events, rng, n_boot=500):
    """Replicas of (varied EEC / nominal EEC) per bin, from a Poisson bootstrap over events.

    ``nominal`` and ``varied`` are (per-jet histograms, event index of each jet). Both are
    re-weighted with the *same* event weights, so the strong correlation between a
    variation and the nominal (same events) is kept — which is what makes a 1% shift
    measurable at all.
    """
    w = rng.poisson(1.0, size=(n_boot, n_events)).astype(float)

    def mean(hist, evt):
        wj = w[:, evt]
        return (wj @ hist) / wj.sum(axis=1, keepdims=True)

    with np.errstate(divide="ignore", invalid="ignore"):
        return mean(*varied) / mean(*nominal)


def run_aod(particles):
    print("\n=== AOD :: full precision, own tracking systematics ===")
    n_events = len(particles)
    rng = np.random.default_rng(7)
    u = ak.unflatten(rng.uniform(0, 1, int(ak.sum(ak.num(particles)))),
                     ak.num(particles))            # common random numbers

    variations = [
        ("nominal",              1.00, 0.0000, "k"),
        ("track eff -5%",        0.95, 0.0000, "#4C72B0"),
        ("track eff -10%",       0.90, 0.0000, "#DD8452"),
        ("angular smear 2 mrad", 1.00, 0.0020, "#C44E52"),
    ]

    # coarse bins: the systematic must be separated from statistical noise
    bins = np.logspace(-3, 0, 25)

    curves, per_jet = {}, {}
    for label, eff, smear, _ in variations:
        p = apply_tracking(particles, eff, smear, u=u, rng=np.random.default_rng(8))
        jets, c = cluster(p)
        cflat = ak.flatten(c, axis=1)
        edges, hist = eec_mod.eec_per_jet(cflat, bins=bins)
        evt = np.repeat(np.arange(n_events), ak.to_numpy(ak.num(jets)))
        per_jet[label] = (hist, evt)
        cen, val = eec_mod.eec(cflat, bins=bins)
        curves[label] = (cen, val)
        print(f"  {label:<22} jets = {len(cflat):5d}  "
              f"norm = {eec_mod.normalisation_check(cflat):.4f}")

    cen, nom = curves["nominal"]
    small = (cen < 0.01) & (nom > 0)   # the small-angle region, where tracking dominates
    large = cen > 0.1
    print("\n  fractional shift w.r.t. nominal (± stat., bootstrap over events):")
    shifts, bands = {}, {}
    boot_rng = np.random.default_rng(9)
    for label, _, _, _ in variations[1:]:
        _, v = curves[label]
        reps = _bootstrap_ratio(per_jet["nominal"], per_jet[label], n_events, boot_rng)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(nom > 0, v / nom, np.nan)
        s = np.nanmean(ratio[small]) - 1
        l = np.nanmean(ratio[large]) - 1
        s_err = np.nanstd(np.nanmean(reps[:, small], axis=1))
        l_err = np.nanstd(np.nanmean(reps[:, large], axis=1))
        shifts[label] = (s, s_err, l, l_err)
        bands[label] = (ratio, np.nanstd(reps, axis=0))
        print(f"    {label:<22} small angle {100*s:+5.1f} ± {100*s_err:3.1f}%   "
              f"large angle {100*l:+5.1f} ± {100*l_err:3.1f}%")
    worst = max(abs(s) / se for s, se, _, _ in shifts.values())
    print(f"  largest small-angle shift = {worst:.1f} sigma: no variation is resolved "
          f"with {n_events} events")

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
        ratio, band = bands[label]
        ax.plot(cen, ratio, lw=2, color=col, label=label)
        ax.fill_between(cen, ratio - band, ratio + band, color=col, alpha=0.18, lw=0)
    ax.axhline(1.0, color="k", lw=1)
    ax.axvspan(1e-3, 0.01, color="grey", alpha=0.15)
    ax.text(3e-3, 1.22, "small angle:\nwhere tracking matters", fontsize=10, ha="center")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1)
    ax.set_xlabel(r"$R_L$"); ax.set_ylabel("ratio to nominal (band: stat.)")
    ax.set_ylim(0.75, 1.3)
    ax.set_title("Systematic band — only derivable at AOD", fontsize=14)
    ax.legend(fontsize=10, loc="lower right")

    finish(fig, "tier_aod.png")
    return {"shifts": shifts,
            "n_jets": {label: len(per_jet[label][1]) for label in per_jet}}


# ======================================================================================
def main(argv=None) -> int:
    global FIGDIR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["nanoaod", "pfnano", "miniaod", "aod", "all"],
                    default="all")
    ap.add_argument("--n-events", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--figdir", type=pathlib.Path, default=FIGDIR,
                    help="where to write the figures (default: results/figures)")
    args = ap.parse_args(argv)
    FIGDIR = args.figdir

    cfg = ToyConfig()
    print(f"Generating {args.n_events} toy events (seed={args.seed}) ...")
    particles, truth = generate(args.n_events, seed=args.seed, cfg=cfg,
                                return_truth=True)
    print(f"  <particles/event> = {ak.mean(ak.num(particles)):.0f}")

    runners = {"nanoaod": lambda: run_nanoaod(particles, truth),
               "pfnano": lambda: run_pfnano(particles),
               "miniaod": lambda: run_miniaod(particles, cfg),
               "aod": lambda: run_aod(particles)}
    todo = list(runners) if args.tier == "all" else [args.tier]
    for t in todo:
        runners[t]()

    print("\n=== RECO ===")
    print("  Not released as Open Data and ~3 MB/event. Contents are hits, clusters and")
    print("  reconstruction intermediates — the physics it enables is detector physics")
    print("  (calibration, alignment, tracking from first principles), not jet physics.")
    print("  Nothing to run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
