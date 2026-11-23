"""First real-data EEC: CMS 2011A Jet primary dataset in MOD HDF5 format.

Data: CMS 2011A Open Data, Jet primary dataset, reprocessed into MOD HDF5 by Komiske, Mastandrea,
Metodiev, Naik, Thaler (arXiv:1908.08542). Zenodo 10.5281/zenodo.3340205, CC BY 4.0. 18 files of
~111 MB, read with ``energyflow``. Jets are anti-kT R = 0.5, the two hardest per event, with
PF candidates (PFCs).

Selection (as in arXiv:1908.08542):
  * JEC-corrected jet pT > 375 GeV (HLT_Jet300 fully efficient), |eta_jet| < 1.9
  * jet quality >= 2 ("medium"), to reject fake jets
  * PFCs: charged hadron subtraction (CHS) -- drop charged PFCs from pileup vertices

Checks performed:
  1. EEC normalisation == 1 on real constituents (same unit test as on Pythia)
  2. EEC in jet-pT slices, against R_L and R_L * pT
  3. comparison with the matched Pythia 8 baseline (``--preset mod2011`` of
     validate_eec_pythia.py), if its table is present

Usage:
  python src/analyze_mod.py --n-files 1                      # downloads file 0 (~111 MB)
  python src/analyze_mod.py --files path/to/CMS_..._0_compressed.h5

**This is detector-level data compared with generator-level MC.** No unfolding, no systematic
uncertainties. Shape differences at small angle are expected from detector resolution and from
the PFC reconstruction thresholds. It is a first look, not a measurement.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from jeteec import eec as eec_mod  # noqa: E402
from jeteec import plot as plot_mod  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "mod"
PT_EDGES = (375.0, 450.0, 550.0, 700.0, 1000.0)
# The data EEC is nearly flat from R_L ~ 0.02 to the jet edge and has a second, large-angle
# maximum. The global maximum then lands at R_L ~ 0.3 and says nothing about the
# hadronisation turnover. The turnover is therefore searched for below this R_L.
TURNOVER_MAX_RL = 0.1
REGIONS = (("R_L < 0.01", 0.0, 0.01), ("0.01 < R_L < 0.1", 0.01, 0.1),
           ("0.1 < R_L < 0.5", 0.1, 0.5))
SELECTION = "corr_jet_pts > 375 & abs_jet_eta < 1.9 & quality >= 2"


def pfcs_to_constituents(pfcs, chs: bool = True, charged_only: bool = False,
                         pt_min: float = 0.0, drop_pids=()) -> ak.Array:
    """Convert MOD per-jet PFC arrays into a jagged Momentum4D array, one sublist per jet.

    MOD stores (pt, y, phi, m, pid, vertex) with *rapidity* y. The EEC here uses
    pseudorapidity, like the Pythia baseline, so eta is computed from (pt, y, m). For
    massive candidates y != eta, and mixing them would shift small-angle separations.

    ``drop_pids`` removes candidates by |PDG id|, e.g. ``(130,)`` for PF neutral hadrons.
    It is a diagnostic: it shows which PFC type carries a feature of the EEC.
    """
    import vector

    vector.register_awkward()
    counts, rows = [], []
    for p in pfcs:
        p = np.asarray(p, dtype=float)
        keep = np.ones(len(p), dtype=bool)
        charged = _is_charged(p[:, 4])
        if chs:
            keep &= ~(charged & (p[:, 5] > 0))  # pileup vertex: id >= 1
        if charged_only:
            keep &= charged
        if pt_min > 0:
            keep &= p[:, 0] >= pt_min
        if len(drop_pids):
            keep &= ~np.isin(np.abs(p[:, 4]), drop_pids)
        p = p[keep]
        counts.append(len(p))
        rows.append(p)
    flat = np.concatenate(rows) if rows else np.zeros((0, 6))
    pt, y, phi, m = flat[:, 0], flat[:, 1], flat[:, 2], flat[:, 3]
    mt = np.sqrt(pt**2 + m**2)
    eta = np.arcsinh(mt * np.sinh(y) / pt)
    phi = (phi + np.pi) % (2 * np.pi) - np.pi
    cnt = np.asarray(counts)
    return ak.zip(
        {"pt": ak.unflatten(pt, cnt), "eta": ak.unflatten(eta, cnt),
         "phi": ak.unflatten(phi, cnt), "mass": ak.unflatten(m, cnt)},
        with_name="Momentum4D",
    )


def _is_charged(pids):
    """Charge mask from PDG ids, via energyflow (the same rule the MOD authors use)."""
    from energyflow.utils import ischrgd

    return np.asarray(ischrgd(pids), dtype=bool)


def _turnover(centres, values, rl_max=TURNOVER_MAX_RL):
    """R_L of the EEC maximum below ``rl_max`` (the hadronisation turnover)."""
    m = centres < rl_max
    return float(centres[m][np.argmax(values[m])])


def _half_rise(centres, values, rl_max=TURNOVER_MAX_RL):
    """R_L where the small-angle side first reaches half the turnover value (log-interpolated).

    A sharper scaling probe than the turnover itself when the maximum is broad.
    """
    m = centres < rl_max
    c, v = centres[m], values[m]
    i_max = int(np.argmax(v))
    half = 0.5 * v[i_max]
    i = int(np.argmax(v[: i_max + 1] >= half))
    if i == 0:
        return float(c[0])
    x0, x1 = np.log(c[i - 1]), np.log(c[i])
    return float(np.exp(x0 + (half - v[i - 1]) / (v[i] - v[i - 1]) * (x1 - x0)))


def load_data(files, n_files):
    import energyflow as ef

    if files:
        dsets = [ef.mod.MODDataset(str(f), SELECTION, store_gens=False) for f in files]
        return ef.mod.MODDataset(datasets=dsets) if len(dsets) > 1 else dsets[0]
    RAW.mkdir(parents=True, exist_ok=True)
    return ef.mod.load(SELECTION, amount=n_files, cache_dir=str(RAW), collection="CMS2011AJets",
                       dataset="cms", validate_files=True, store_gens=False, verbose=1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", nargs="*", help="local MOD .h5 files (skip download)")
    ap.add_argument("--n-files", type=int, default=1, help="files to download (18 total)")
    ap.add_argument("--chunk-size", type=int, default=2000)
    ap.add_argument("--pythia-table", default=str(REPO / "results" / "tables"
                                                  / "eec_pythia_mod2011.npz"))
    args = ap.parse_args()

    print("[1/4] loading MOD data ...")
    ds = load_data(args.files, args.n_files)
    print(f"      {len(ds)} jets after selection: {SELECTION}")
    pt = np.asarray(ds.corr_jet_pts, dtype=float)
    w = np.asarray(ds.weights, dtype=float)
    print(f"      <corrected pT> = {np.average(pt, weights=w):.1f} GeV, "
          f"weights: min {w.min():.3g} max {w.max():.3g} (n_eff {w.sum()**2/(w**2).sum():.0f})")
    print(f"      npv: mean {np.mean(ds.npvs):.1f}  (pileup; CHS removes charged PU)")

    print("[2/4] building constituents (CHS) ...")
    const = pfcs_to_constituents(ds.pfcs, chs=True)
    const_ch = pfcs_to_constituents(ds.pfcs, chs=True, charged_only=True)
    const_noh0 = pfcs_to_constituents(ds.pfcs, chs=True, drop_pids=(130,))
    print(f"      <n_const> all = {np.average(ak.to_numpy(ak.num(const)), weights=w):.1f}, "
          f"charged = {np.average(ak.to_numpy(ak.num(const_ch)), weights=w):.1f}")

    print("[3/4] EEC ...")
    norm = eec_mod.normalisation_check(const, chunk_size=args.chunk_size)
    print(f"      normalisation = {norm:.6f} (expect 1.0)")
    ok = abs(norm - 1.0) < 1e-6

    res = eec_mod.eec_with_errors(const, chunk_size=args.chunk_size, jet_weights=w)
    res_ch = eec_mod.eec_with_errors(const_ch[ak.num(const_ch) > 0],
                                     chunk_size=args.chunk_size,
                                     jet_weights=w[ak.to_numpy(ak.num(const_ch)) > 0])
    res_noh0 = eec_mod.eec_with_errors(const_noh0[ak.num(const_noh0) > 0],
                                       chunk_size=args.chunk_size,
                                       jet_weights=w[ak.to_numpy(ak.num(const_noh0)) > 0])
    slices = eec_mod.eec_in_pt_bins(const, pt, PT_EDGES, chunk_size=args.chunk_size,
                                    jet_weights=w)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_mod.use_style("CMS")
    outdir = REPO / "results" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    plotted = []
    print("      pT slice [GeV]   jets   peak R_L   turnover R_L (x<pT> GeV)"
          "   half-rise R_L (x<pT> GeV)   median rel. err")
    for lo, hi, r in slices:
        if r is None:
            print(f"      {lo:5.0f}-{hi:<5.0f}        0   (empty)")
            continue
        sel = (pt >= lo) & (pt < hi)
        mean_pt = float(np.average(pt[sel], weights=w[sel]))
        good = r.values > 0
        pk = r.centres[np.argmax(r.values)]
        to = _turnover(r.centres, r.values)
        hr = _half_rise(r.centres, r.values)
        rel = float(np.median(r.errors[good] / r.values[good]))
        print(f"      {lo:5.0f}-{hi:<5.0f} {r.n_jets:8d}   {pk:8.3f}   {to:8.3f} ({to * mean_pt:5.1f})"
              f"          {hr:8.4f} ({hr * mean_pt:5.2f})          {rel:6.1%}")
        v, e = r.values[good], r.errors[good]
        lo_edge, _ = plot_mod.log_band(v, e)
        plotted.append(v)
        for ax, x in ((ax1, r.centres), (ax2, r.centres * mean_pt)):
            ax.errorbar(x[good], v, yerr=(v - lo_edge, e), fmt="o-", ms=3, lw=1,
                        label=f"{lo:.0f}-{hi:.0f} GeV")
    for ax, xl in ((ax1, r"$R_L$"), (ax2, r"$R_L \cdot \langle p_T \rangle$ [GeV]")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
        ax.legend(title="CMS 2011A Open Data, jet $p_T$", fontsize=10, title_fontsize=11,
                  loc="lower center")
    ax1.set_ylabel(r"$\frac{1}{\sigma}\frac{d\sigma}{d\log R_L}$")
    plot_mod.fit_log_ylim(ax1, *plotted)
    fig.tight_layout()
    out = outdir / "eec_mod2011_data_pt_slices.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")

    tables = REPO / "results" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    np.savez(tables / "eec_mod2011_data.npz", edges=res.edges, values=res.values,
             errors=res.errors, values_charged=res_ch.values, errors_charged=res_ch.errors,
             values_no_h0=res_noh0.values, errors_no_h0=res_noh0.errors,
             n_jets=res.n_jets)

    print("[4/4] comparison with the Pythia 8 baseline ...")
    dlog = np.diff(np.log(res.edges))

    def shape(v):
        return v / np.sum(v * dlog)

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.5, 7.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    d_shape = shape(res.values)
    k = d_shape[res.values > 0] / res.values[res.values > 0]
    k = float(np.median(k))
    lo_edge, hi_edge = plot_mod.log_band(res.values, res.errors)
    ax.plot(res.centres, d_shape, "ko", ms=3, label="CMS 2011A data (all PFCs, CHS)")
    ax.fill_between(res.centres, lo_edge * k, hi_edge * k, color="k", alpha=0.2)
    ax.plot(res.centres, shape(res_ch.values), "o", ms=3, color="grey",
            label="CMS 2011A data (charged PFCs)")
    ax.plot(res.centres, shape(res_noh0.values), "--", lw=1.2, color="C0",
            label="CMS 2011A data (no neutral-hadron PFCs)")
    print("      data shape vs all PFCs:          " + "   ".join(n for n, _, _ in REGIONS))
    for name, v in (("charged PFCs only", res_ch.values),
                    ("without neutral hadrons", res_noh0.values)):
        ratio = shape(v) / d_shape
        cells = [np.mean(ratio[(res.centres > a) & (res.centres < b) & np.isfinite(ratio)])
                 for _, a, b in REGIONS]
        print(f"      {name:24s}         " + "   ".join(f"{c:10.2f}" for c in cells))
    curves = [d_shape]
    py_path = pathlib.Path(args.pythia_table)
    if py_path.exists():
        py = np.load(py_path)
        if not np.allclose(py["edges"], res.edges):
            print("      WARNING: Pythia table has different binning; skipping comparison")
        else:
            p_shape = shape(py["values"])
            curves.append(p_shape)
            ax.plot(res.centres, p_shape, "-", lw=2, color="C1",
                    label="Pythia 8 truth, matched selection")
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = d_shape / p_shape
            axr.plot(res.centres, ratio, "ko", ms=3)
            for name, a, b in REGIONS:
                good = (res.centres > a) & (res.centres < b) & np.isfinite(ratio)
                print(f"      data / Pythia shape, {name:17s}: {np.mean(ratio[good]):.2f}")
            pk_d = res.centres[np.argmax(res.values)]
            pk_p = res.centres[np.argmax(py["values"])]
            print(f"      EEC peak R_L: data {pk_d:.3f}, Pythia {pk_p:.3f}")
            print(f"      half-rise R_L: data {_half_rise(res.centres, res.values):.4f}, "
                  f"Pythia {_half_rise(res.centres, py['values']):.4f}")
            if "slice_values" in py.files:
                print("      half-rise R_L*<pT> [GeV] by slice:   data   Pythia")
                d_by = {(lo, hi): (r, pt[(pt >= lo) & (pt < hi)]) for lo, hi, r in slices}
                for lo, hi, mpt, v in zip(py["slice_lo"], py["slice_hi"],
                                          py["slice_mean_pt"], py["slice_values"]):
                    r, ptsel = d_by.get((float(lo), float(hi)), (None, None))
                    if r is None:
                        continue
                    hd = _half_rise(r.centres, r.values) * float(np.mean(ptsel))
                    hp = _half_rise(res.centres, v) * float(mpt)
                    print(f"        {lo:5.0f}-{hi:<5.0f}                       {hd:5.2f}   {hp:5.2f}")
            print(f"      turnover R_L (R_L < {TURNOVER_MAX_RL}): "
                  f"data {_turnover(res.centres, res.values):.3f}, "
                  f"Pythia {_turnover(res.centres, py['values']):.3f}")
    else:
        print(f"      no Pythia table at {py_path}; run "
              "validate_eec_pythia.py --preset mod2011 first")
    axr.axhline(1, color="grey", ls=":")
    axr.set_ylim(0, 2)
    axr.set_xscale("log")
    axr.set_xlabel(r"$R_L$")
    axr.set_ylabel("data / Pythia")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("EEC, unit area in range")
    plot_mod.fit_log_ylim(ax, *curves)
    ax.legend(fontsize=10)
    ax.set_title(r"$p_T > 375$ GeV, $|\eta| < 1.9$, anti-$k_T$ R = 0.5", fontsize=12)
    out = outdir / "eec_mod2011_data_vs_pythia.png"
    plot_mod.save(fig, str(out))
    plt.close(fig)
    print(f"      figure -> {out.relative_to(REPO)}")

    print("\nRESULT:", "PASS" if ok else "FAIL (normalisation off)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
