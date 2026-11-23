"""Detector check of the CMS 2011A EEC with the MOD SIM sample (Pythia 6 + CMS GEANT4 simulation).

Every SIM jet carries two particle lists:
  * ``pfcs`` -- detector-level PF candidates, reconstructed exactly like the data;
  * ``gens`` -- generator-level particles of the matched GEN jet (``gen_jet_pt = -1`` if none).
So the EEC can be computed at reco and at truth level *on the same jets*. Three comparisons:

  data / SIM-reco          ~ 1  -> the data features are reproduced by Pythia 6 + detector
  SIM-reco / SIM-gen             -> the detector effect itself (same jets, same generator)
  SIM-gen / Pythia 8 truth       -> Pythia 6 vs Pythia 8 at truth level

Weights: SIM is generated in 8 hard-parton pT-hat slices (SIM170 ... SIM1800). Each jet has a
cross-section weight; when only some files of a slice are used they are rescaled so that the
slice keeps its total cross section (``total_weights`` in energyflow's dataset info, the same
rule as ``energyflow.mod.load``). Slices are then added with these weights. n_eff is reported.

Statistical errors: jet-level bootstrap (Poisson(1) weights per jet, ``--n-boot`` replicas,
fixed ``--seed``). Reco and gen of the same jet use the same replica weights, so the correlation
between them is kept in the reco/gen ratio. The Pythia 8 baseline enters with its jet-level SEM
(Gaussian replicas).

Usage:
  python src/validate_eec_pythia.py --preset mod2011 --n-events 40000 --seed 12345
  python src/analyze_mod.py --n-files 3
  python src/analyze_mod_sim.py --files SIM170:1 SIM300:6 SIM470:6 SIM600:5 SIM800:3 \
      SIM1000:2 SIM1400:2 SIM1800:2 --data-files 3

Per-file per-jet EEC histograms are cached in data/interim/mod_sim/ and data/interim/mod_data/
(not committed); delete them to recompute.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from analyze_mod import (PT_EDGES, REGIONS, SELECTION, TURNOVER_MAX_RL,  # noqa: E402
                         _half_rise, _turnover, pfcs_to_constituents)
from jeteec import eec as eec_mod  # noqa: E402
from jeteec import plot as plot_mod  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "mod" / "datasets" / "CMS2011AJets"
CACHE = REPO / "data" / "interim"
BINS = np.logspace(-3, 0, 51)
SIM_SLICES = ("SIM170", "SIM300", "SIM470", "SIM600", "SIM800", "SIM1000", "SIM1400", "SIM1800")
# Neutral hadrons: PF neutral hadron (130) at reco; K0L (130) and neutrons (2112) at gen level.
# The MOD gen lists contain no K0S, Lambda or neutrinos (checked on SIM170 file 0).
NEUTRAL_HADRONS = (130, 2112)
# Finer large-angle ranges used in docs/mod_data_results.md for the second maximum.
WIDE = (("0.1-0.2", 0.1, 0.2), ("0.2-0.45", 0.2, 0.45))
VARIANTS = ("all", "ch", "noh0")


# --------------------------------------------------------------------------------------------
# small helpers (tested in tests/test_mod_sim.py)
# --------------------------------------------------------------------------------------------

def slice_weight_scale(total_weight: float, file_orig_weights) -> float:
    """Factor that restores a pT-hat slice's total cross section when only some files are used.

    ``total_weight`` is the slice total over all its files (energyflow dataset info);
    ``file_orig_weights`` are the summed jet weights of the files actually read, before any
    selection. Same rule as ``energyflow.mod.load``.
    """
    s = float(np.sum(file_orig_weights))
    if s <= 0:
        raise ValueError("files carry no weight")
    return float(total_weight) / s


def n_eff(w) -> float:
    """Effective number of jets (sum w)^2 / sum w^2."""
    w = np.asarray(w, dtype=float)
    return float(w.sum() ** 2 / (w**2).sum())


def region_mean_ratio(num, den, centres, lo, hi) -> float:
    """Mean over bins with lo < centre < hi of num/den (the convention of analyze_mod.py)."""
    num, den, centres = (np.asarray(a, dtype=float) for a in (num, den, centres))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = num / den
    m = (centres > lo) & (centres < hi) & np.isfinite(r)
    return float(np.mean(r[m]))


def region_mean(values, centres, lo, hi) -> float:
    m = (centres > lo) & (centres < hi)
    return float(np.mean(np.asarray(values)[m]))


def boot_means(hists, w, n_boot, rng, chunk=20000):
    """Bootstrap replicas of the weighted jet mean of several per-jet histogram arrays.

    ``hists``: list of (n_jets, n_bins) arrays for the *same* jets (e.g. reco and gen), which
    share the Poisson(1) replica weights. Returns a list of (n_boot, n_bins) arrays.
    """
    n = len(w)
    acc = [np.zeros((n_boot, h.shape[1])) for h in hists]
    sw = np.zeros(n_boot)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        p = rng.poisson(1.0, size=(n_boot, hi - lo)) * w[lo:hi]
        sw += p.sum(axis=1)
        for a, h in zip(acc, hists):
            a += p @ np.asarray(h[lo:hi], dtype=float)
    return [a / sw[:, None] for a in acc]


# --------------------------------------------------------------------------------------------
# per-file processing (cached)
# --------------------------------------------------------------------------------------------

def _hist(const):
    keep = ak.to_numpy(ak.num(const)) > 0
    out = np.zeros((len(const), len(BINS) - 1), dtype=np.float32)
    if keep.any():
        _, h = eec_mod.per_jet_histograms(const[keep], BINS, chunk_size=2000, dtype=np.float32)
        out[keep] = h
    return out


def _sets(parts):
    """All, charged-only and no-neutral-hadron constituents from MOD-format particle lists."""
    return {
        "all": pfcs_to_constituents(parts, chs=True),
        "ch": pfcs_to_constituents(parts, chs=True, charged_only=True),
        "noh0": pfcs_to_constituents(parts, chs=True, drop_pids=NEUTRAL_HADRONS),
    }


def _h0_fraction(parts):
    out = np.zeros(len(parts))
    for i, p in enumerate(parts):
        p = np.asarray(p)
        if len(p):
            out[i] = p[np.isin(np.abs(p[:, 4]), NEUTRAL_HADRONS), 0].sum() / p[:, 0].sum()
    return out


def process_file(path: str) -> str:
    """Per-jet EEC histograms of one MOD file -> cached .npz. Returns the cache path."""
    import energyflow as ef

    path = pathlib.Path(path)
    sim = path.name.startswith("SIM")
    out = CACHE / ("mod_sim" if sim else "mod_data") / (path.stem + ".npz")
    if out.exists():
        return str(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds = ef.mod.MODDataset(str(path), SELECTION, store_gens=sim)
    rec = {
        "corr_pt": np.asarray(ds.corr_jet_pts, float),
        "weight": np.asarray(ds.weights, float),
        "orig_total_weight": float(ds._orig_total_weight),
        "npv": np.asarray(ds.npvs, int),
        # ds.quality is the *column index* (5) of quality in jets_i; the values are ds.qualitys
        "quality": np.asarray(ds.qualitys, int),
        "h0_frac_reco": _h0_fraction(ds.pfcs),
    }
    sets = _sets(ds.pfcs)
    rec["n_reco"] = ak.to_numpy(ak.num(sets["all"]))
    for k, c in sets.items():
        rec[f"H_reco_{k}"] = _hist(c)
    if sim:
        gen_pt = np.asarray(ds.gen_jet_pts, float)
        n_gens = np.array([len(g) for g in ds.gens])
        rec["gen_pt"] = gen_pt
        rec["matched"] = (gen_pt > 0) & (n_gens > 0)
        rec["hard_pid"] = np.asarray(ds.hard_pids, int)
        rec["h0_frac_gen"] = _h0_fraction(ds.gens)
        gsets = _sets(ds.gens)  # vertex column is 0 at gen level, so CHS is a no-op
        rec["n_gen"] = ak.to_numpy(ak.num(gsets["all"]))
        for k, c in gsets.items():
            rec[f"H_gen_{k}"] = _hist(c)
    np.savez(out, **rec)
    return str(out)


def _sim_files(spec):
    files = []
    for item in spec:
        name, n = item.split(":")
        sub = RAW / f"{name}_Jet300_pT375-infGeV"
        for i in range(int(n)):
            f = sub / f"{name}_Jet300_pT375-infGeV_{i}_compressed.h5"
            if not f.exists():
                raise FileNotFoundError(f"{f} (download it first, see docs/mod_sim_results.md)")
            files.append((name, f))
    return files


def _concat(recs, keys):
    return {k: np.concatenate([r[k] for r in recs]) for k in keys}


def load_sim(spec, workers):
    import energyflow as ef

    info = ef.datasets.mod._get_dataset_info("CMS2011AJets")["total_weights"]
    files = _sim_files(spec)
    with ProcessPoolExecutor(workers) as ex:
        caches = list(ex.map(process_file, [str(f) for _, f in files]))
    by_slice: dict[str, list] = {}
    for (name, f), c in zip(files, caches):
        by_slice.setdefault(name, []).append((f, dict(np.load(c))))
    recs, table = [], []
    for name in SIM_SLICES:
        if name not in by_slice:
            continue
        items = by_slice[name]
        scale = slice_weight_scale(info[f"{name}_Jet300_pT375-infGeV"],
                                   [r["orig_total_weight"] for _, r in items])
        for _, r in items:
            r["w"] = r["weight"] * scale
            r["slice"] = np.full(len(r["w"]), SIM_SLICES.index(name))
            recs.append(r)
        w = np.concatenate([r["w"] for _, r in items])
        size = sum(f.stat().st_size for f, _ in items) / 1e6
        table.append((name, len(items), size, len(w), w.sum(), n_eff(w)))
    keys = [k for k in recs[0] if k != "orig_total_weight"]
    return _concat(recs, keys), table


def load_data(n_files, workers):
    sub = RAW / "CMS_Jet300_pT375-infGeV"
    files = [sub / f"CMS_Jet300_pT375-infGeV_{i}_compressed.h5" for i in range(n_files)]
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"{f}; run analyze_mod.py --n-files {n_files} first")
    with ProcessPoolExecutor(workers) as ex:
        caches = list(ex.map(process_file, [str(f) for f in files]))
    recs = [dict(np.load(c)) for c in caches]
    d = _concat(recs, [k for k in recs[0] if k != "orig_total_weight"])
    d["w"] = d["weight"]
    return d


# --------------------------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------------------------

class Sample:
    """Per-jet histograms + weights -> nominal shape and bootstrap replicas of the shape."""

    def __init__(self, label, hists: dict, w, n_boot, rng, pt=None):
        self.label, self.w, self.pt = label, w, pt
        dlog = np.diff(np.log(BINS))
        keys = list(hists)
        means = {k: (w @ np.asarray(hists[k], float)) / w.sum() for k in keys}
        reps = boot_means([hists[k] for k in keys], w, n_boot, rng)
        self.values = {k: means[k] / dlog for k in keys}
        self.shape = {k: eec_mod.shape(self.values[k], BINS) for k in keys}
        self.boot = {k: np.array([eec_mod.shape(r / dlog, BINS) for r in reps[i]])
                     for i, k in enumerate(keys)}
        self.boot_raw = {k: reps[i] / dlog for i, k in enumerate(keys)}
        self.n_jets, self.n_eff = len(w), n_eff(w)


class Fixed:
    """A sample known only by its mean and SEM (the Pythia 8 table): Gaussian replicas."""

    def __init__(self, label, values, errors, n_boot, rng):
        self.label = label
        self.values = {"all": values}
        self.shape = {"all": eec_mod.shape(values, BINS)}
        reps = values + errors * rng.standard_normal((n_boot, len(values)))
        self.boot = {"all": np.array([eec_mod.shape(r, BINS) for r in reps])}
        self.boot_raw = {"all": reps}


def ratio_with_err(a, ka, b, kb, lo, hi):
    c = np.sqrt(BINS[:-1] * BINS[1:])
    nom = region_mean_ratio(a.shape[ka], b.shape[kb], c, lo, hi)
    reps = [region_mean_ratio(x, y, c, lo, hi) for x, y in zip(a.boot[ka], b.boot[kb])]
    return nom, float(np.std(reps, ddof=1))


def wide_with_err(a, k):
    c = np.sqrt(BINS[:-1] * BINS[1:])
    out = []
    for _, lo, hi in WIDE:
        out.append((region_mean(a.shape[k], c, lo, hi),
                    float(np.std([region_mean(x, c, lo, hi) for x in a.boot[k]], ddof=1))))
    rise = out[1][0] / out[0][0]
    reps = [region_mean(x, c, *WIDE[1][1:]) / region_mean(x, c, *WIDE[0][1:]) for x in a.boot[k]]
    return out, (rise, float(np.std(reps, ddof=1)))


def half_rise_with_err(values, boot):
    c = np.sqrt(BINS[:-1] * BINS[1:])
    nom = _half_rise(c, values)
    return nom, float(np.std([_half_rise(c, b) for b in boot], ddof=1))


def fmt(x, e, nd=2):
    return f"{x:.{nd}f} ± {e:.{nd}f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", nargs="+", required=True,
                    help="SIM slices to read as NAME:N_FILES, e.g. SIM300:6")
    ap.add_argument("--data-files", type=int, default=3)
    ap.add_argument("--n-boot", type=int, default=100)
    ap.add_argument("--seed", type=int, default=2026, help="bootstrap seed")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--kfactors", action="store_true",
                    help="multiply SIM weights by energyflow.mod.kfactors('sim', pT, npv)")
    ap.add_argument("--pythia-table", default=str(REPO / "results" / "tables"
                                                  / "eec_pythia_mod2011.npz"))
    ap.add_argument("--tag", default="", help="suffix for output files")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    centres = np.sqrt(BINS[:-1] * BINS[1:])
    summary: dict = {"seed": args.seed, "n_boot": args.n_boot, "files": args.files,
                     "kfactors": args.kfactors}

    print("[1/5] SIM files -> per-jet EEC histograms (cached) ...")
    sim, table = load_sim(args.files, args.workers)
    print("      slice     files    MB     jets    sum w [nb]    n_eff   share of sum w")
    tot = sum(t[4] for t in table)
    for name, nf, mb, nj, sw, ne in table:
        print(f"      {name:8s} {nf:5d} {mb:7.0f} {nj:8d} {sw:12.4g} {ne:9.0f}   {sw / tot:7.1%}")
    summary["slices"] = [dict(zip(("name", "files", "MB", "jets", "sum_w", "n_eff"), t))
                         for t in table]
    w = sim["w"].copy()
    if args.kfactors:
        import energyflow as ef
        k = ef.mod.kfactors("sim", sim["corr_pt"], sim["npv"])
        print(f"      k-factors: mean {np.average(k, weights=w):.3f}, "
              f"range {k.min():.3f}-{k.max():.3f}")
        w = w * k
    m = sim["matched"]
    print(f"      total: {len(w)} jets, n_eff {n_eff(w):.0f}; matched to a GEN jet: "
          f"{np.average(m, weights=w):.4%}")
    print(f"      quality values in SIM: {dict(zip(*np.unique(sim['quality'], return_counts=True)))}")
    print(f"      <corr pT> {np.average(sim['corr_pt'], weights=w):.1f} GeV, <npv> "
          f"{np.average(sim['npv'], weights=w):.2f}, <n_reco> {np.average(sim['n_reco'], weights=w):.1f}, "
          f"<n_gen> {np.average(sim['n_gen'][m], weights=w[m]):.1f}")
    print(f"      neutral-hadron pT fraction: reco {np.average(sim['h0_frac_reco'], weights=w):.3f}, "
          f"gen {np.average(sim['h0_frac_gen'][m], weights=w[m]):.3f}")
    summary.update(n_jets=int(len(w)), n_eff=n_eff(w), matched=float(np.average(m, weights=w)),
                   quality=sorted(set(sim["quality"].tolist())))

    print("[2/5] data (cached per-jet histograms) ...")
    data = load_data(args.data_files, args.workers)
    print(f"      {len(data['w'])} jets, <corr pT> {np.average(data['corr_pt'], weights=data['w']):.1f}"
          f" GeV, <npv> {data['npv'].mean():.2f}, <n_reco> {data['n_reco'].mean():.1f}, "
          f"neutral-hadron pT fraction {data['h0_frac_reco'].mean():.3f}")

    print(f"[3/5] shapes and bootstrap ({args.n_boot} replicas, seed {args.seed}) ...")
    # SIM, matched jets: reco and gen of the same jets share bootstrap weights
    S = Sample("SIM reco", {f"reco_{k}": sim[f"H_reco_{k}"][m] for k in VARIANTS}
               | {f"gen_{k}": sim[f"H_gen_{k}"][m] for k in VARIANTS}, w[m], args.n_boot, rng)
    S_all = Sample("SIM reco (all jets)", {"all": sim["H_reco_all"]}, w, args.n_boot, rng)
    D = Sample("data", {k: data[f"H_reco_{k}"] for k in VARIANTS}, data["w"], args.n_boot, rng)
    P = None
    py = None
    if pathlib.Path(args.pythia_table).exists():
        py = np.load(args.pythia_table)
        P = Fixed("Pythia 8", py["values"], py["errors"], args.n_boot, rng)
    print(f"      SIM matched-jet n_eff {S.n_eff:.0f}; data n_eff {D.n_eff:.0f}")
    print("      SIM-reco all jets / matched jets, max |ratio-1| over bins: "
          f"{np.max(np.abs(S_all.shape['all'] / S.shape['reco_all'] - 1)):.4f}")

    print("[4/5] region ratios of unit-area shapes (± bootstrap std) ...")
    comps = [("data / SIM-reco", D, "all", S, "reco_all"),
             ("SIM-reco / SIM-gen", S, "reco_all", S, "gen_all"),
             ("SIM-gen / Pythia 8", S, "gen_all", P, "all"),
             ("SIM-reco / Pythia 8", S, "reco_all", P, "all"),
             ("data / Pythia 8", D, "all", P, "all"),
             ("data / SIM-reco, charged", D, "ch", S, "reco_ch"),
             ("SIM-reco / SIM-gen, charged", S, "reco_ch", S, "gen_ch"),
             ("data / SIM-reco, no h0", D, "noh0", S, "reco_noh0"),
             ("SIM-reco / SIM-gen, no h0", S, "reco_noh0", S, "gen_noh0"),
             ("SIM-reco: charged / all", S, "reco_ch", S, "reco_all"),
             ("data: charged / all", D, "ch", D, "all"),
             ("SIM-reco: no h0 / all", S, "reco_noh0", S, "reco_all"),
             ("data: no h0 / all", D, "noh0", D, "all"),
             ("SIM-gen: no h0 / all", S, "gen_noh0", S, "gen_all")]
    summary["regions"] = {}
    print("      " + " " * 30 + "".join(f"{n:>20s}" for n, _, _ in REGIONS))
    for name, a, ka, b, kb in comps:
        if a is None or b is None:
            continue
        cells = [ratio_with_err(a, ka, b, kb, lo, hi) for _, lo, hi in REGIONS]
        summary["regions"][name] = cells
        print(f"      {name:30s}" + "".join(f"{fmt(*c):>20s}" for c in cells))

    print("      large-angle shape means (second maximum): 0.1-0.2, 0.2-0.45, ratio")
    summary["wide"] = {}
    rows = [("data, all", D, "all"), ("data, no h0", D, "noh0"),
            ("SIM-reco, all", S, "reco_all"), ("SIM-reco, no h0", S, "reco_noh0"),
            ("SIM-gen, all", S, "gen_all"), ("SIM-gen, no h0", S, "gen_noh0")]
    if P is not None:
        rows.append(("Pythia 8", P, "all"))
    for name, a, k in rows:
        cells, rise = wide_with_err(a, k)
        summary["wide"][name] = cells + [rise]
        print(f"      {name:18s} " + "  ".join(fmt(*c, 3) for c in cells) + f"   {fmt(*rise, 3)}")

    print("[5/5] pT slices: half-rise R_L * <pT> [GeV] (± bootstrap std)")
    hr_rows = []
    lo_hi = list(zip(PT_EDGES[:-1], PT_EDGES[1:]))
    for lo, hi in lo_hi:
        row = {"lo": lo, "hi": hi}
        sd = (data["corr_pt"] >= lo) & (data["corr_pt"] < hi)
        Ds = Sample("", {"all": data["H_reco_all"][sd]}, data["w"][sd], args.n_boot, rng)
        mpt = np.average(data["corr_pt"][sd], weights=data["w"][sd])
        h, e = half_rise_with_err(Ds.values["all"], Ds.boot_raw["all"])
        row["data"] = (h * mpt, e * mpt, int(sd.sum()), float(sd.sum()))
        row["data_curve"] = Ds.values["all"]
        # SIM, binned in reco pT: reco and gen of the same jets
        ss = m & (sim["corr_pt"] >= lo) & (sim["corr_pt"] < hi)
        Ss = Sample("", {"reco": sim["H_reco_all"][ss], "gen": sim["H_gen_all"][ss]}, w[ss],
                    args.n_boot, rng)
        mpt_r = np.average(sim["corr_pt"][ss], weights=w[ss])
        mpt_g = np.average(sim["gen_pt"][ss], weights=w[ss])
        h, e = half_rise_with_err(Ss.values["reco"], Ss.boot_raw["reco"])
        row["reco"] = (h * mpt_r, e * mpt_r, int(ss.sum()), Ss.n_eff)
        h, e = half_rise_with_err(Ss.values["gen"], Ss.boot_raw["gen"])
        row["gen_recobin"] = (h * mpt_g, e * mpt_g, int(ss.sum()), Ss.n_eff)
        row["reco_curve"], row["gen_curve"] = Ss.values["reco"], Ss.values["gen"]
        row["mpt_reco"], row["mpt_gen"], row["mpt_data"] = mpt_r, mpt_g, mpt
        # SIM gen, binned in gen-jet pT (truth-level scaling)
        sg = m & (sim["gen_pt"] >= lo) & (sim["gen_pt"] < hi)
        Sg = Sample("", {"gen": sim["H_gen_all"][sg]}, w[sg], args.n_boot, rng)
        mpt_gg = np.average(sim["gen_pt"][sg], weights=w[sg])
        h, e = half_rise_with_err(Sg.values["gen"], Sg.boot_raw["gen"])
        row["gen_genbin"] = (h * mpt_gg, e * mpt_gg, int(sg.sum()), Sg.n_eff)
        if py is not None and "slice_values" in py.files:
            i = [k for k, (a, b) in enumerate(zip(py["slice_lo"], py["slice_hi"]))
                 if (a, b) == (lo, hi)]
            if i:
                i = i[0]
                v = py["slice_values"][i]
                mp = float(py["slice_mean_pt"][i])
                if "slice_errors" in py.files:
                    reps = v + py["slice_errors"][i] * rng.standard_normal((args.n_boot, len(v)))
                    h, e = half_rise_with_err(v, reps)
                else:
                    h, e = _half_rise(centres, v), float("nan")
                row["pythia8"] = (h * mp, e * mp, 0,
                                  float(py["slice_n_eff"][i]) if "slice_n_eff" in py.files else 0)
        hr_rows.append(row)
    print("      pT [GeV]    data             SIM-reco         SIM-gen (reco-pT bin)  "
          "SIM-gen (gen-pT bin)  Pythia 8        SIM n_eff")
    for r in hr_rows:
        p8 = fmt(*r["pythia8"][:2]) if "pythia8" in r else "-"
        print(f"      {r['lo']:4.0f}-{r['hi']:<5.0f} {fmt(*r['data'][:2]):16s} {fmt(*r['reco'][:2]):16s} "
              f"{fmt(*r['gen_recobin'][:2]):22s} {fmt(*r['gen_genbin'][:2]):21s} {p8:15s} "
              f"{r['reco'][3]:.0f}")
    summary["half_rise"] = [{k: v for k, v in r.items() if not k.endswith("curve")}
                            for r in hr_rows]
    for name, a, k in (("data", D, "all"), ("SIM-reco", S, "reco_all"),
                       ("SIM-gen", S, "gen_all")) + ((("Pythia 8", P, "all"),) if P else ()):
        print(f"      turnover R_L ({name}): {_turnover(centres, a.values[k]):.3f}, half-rise R_L "
              f"{_half_rise(centres, a.values[k]):.4f}")

    tables = REPO / "results" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    tag = args.tag
    with open(tables / f"mod_sim_summary{tag}.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    np.savez(tables / f"eec_mod2011_sim{tag}.npz", edges=BINS,
             **{f"sim_{k}": v for k, v in S.values.items()},
             **{f"data_{k}": v for k, v in D.values.items()})

    _plots(D, S, P, hr_rows, tag)
    return 0


def _plots(D, S, P, hr_rows, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_mod.use_style("CMS")
    outdir = REPO / "results" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    c = np.sqrt(BINS[:-1] * BINS[1:])

    def band(ax, a, k, num_den=None, **kw):
        if num_den is None:
            reps, nom = a.boot[k], a.shape[k]
        else:
            (b, kb) = num_den
            reps = a.boot[k] / b.boot[kb]
            nom = a.shape[k] / b.shape[kb]
        lo, hi = np.percentile(reps, [16, 84], axis=0)
        line, = ax.plot(c, nom, **kw)
        ax.fill_between(c, lo, hi, color=line.get_color(), alpha=0.25, lw=0)

    # 1) shapes and the three ratios
    fig, axs = plt.subplots(4, 1, figsize=(7, 10), sharex=True,
                            gridspec_kw={"height_ratios": [3, 1, 1, 1]})
    band(axs[0], D, "all", color="k", marker="o", ms=3, ls="", label="CMS 2011A data")
    band(axs[0], S, "reco_all", color="C0", lw=1.8, label="SIM reco (Pythia 6 + CMS sim.)")
    band(axs[0], S, "gen_all", color="C2", lw=1.8, ls="--", label="SIM gen (same jets)")
    if P is not None:
        band(axs[0], P, "all", color="C1", lw=1.5, ls=":", label="Pythia 8 truth")
    axs[0].set_yscale("log")
    axs[0].set_ylabel("EEC, unit area in range")
    plot_mod.fit_log_ylim(axs[0], D.shape["all"], S.shape["gen_all"])
    axs[0].legend(fontsize=10, loc="lower center")
    axs[0].set_title(r"$p_T > 375$ GeV, $|\eta| < 1.9$, anti-$k_T$ R = 0.5", fontsize=12)
    for ax, (a, k, b, kb, lab) in zip(axs[1:], [
            (D, "all", S, "reco_all", "data /\nSIM reco"),
            (S, "reco_all", S, "gen_all", "SIM reco /\nSIM gen"),
            (S, "gen_all", P, "all", "SIM gen /\nPythia 8")]):
        if b is None:
            continue
        band(ax, a, k, (b, kb), color="k", lw=1.2)
        ax.axhline(1, color="grey", ls=":")
        ax.set_ylim(0.2, 1.8)
        ax.set_ylabel(lab, fontsize=11)
    axs[-1].set_xscale("log")
    axs[-1].set_xlabel(r"$R_L$")
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / f"eec_mod2011_sim_shapes{tag}.png"))
    plt.close(fig)

    # 2) neutral hadrons at large angle
    fig, ax = plt.subplots(figsize=(7, 5))
    for a, k, lab, kw in ((D, "all", "data, all PFCs", dict(color="k", marker="o", ms=3, ls="")),
                          (D, "noh0", "data, no neutral hadrons", dict(color="k", ls="--")),
                          (S, "reco_all", "SIM reco, all", dict(color="C0")),
                          (S, "reco_noh0", "SIM reco, no neutral hadrons", dict(color="C0", ls="--")),
                          (S, "gen_all", "SIM gen, all", dict(color="C2")),
                          (S, "gen_noh0", "SIM gen, no K0L/n", dict(color="C2", ls="--"))):
        band(ax, a, k, label=lab, **kw)
    ax.set_xscale("log")
    ax.set_xlim(0.02, 1)
    ax.set_ylim(0.1, 0.35)
    ax.set_xlabel(r"$R_L$")
    ax.set_ylabel("EEC, unit area in range")
    ax.legend(fontsize=9, ncol=2, loc="lower center")
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / f"eec_mod2011_sim_neutral_hadrons{tag}.png"))
    plt.close(fig)

    # 3) pT slices vs R_L * pT, and the half-rise vs pT
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17, 5.2))
    for i, r in enumerate(hr_rows):
        col = f"C{i}"
        lab = f"{r['lo']:.0f}-{r['hi']:.0f} GeV"
        a1.plot(c * r["mpt_reco"], eec_mod.shape(r["reco_curve"], BINS), color=col, label=lab)
        a1.plot(c * r["mpt_data"], eec_mod.shape(r["data_curve"], BINS), "o", ms=2.5, color=col)
        a2.plot(c * r["mpt_gen"], eec_mod.shape(r["gen_curve"], BINS), color=col, label=lab)
    for ax, t in ((a1, "SIM reco (lines), data (points)"), (a2, "SIM gen, same jets")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.3, 300)
        ax.set_ylim(0.01, 0.5)
        ax.set_xlabel(r"$R_L \cdot \langle p_T \rangle$ [GeV]")
        ax.set_title(t, fontsize=12)
        ax.legend(fontsize=9, loc="lower right")
    a1.set_ylabel("EEC, unit area in range")
    x = np.array([(r["lo"] + r["hi"]) / 2 for r in hr_rows])
    for k, lab, kw in (("data", "data", dict(color="k", marker="o")),
                       ("reco", "SIM reco", dict(color="C0", marker="s")),
                       ("gen_recobin", "SIM gen (reco-pT bins)", dict(color="C2", marker="^")),
                       ("gen_genbin", "SIM gen (gen-pT bins)", dict(color="C2", marker="v", ls="--")),
                       ("pythia8", "Pythia 8 truth", dict(color="C1", marker="D", ls=":"))):
        pts = [(xx, r[k][0], r[k][1]) for xx, r in zip(x, hr_rows) if k in r]
        if pts:
            xx, y, e = map(np.array, zip(*pts))
            a3.errorbar(xx, y, yerr=e, label=lab, ms=5, capsize=2, **kw)
    a3.set_xlabel(r"jet $p_T$ slice centre [GeV]")
    a3.set_ylabel(r"half-rise $R_L \cdot \langle p_T \rangle$ [GeV]")
    a3.legend(fontsize=9)
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / f"eec_mod2011_sim_pt_slices{tag}.png"))
    plt.close(fig)
    print("      figures -> results/figures/eec_mod2011_sim_*.png")


if __name__ == "__main__":
    raise SystemExit(main())
