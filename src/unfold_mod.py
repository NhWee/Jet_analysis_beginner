"""Unfolding of the CMS 2011A EEC (MOD) to particle level with the MOD SIM response.

Main observable: the EEC of **charged particles with pT > 1 GeV** inside anti-kT R = 0.5 jets
(pT > 375 GeV, |eta| < 1.9, medium quality), normalised per jet with the scalar pT sum of
those particles (``jeteec.eec`` convention) and shown as a unit-area shape per jet-pT bin.
The same at reco level uses charged PFCs after CHS. The pT > 1 GeV threshold follows
arXiv:2201.07800 (same data, "to minimize acceptance effects") and arXiv:2402.13864. All PFCs
with pT > 1 GeV are unfolded the same way for comparison.

Method (details and choices: docs/mod_unfolding_results.md):
  * particle matching reco <-> gen (same charge, dR < 0.05, closest pT; arXiv:2402.13864);
  * pair-level response in (jet pT bin x R_L bin), reco jet pT = JEC-corrected pT, gen jet pT
    = matched GEN jet pT, so jet-pT migration is part of the response;
  * fakes (reco pairs with an unmatched particle) are subtracted as a fraction, misses enter
    through the efficiency; D'Agostini iterations (``jeteec.unfold``), Pythia 6 SIM-gen prior.
Each jet goes to one of ``N_GROUPS`` random sub-samples, so the response can be split
(closure) and resampled (its statistical uncertainty).

Usage (inputs as in docs/mod_sim_results.md; per-file results cached in data/interim/unfold/):
  for s in 12345 12346 12347 12348; do python src/validate_eec_pythia.py --preset mod2011 \
      --n-events 50000 --seed $s --charged-only --const-pt-min 1 --tag mod2011_ch1_s$s; done
  python src/unfold_mod.py --files SIM170:1 SIM300:6 SIM470:6 SIM600:5 SIM800:3 SIM1000:2 \
      SIM1400:2 SIM1800:2 --data-files 3
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import awkward as ak
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from analyze_mod import PT_EDGES, SELECTION, _half_rise  # noqa: E402
from analyze_mod_sim import RAW, SIM_SLICES, _sim_files, n_eff, slice_weight_scale  # noqa: E402
from jeteec import eec as eec_mod  # noqa: E402
from jeteec import matching  # noqa: E402
from jeteec import plot as plot_mod  # noqa: E402
from jeteec import unfold  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
CACHE = REPO / "data" / "interim" / "unfold"
BINS = np.logspace(-3, 0, 51)          # same R_L binning as analyze_mod_sim.py
NR = len(BINS) - 1
# reco jet bins: the four analysis bins + overflow (>= 1000 GeV)
RECO_JET_EDGES = np.array(list(PT_EDGES) + [np.inf])
# gen jet bins: underflow (< 375 GeV, jets that pass at reco only) + four bins + overflow
GEN_JET_EDGES = np.array([0.0] + list(PT_EDGES) + [np.inf])
N_RECO_JET, N_GEN_JET = len(RECO_JET_EDGES) - 1, len(GEN_JET_EDGES) - 1
N_GROUPS = 20
PT_MIN = 1.0                            # GeV, particle threshold at both levels
TRACK_DROP = 0.03                       # tracking-efficiency variation (see docs)
JER_WIDEN = 1.10                        # jet-pT resolution variation (see docs)
VARIANTS = ("ch1", "ch1_trk", "ch1_jer", "all1")
MAX_DR = 0.05


def jet_bin(pt, edges):
    """Index of pt in ``edges`` (-1 outside), lower edge inclusive."""
    b = np.searchsorted(edges, pt, side="right") - 1
    b[(b < 0) | (b >= len(edges) - 1)] = -1
    return b


def mod_to_particles(parts, chs=True, charged_only=False, pt_min=0.0, drop_frac=0.0, rng=None):
    """MOD (pt, y, phi, m, pid, vertex) per jet -> awkward (pt, eta, phi, charge) per jet.

    Rapidity is converted to pseudorapidity as in ``analyze_mod.pfcs_to_constituents``.
    CHS removes charged candidates from pileup vertices (vertex >= 1). ``drop_frac`` removes
    that fraction of the charged particles at random (tracking-efficiency variation).
    """
    from energyflow.utils import pids2chrgs

    counts = np.array([len(p) for p in parts])
    flat = np.concatenate([np.asarray(p, float) for p in parts]) if len(parts) else np.zeros((0, 6))
    q = np.asarray(pids2chrgs(flat[:, 4]), float) if len(flat) else np.zeros(0)
    keep = np.ones(len(flat), bool)
    charged = q != 0
    if chs:
        keep &= ~(charged & (flat[:, 5] > 0))
    if charged_only:
        keep &= charged
    keep &= flat[:, 0] >= pt_min
    if drop_frac > 0:
        keep &= ~(charged & (rng.random(len(flat)) < drop_frac))
    jet = np.repeat(np.arange(len(counts)), counts)
    counts = np.bincount(jet[keep], minlength=len(counts))
    f = flat[keep]
    pt, y, phi, m = f[:, 0], f[:, 1], f[:, 2], f[:, 3]
    mt = np.sqrt(pt**2 + m**2)
    eta = np.arcsinh(mt * np.sinh(y) / pt) if len(pt) else pt
    phi = (phi + np.pi) % (2 * np.pi) - np.pi
    return ak.zip({k: ak.unflatten(v, counts) for k, v in
                   (("pt", pt), ("eta", eta), ("phi", phi), ("charge", q[keep]))})


def _seed(path, salt):
    return int(hashlib.sha256(f"{pathlib.Path(path).name}:{salt}".encode()).hexdigest()[:12], 16)


def process_file(path: str, chunk: int = 5000) -> str:
    """Pair-level response (SIM) or per-jet reco histograms (data) of one MOD file, cached."""
    import energyflow as ef

    path = pathlib.Path(path)
    sim = path.name.startswith("SIM")
    out = CACHE / (path.stem + ".npz")
    if out.exists():
        return str(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds = ef.mod.MODDataset(str(path), SELECTION, store_gens=sim)
    rec: dict = {"corr_pt": np.asarray(ds.corr_jet_pts, float),
                 "weight": np.asarray(ds.weights, float),
                 "orig_total_weight": float(ds._orig_total_weight)}
    pfcs = list(ds.pfcs)
    if not sim:
        for v, kw in (("ch1", dict(charged_only=True)), ("all1", {})):
            hs = []
            for lo in range(0, len(pfcs), chunk):
                p = mod_to_particles(pfcs[lo:lo + chunk], pt_min=PT_MIN, **kw)
                hs.append(matching.per_jet_hist(p, BINS).astype(np.float32))
            rec[f"H_{v}"] = np.concatenate(hs)
        np.savez(out, **rec)
        return str(out)

    gens = list(ds.gens)
    gen_pt = np.asarray(ds.gen_jet_pts, float)
    matched = (gen_pt > 0) & np.array([len(g) > 0 for g in gens])
    grp = np.random.default_rng(_seed(path, "groups")).integers(0, N_GROUPS, len(gen_pt))
    drop_rng = np.random.default_rng(_seed(path, "trackdrop"))
    rjb = jet_bin(rec["corr_pt"], RECO_JET_EDGES)
    gjb = jet_bin(gen_pt, GEN_JET_EDGES)
    # JER variation: reco pT pulled away from gen pT by JER_WIDEN (resolution 10% worse)
    rjb_jer = jet_bin(gen_pt + JER_WIDEN * (rec["corr_pt"] - gen_pt), RECO_JET_EDGES)
    rec.update(gen_pt=gen_pt, matched=matched, group=grp)
    idx = np.nonzero(matched)[0]
    acc: dict = {}
    for lo in range(0, len(idx), chunk):
        sel = idx[lo:lo + chunk]
        pf = [pfcs[i] for i in sel]
        gn = [gens[i] for i in sel]
        for v in VARIANTS:
            ch = v.startswith("ch")
            g = mod_to_particles(gn, chs=False, charged_only=ch, pt_min=PT_MIN)
            r = mod_to_particles(pf, charged_only=ch, pt_min=PT_MIN,
                                 drop_frac=TRACK_DROP if v == "ch1_trk" else 0.0, rng=drop_rng)
            rb_ = rjb_jer[sel] if v == "ch1_jer" else rjb[sel]
            res = matching.pair_response(g, r, gjb[sel], rb_, rec["weight"][sel], grp[sel],
                                         N_GROUPS, BINS, N_GEN_JET, N_RECO_JET, MAX_DR)
            # gen pairs split by the jet's reco bin (for the jet-migration study)
            gj, _, _, gR, gw, _ = matching._pairs(g)
            b = matching.rl_bin(gR, BINS)
            ok = (b >= 0) & (rb_[gj] >= 0)
            res["T_by_reco"] = np.bincount(
                (rb_[gj[ok]] * N_GEN_JET + gjb[sel][gj[ok]]) * NR + b[ok],
                weights=(gw * rec["weight"][sel][gj])[ok],
                minlength=N_RECO_JET * N_GEN_JET * NR).reshape(N_RECO_JET, N_GEN_JET * NR)
            # matched-particle fractions (pT-weighted), a sanity number for the docs
            g2r, rm = matching.match_particles(g, r, MAX_DR)
            gpt = ak.to_numpy(ak.flatten(g.pt)).astype(float)
            rpt = ak.to_numpy(ak.flatten(r.pt)).astype(float)
            res["frac"] = np.array([np.sum(g2r >= 0), len(g2r), gpt[g2r >= 0].sum(), gpt.sum(),
                                    rm.sum(), len(rm), rpt[rm].sum(), rpt.sum()], float)
            for k, a in res.items():
                key = f"{v}_{k}"
                acc[key] = acc[key] + a if key in acc else a
    rec.update(acc)
    np.savez(out, **rec)
    return str(out)


def _data_files(n):
    sub = RAW / "CMS_Jet300_pT375-infGeV"
    files = [sub / f"CMS_Jet300_pT375-infGeV_{i}_compressed.h5" for i in range(n)]
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f)
    return files


def load(spec, n_data, workers):
    """Process (cached) and combine SIM files with slice weights; read data per-jet histograms."""
    import energyflow as ef

    info = ef.datasets.mod._get_dataset_info("CMS2011AJets")["total_weights"]
    sim_files = _sim_files(spec)
    data_files = _data_files(n_data)
    with ProcessPoolExecutor(workers) as ex:
        caches = list(ex.map(process_file, [str(f) for _, f in sim_files]
                             + [str(f) for f in data_files]))
    sim_c, data_c = caches[:len(sim_files)], caches[len(sim_files):]
    by_slice: dict = {}
    for (name, _), c in zip(sim_files, sim_c):
        by_slice.setdefault(name, []).append(dict(np.load(c)))
    S: dict = {}
    jets = {"corr_pt": [], "gen_pt": [], "w": [], "matched": [], "group": []}
    for name in SIM_SLICES:
        if name not in by_slice:
            continue
        items = by_slice[name]
        scale = slice_weight_scale(info[f"{name}_Jet300_pT375-infGeV"],
                                   [r["orig_total_weight"] for r in items])
        for r in items:
            for k, a in r.items():
                if k.startswith(VARIANTS):
                    a = a if k.endswith("_frac") else a * scale
                    S[k] = S[k] + a if k in S else a
            for k in ("corr_pt", "gen_pt", "matched", "group"):
                jets[k].append(r[k])
            jets["w"].append(r["weight"] * scale)
    S["jets"] = {k: np.concatenate(v) for k, v in jets.items()}
    recs = [dict(np.load(c)) for c in data_c]
    D = {k: np.concatenate([r[k] for r in recs]) for k in ("corr_pt", "weight", "H_ch1", "H_all1")}
    return S, D



# --------------------------------------------------------------------------------------------
# analysis helpers
# --------------------------------------------------------------------------------------------

DLOG = np.diff(np.log(BINS))
CENTRES = np.sqrt(BINS[:-1] * BINS[1:])
GEN_REPORT = range(1, 5)       # gen jet bins reported: 375-450, 450-550, 550-700, 700-1000
REGIONS = (("R_L < 0.01", 0.0, 0.01), ("0.01-0.1", 0.01, 0.1), ("0.1-0.5", 0.1, 0.5))


def gen_shapes(t):
    """Unit-area shapes (4, NR) of the reported gen jet bins from a flat gen histogram."""
    t = np.asarray(t, float).reshape(N_GEN_JET, NR)
    return np.array([eec_mod.shape(t[k] / DLOG, BINS) for k in GEN_REPORT])


def reco_shapes(r):
    """Unit-area shapes (4, NR) of the four analysis reco jet bins (overflow not shown)."""
    r = np.asarray(r, float).reshape(N_RECO_JET, NR)
    return np.array([eec_mod.shape(r[k] / DLOG, BINS) for k in range(4)])


def region_means(shape_ratio):
    """Mean over bins of a ratio (4, NR) in the three R_L regions -> (4, 3)."""
    out = np.zeros((shape_ratio.shape[0], len(REGIONS)))
    for j, (_, lo, hi) in enumerate(REGIONS):
        m = (CENTRES > lo) & (CENTRES < hi)
        out[:, j] = np.nanmean(shape_ratio[:, m], axis=1)
    return out


def half_rises(shapes, mean_pt):
    """Half-rise R_L * <pT> [GeV] per pT bin."""
    return np.array([_half_rise(CENTRES, s) * m for s, m in zip(shapes, mean_pt)])


class Response:
    """Response of one variant from summed (or resampled) group histograms."""

    def __init__(self, S, v, groups=None):
        g = np.arange(N_GROUPS) if groups is None else np.asarray(groups)
        self.T = S[f"{v}_T"][g].sum(0)
        self.M = S[f"{v}_M"][g].sum(0)
        self.R = S[f"{v}_R"][g].sum(0)
        self.A = unfold.response_from_counts(self.M, self.T)
        matched_reco = self.M.sum(1)
        self.fake = np.zeros_like(self.R)
        ok = self.R > 0
        self.fake[ok] = 1 - matched_reco[ok] / self.R[ok]

    def unfold(self, d, n_iter, prior=None, return_all=False):
        return unfold.dagostini(d, self.A, self.T if prior is None else prior, n_iter,
                                fake_fraction=self.fake, return_all=return_all)


def data_hist(D, v, w=None, pt_scale=1.0):
    """Flat reco histogram of the data (reco jet bin x R_L), optionally with jet weights."""
    rjb = jet_bin(D["corr_pt"] * pt_scale, RECO_JET_EDGES)
    w = D["weight"] if w is None else w
    out = np.zeros((N_RECO_JET, NR))
    for k in range(N_RECO_JET):
        s = rjb == k
        out[k] = np.asarray(w[s], float) @ np.asarray(D[f"H_{v}"][s], float)
    return out.ravel()


def data_boot(D, v, n_boot, rng, chunk=20000):
    """Poisson(1) jet-level bootstrap replicas of the data reco histogram, (n_boot, n_reco)."""
    rjb = jet_bin(D["corr_pt"], RECO_JET_EDGES)
    H = D[f"H_{v}"]
    out = np.zeros((n_boot, N_RECO_JET, NR))
    for lo in range(0, len(rjb), chunk):
        hi = min(lo + chunk, len(rjb))
        p = rng.poisson(1.0, size=(n_boot, hi - lo)) * D["weight"][lo:hi]
        for k in range(N_RECO_JET):
            s = rjb[lo:hi] == k
            if s.any():
                out[:, k] += p[:, s] @ np.asarray(H[lo:hi][s], float)
    return out.reshape(n_boot, -1)


def mean_gen_pt(S):
    j = S["jets"]
    gb = jet_bin(j["gen_pt"], GEN_JET_EDGES)
    m = j["matched"]
    return np.array([np.average(j["gen_pt"][m & (gb == k)], weights=j["w"][m & (gb == k)])
                     for k in GEN_REPORT])


def mean_reco_pt(D):
    rb = jet_bin(D["corr_pt"], RECO_JET_EDGES)
    return np.array([np.average(D["corr_pt"][rb == k], weights=D["weight"][rb == k])
                     for k in range(4)])


def pythia8_shapes(paths):
    """Pythia 8 per-slice EEC, errors (SEM), <pT> and n_eff, averaged over independent runs.

    Several equal-size runs (different seeds) are averaged with equal weight; the error of the
    mean is sqrt(sum sigma^2) / n_runs. (One 200k-event run does not fit in memory.)
    """
    vs, es, pts, ns = [], [], [], []
    for path in paths:
        py = np.load(path)
        if not np.allclose(py["edges"], BINS):
            raise ValueError("Pythia 8 table has a different R_L binning")
        lo = list(py["slice_lo"])
        idx = [lo.index(e) for e in PT_EDGES[:-1]]
        vs.append(py["slice_values"][idx])
        es.append(py["slice_errors"][idx])
        pts.append(py["slice_mean_pt"][idx])
        ns.append(py["slice_n_eff"][idx])
    n = len(vs)
    return (np.mean(vs, axis=0), np.sqrt(np.sum(np.square(es), axis=0)) / n,
            np.mean(pts, axis=0), np.sum(ns, axis=0))


def choose_iterations(t_all, boot_all, n_max):
    """First n >= 2 where p(chi2 between t_n and t_{n-1}) > 0.05, sigma from the bootstrap."""
    rows = []
    sel = np.concatenate([np.arange(k * NR, (k + 1) * NR) for k in GEN_REPORT])
    for n in range(2, n_max + 1):
        sig = np.std(boot_all[:, n - 1, sel], axis=0, ddof=1)
        p, chi2 = unfold.iteration_pvalue(t_all[n - 1][sel], t_all[n - 2][sel], sig)
        rows.append((n, chi2, p))
    chosen = next((n for n, _, p in rows if p > 0.05), n_max)
    return chosen, rows


def fmt(x, e, nd=3):
    return f"{x:.{nd}f} ± {e:.{nd}f}"


# --------------------------------------------------------------------------------------------
# main analysis
# --------------------------------------------------------------------------------------------

def run_variant(S, D, v, args, rng, py8, mean_pt_gen, summary, label):
    """Unfold one particle definition with validation and systematics. Returns a result dict."""
    out: dict = {"label": label}
    nom = Response(S, v)
    d = data_hist(D, v)
    n_max = args.n_iter_max
    print(f"  [{v}] fake fraction (reco-weighted, all reco bins): "
          f"{1 - nom.M.sum() / nom.R.sum():.3f}; gen weight in matched pairs "
          f"{S[f'{v}_Mg'].sum() / nom.T.sum():.3f}; <eff_j> "
          f"{np.average(nom.A.sum(0), weights=nom.T):.3f}")

    # data: nominal + bootstrap at every iteration
    t_nom, t_all = nom.unfold(d, n_max, return_all=True)
    db = data_boot(D, v, args.n_boot, rng)
    boot_all = np.array([nom.unfold(x, n_max, return_all=True)[1] for x in db])  # (B, n, G)
    n_it = args.n_iter or choose_iterations(t_all, boot_all, n_max)[0]
    _, rows = choose_iterations(t_all, boot_all, n_max)
    print(f"  [{v}] iterations: chosen n = {n_it} (n: chi2/ndf={4 * NR}, p) " +
          ", ".join(f"{n}: {c:.1f}, {p:.2f}" for n, c, p in rows[:8]))
    out["iter_rows"] = rows
    out["n_iter"] = n_it
    t = t_all[n_it - 1]
    shp = gen_shapes(t)
    boot_shp = np.array([gen_shapes(b[n_it - 1]) for b in boot_all])
    stat = np.std(boot_shp, axis=0, ddof=1)
    hr = half_rises(shp, mean_pt_gen)
    hr_boot = np.array([half_rises(b, mean_pt_gen) for b in boot_shp])
    out.update(shape=shp, stat=stat, hr=hr, hr_stat=np.std(hr_boot, axis=0, ddof=1))

    # detector level and bin-by-bin
    reco_shp = reco_shapes(d)
    out["reco_shape"] = reco_shp
    sim_gen = gen_shapes(nom.T)
    sim_reco = reco_shapes(nom.R)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.where(sim_reco > 0, sim_gen / sim_reco, 0.0)
    bbb = np.array([eec_mod.shape(x, BINS) for x in reco_shp * C])
    out.update(bbb=bbb, C=C, sim_gen=sim_gen, sim_reco=sim_reco,
               hr_bbb=half_rises(bbb, mean_pt_gen), hr_sim_gen=half_rises(sim_gen, mean_pt_gen))

    # ---- closure: response from even groups, pseudo-data and truth from odd groups ----
    even, odd = np.arange(0, N_GROUPS, 2), np.arange(1, N_GROUPS, 2)
    r_even, r_odd = Response(S, v, even), Response(S, v, odd)
    clo = {}
    for n in (1, 2, 3, 4, 6, 8, n_it):
        u = gen_shapes(r_even.unfold(r_odd.R, n))
        with np.errstate(divide="ignore", invalid="ignore"):
            clo[n] = u / gen_shapes(r_odd.T)
    # stat error of the closure ratio: resample the odd groups (truth and reco together)
    cb = []
    for _ in range(args.n_boot_resp):
        g = rng.choice(odd, size=len(odd))
        ro = Response(S, v, g)
        with np.errstate(divide="ignore", invalid="ignore"):
            cb.append(gen_shapes(r_even.unfold(ro.R, n_it)) / gen_shapes(ro.T))
    clo_err = np.std(cb, axis=0, ddof=1)
    out.update(closure=clo[n_it], closure_err=clo_err, closure_by_n=clo)
    out["closure_regions"] = region_means(clo[n_it] - 1)
    pulls = (clo[n_it] - 1) / np.where(clo_err > 0, clo_err, np.inf)
    out["closure_chi2"] = (float(np.sum(pulls**2)), int(np.sum(clo_err > 0)))
    hr_clo = half_rises(gen_shapes(r_even.unfold(r_odd.R, n_it)), mean_pt_gen)
    out["closure_hr"] = (hr_clo, half_rises(gen_shapes(r_odd.T), mean_pt_gen))
    print(f"  [{v}] closure (n={n_it}), mean (unfolded/truth - 1) per region and pT bin:")
    for k, (lo, hi) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:])):
        print(f"        {lo:.0f}-{hi:.0f}: " + "  ".join(f"{x:+.4f}" for x in out["closure_regions"][k])
              + f"   half-rise {hr_clo[k]:.3f} vs truth {out['closure_hr'][1][k]:.3f} GeV")
    print(f"        chi2/ndf of closure pulls: {out['closure_chi2'][0]:.1f}/{out['closure_chi2'][1]}")

    # ---- stress test: truth reweighted to the Pythia 8 shape ----
    syst: dict = {}
    if py8 is not None:
        p8v = py8[0]
        T = nom.T.reshape(N_GEN_JET, NR).copy()
        r = np.ones((N_GEN_JET, NR))
        for i, k in enumerate(GEN_REPORT):
            a, b = eec_mod.shape(p8v[i], BINS), eec_mod.shape(T[k] / DLOG, BINS)
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = np.where((a > 0) & (b > 0), a / b, 1.0)
            r[k] = rr
        r[0], r[-1] = r[1], r[-2]          # under/overflow gen bins: nearest reported bin
        T_rw = (T * r).ravel()
        d_rw = (nom.A @ T_rw) / np.where(nom.fake < 1, 1 - nom.fake, 1)
        stress = {}
        for n in (1, 2, 3, 4, 6, 8, n_it):
            stress[n] = gen_shapes(nom.unfold(d_rw, n)) / gen_shapes(T_rw)
        out.update(stress=stress[n_it], stress_by_n=stress,
                   stress_regions=region_means(stress[n_it] - 1),
                   stress_hr=(half_rises(gen_shapes(nom.unfold(d_rw, n_it)), mean_pt_gen),
                              half_rises(gen_shapes(T_rw), mean_pt_gen)))
        print(f"  [{v}] stress test (truth -> Pythia 8 shape), mean (unfolded/true - 1):")
        for k in range(4):
            print(f"        bin {k}: " + "  ".join(f"{x:+.4f}" for x in out["stress_regions"][k])
                  + f"   half-rise {out['stress_hr'][0][k]:.3f} vs {out['stress_hr'][1][k]:.3f}")
        # model dependence of the data result: prior = Pythia 8-reweighted truth
        syst["prior (Pythia 6 -> 8 shape)"] = gen_shapes(nom.unfold(d, n_it, prior=T_rw)) - shp
        syst["stress-test bias"] = (stress[n_it] - 1) * shp

    # ---- tracking efficiency: response with 3% of reco tracks removed ----
    if v.startswith("ch"):
        rt = Response(S, "ch1_trk")
        syst["tracking efficiency (3%)"] = gen_shapes(
            unfold.dagostini(d, rt.A, nom.T, n_it, fake_fraction=rt.fake)) - shp
        rj = Response(S, "ch1_jer")
        syst["jet pT resolution (+10%)"] = gen_shapes(
            unfold.dagostini(d, rj.A, nom.T, n_it, fake_fraction=rj.fake)) - shp
    # ---- jet energy scale: data jet pT x (1 +- 3%) ----
    up = gen_shapes(nom.unfold(data_hist(D, v, pt_scale=1.03), n_it)) - shp
    dn = gen_shapes(nom.unfold(data_hist(D, v, pt_scale=0.97), n_it)) - shp
    syst["jet energy scale (±3%)"] = np.where(np.abs(up) > np.abs(dn), up, dn)
    # ---- response statistics: resample the 20 groups ----
    rb = []
    for _ in range(args.n_boot_resp):
        g = rng.choice(N_GROUPS, size=N_GROUPS)
        rb.append(gen_shapes(Response(S, v, g).unfold(d, n_it)))
    syst["response statistics"] = np.std(rb, axis=0, ddof=1)
    # ---- jets of gen bin 375-450 lost below the 375 GeV reco cut (absent from the SIM files) --
    # proxy: gen 450-550 jets reconstructed one bin lower, relative to those in their own bin,
    # times their share of that bin's pair weight; applied to the first bin only
    Tb = S[f"{v}_T_by_reco"].reshape(N_RECO_JET, N_GEN_JET, NR)
    a = eec_mod.shape(Tb[0, 2] / DLOG, BINS)
    b = eec_mod.shape(Tb[1, 2] / DLOG, BINS)
    share = Tb[0, 2].sum() / Tb[:, 2].sum()
    loss = np.zeros_like(shp)
    with np.errstate(divide="ignore", invalid="ignore"):
        loss[0] = np.nan_to_num(share * (a / b - 1)) * shp[0]
    syst["jets lost below 375 GeV (bin 1 only)"] = loss
    # ---- regularisation: one more iteration ----
    syst["iterations (n -> n+1)"] = gen_shapes(t_all[n_it]) - shp

    out["syst"] = syst
    # half-rise systematics: from shifted shapes
    out["hr_syst"] = {k: half_rises(np.clip(shp + s, 1e-12, None), mean_pt_gen) - hr
                      for k, s in syst.items() if k != "response statistics"}
    out["hr_syst"]["response statistics"] = np.std(
        [half_rises(x, mean_pt_gen) for x in rb], axis=0, ddof=1)
    tot_syst = np.sqrt(sum(s**2 for s in syst.values()))
    out["tot_syst"] = tot_syst
    out["tot"] = np.sqrt(tot_syst**2 + stat**2)
    out["hr_tot_syst"] = np.sqrt(sum(s**2 for s in out["hr_syst"].values()))
    print(f"  [{v}] half-rise R_L*<pT> [GeV] (stat, syst):")
    for k in range(4):
        print(f"        bin {k}: {hr[k]:.3f} ± {out['hr_stat'][k]:.3f} ± {out['hr_tot_syst'][k]:.3f}  "
              f"(bin-by-bin {out['hr_bbb'][k]:.3f}, Pythia 6 truth {out['hr_sim_gen'][k]:.3f})  | " +
              ", ".join(f"{n.split(' (')[0]} {x[k]:+.3f}" for n, x in out["hr_syst"].items()))
    summary[v] = _json_variant(out)
    return out


def _json_variant(o):
    keep = {}
    for k in ("n_iter", "iter_rows", "hr", "hr_stat", "hr_bbb", "hr_sim_gen", "closure_regions",
              "closure_chi2", "closure_hr", "stress_regions", "stress_hr", "hr_tot_syst"):
        if k in o:
            keep[k] = o[k]
    keep["hr_syst"] = o["hr_syst"]
    keep["syst_regions"] = {k: region_means(np.abs(s) / o["shape"]) for k, s in o["syst"].items()}
    keep["stat_regions"] = region_means(o["stat"] / o["shape"])
    return keep


def _tojson(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, dict):
        return {str(k): _tojson(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_tojson(v) for v in x]
    return x


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--files", nargs="+", help="SIM slices as NAME:N_FILES (required to analyse)")
    ap.add_argument("--data-files", type=int, default=3)
    ap.add_argument("--n-boot", type=int, default=100, help="data bootstrap replicas")
    ap.add_argument("--n-boot-resp", type=int, default=50, help="response-group resamples")
    ap.add_argument("--n-iter", type=int, default=0, help="fix iterations (0 = stopping rule)")
    ap.add_argument("--n-iter-max", type=int, default=20)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pythia-ch1", nargs="+", default=[
        str(REPO / "results" / "tables" / f"eec_pythia_mod2011_ch1_s{s}.npz")
        for s in (12345, 12346, 12347, 12348)],
        help="Pythia 8 charged pT > 1 GeV tables (averaged)")
    ap.add_argument("--process-only", nargs="+", help="just fill the cache for these files")
    args = ap.parse_args()
    if args.process_only:
        for f in args.process_only:
            print(process_file(f))
        return 0
    if not args.files:
        ap.error("--files is required")
    rng = np.random.default_rng(args.seed)

    print("[1/4] SIM and data -> pair-level response / per-jet histograms (cached) ...")
    S, D = load(args.files, args.data_files, args.workers)
    j = S["jets"]
    print(f"      SIM: {len(j['w'])} jets ({j['matched'].sum()} matched), n_eff {n_eff(j['w']):.0f};"
          f" data: {len(D['weight'])} jets")
    for v in VARIANTS:
        f = S[f"{v}_frac"]
        print(f"      {v:8s} matched gen particles {f[0] / f[1]:.3f} (pT-weighted {f[2] / f[3]:.3f}),"
              f" matched reco {f[4] / f[5]:.3f} (pT-weighted {f[6] / f[7]:.3f})")
    mpt_gen = mean_gen_pt(S)
    mpt_reco = mean_reco_pt(D)
    print("      <gen pT> per gen bin: " + ", ".join(f"{x:.1f}" for x in mpt_gen)
          + "; data <reco pT> per reco bin: " + ", ".join(f"{x:.1f}" for x in mpt_reco))
    summary: dict = {"seed": args.seed, "n_boot": args.n_boot, "n_boot_resp": args.n_boot_resp,
                     "files": args.files, "mean_gen_pt": mpt_gen, "mean_reco_pt": mpt_reco,
                     "matching": {v: S[f"{v}_frac"] for v in VARIANTS}}
    py8 = None
    if all(pathlib.Path(x).exists() for x in args.pythia_ch1):
        py8 = pythia8_shapes(args.pythia_ch1)
        p8s = np.array([eec_mod.shape(x, BINS) for x in py8[0]])
        summary["pythia8_hr"] = half_rises(p8s, py8[2])
        reps = [half_rises(np.array([eec_mod.shape(x, BINS) for x in
                                     py8[0] + py8[1] * rng.standard_normal(py8[0].shape)]), py8[2])
                for _ in range(args.n_boot)]
        summary["pythia8_hr_err"] = np.std(reps, axis=0, ddof=1)
        summary["pythia8_n_eff"] = py8[3]
        print("      Pythia 8 (charged, pT > 1 GeV) half-rise: " + ", ".join(
            f"{a:.3f} ± {b:.3f}" for a, b in zip(summary["pythia8_hr"], summary["pythia8_hr_err"])))
    else:
        print(f"      (no Pythia 8 charged table at {args.pythia_ch1}: stress test skipped)")

    # jet-pT migration across a bin boundary (proxy for jets lost below 375 GeV)
    Tb = S["ch1_T_by_reco"].reshape(N_RECO_JET, N_GEN_JET, NR)
    a = eec_mod.shape(Tb[0, 2] / DLOG, BINS)        # gen 450-550, reco 375-450 (migrated down)
    b = eec_mod.shape(Tb[1, 2] / DLOG, BINS)        # gen 450-550, reco 450-550
    summary["migration_proxy"] = region_means((a / b)[None, :] - 1)[0]
    frac_down = Tb[0, 2].sum() / Tb[:, 2].sum()
    print("      gen 450-550 jets reconstructed in 375-450 vs in 450-550, shape ratio - 1 per "
          "region: " + "  ".join(f"{x:+.4f}" for x in summary["migration_proxy"])
          + f" (their share of the bin's pair weight {frac_down:.3f})")
    summary["migration_share"] = frac_down

    print("[2/4] charged particles, pT > 1 GeV (main result) ...")
    res_ch = run_variant(S, D, "ch1", args, rng, py8, mpt_gen, summary, "charged, pT > 1 GeV")
    print("[3/4] all PFCs, pT > 1 GeV (comparison) ...")
    res_all = run_variant(S, D, "all1", args, rng, None, mpt_gen, summary, "all particles, pT > 1 GeV")

    tables = REPO / "results" / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    with open(tables / "mod_unfolding_summary.json", "w") as f:
        json.dump(_tojson(summary), f, indent=1)
    np.savez(tables / "eec_mod2011_unfolded.npz", edges=BINS, pt_edges=np.array(PT_EDGES),
             mean_gen_pt=mpt_gen,
             **{f"{v}_{k}": r[k] for v, r in (("ch1", res_ch), ("all1", res_all))
                for k in ("shape", "stat", "tot_syst", "tot", "reco_shape", "sim_gen", "bbb")})
    print("[4/4] figures ...")
    _plots(res_ch, res_all, py8, mpt_gen, mpt_reco, summary)
    return 0


def _plots(rc, ra, py8, mpt_gen, mpt_reco, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_mod.use_style("CMS")
    outdir = REPO / "results" / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    c = CENTRES
    labs = [f"{lo:.0f}-{hi:.0f} GeV" for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:])]
    p8s = None if py8 is None else np.array([eec_mod.shape(x, BINS) for x in py8[0]])

    # 1) unfolded charged EEC per pT bin, with Pythia 6 / 8 truth and ratio to Pythia 8
    fig, axs = plt.subplots(2, 4, figsize=(18, 7.5), sharex=True,
                            gridspec_kw={"height_ratios": [2.2, 1]})
    for k in range(4):
        ax, axr = axs[0, k], axs[1, k]
        s, tot, st = rc["shape"][k], rc["tot"][k], rc["stat"][k]
        ax.fill_between(c, s - tot, s + tot, color="k", alpha=0.18, lw=0, label="stat. + syst. (quadrature)")
        ax.errorbar(c, s, yerr=st, fmt="o", color="k", ms=3, label="data, unfolded")
        ax.plot(c, rc["reco_shape"][k], color="0.5", ls="--", lw=1, label="data, detector level")
        ax.plot(c, rc["sim_gen"][k], color="C2", lw=1.6, label="Pythia 6 Z2 (SIM gen)")
        if p8s is not None:
            ax.plot(c, p8s[k], color="C1", lw=1.6, ls=":", label="Pythia 8.318")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 0.6)
        ax.set_title(labs[k], fontsize=12)
        ref = p8s[k] if p8s is not None else rc["sim_gen"][k]
        axr.fill_between(c, (s - tot) / ref, (s + tot) / ref, color="k", alpha=0.18, lw=0)
        axr.errorbar(c, s / ref, yerr=st / ref, fmt="o", color="k", ms=2.5)
        axr.plot(c, rc["sim_gen"][k] / ref, color="C2", lw=1.4)
        axr.axhline(1, color="C1", ls=":")
        axr.set_ylim(0.5, 1.5)
        axr.set_xlabel(r"$R_L$")
    axs[0, 0].set_ylabel("charged EEC, unit area", fontsize=14)
    axs[1, 0].set_ylabel("ratio to\nPythia 8" if p8s is not None else "ratio to\nPythia 6",
                         fontsize=13)
    axs[0, 0].legend(fontsize=9, loc="lower center")
    fig.suptitle(r"CMS 2011A, anti-$k_T$ R = 0.5, $|\eta| < 1.9$; charged particles "
                 r"$p_T > 1$ GeV, unfolded to particle level", fontsize=13)
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / "eec_mod2011_unfolded_charged.png"))
    plt.close(fig)

    # 2) half-rise vs pT
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.2))
    x = np.array([(lo + hi) / 2 for lo, hi in zip(PT_EDGES[:-1], PT_EDGES[1:])])
    tot = np.sqrt(rc["hr_stat"] ** 2 + rc["hr_tot_syst"] ** 2)
    a1.errorbar(x, rc["hr"], yerr=tot, fmt="o", color="k", capsize=3, ms=6,
                label="data, unfolded (stat. + syst.)")
    a1.errorbar(x, rc["hr"], yerr=rc["hr_stat"], fmt="none", color="k", capsize=0, lw=3)
    a1.plot(x + 8, rc["hr_bbb"], "s", color="C0", mfc="none", label="data, bin-by-bin")
    a1.plot(x, rc["hr_sim_gen"], "^-", color="C2", label="Pythia 6 Z2 (SIM gen)")
    if "pythia8_hr" in summary:
        a1.errorbar(x, summary["pythia8_hr"], yerr=summary["pythia8_hr_err"], fmt="D:",
                    color="C1", capsize=2, label="Pythia 8.318")
    hr_reco = half_rises(rc["reco_shape"], mpt_reco)
    a1.plot(x, hr_reco, "v--", color="0.5", label="data, detector level (reco pT)")
    a1.set_xlabel(r"jet $p_T$ bin centre [GeV]")
    a1.set_ylabel(r"half-rise $R_L \cdot \langle p_T \rangle$ [GeV]")
    a1.set_title("charged particles, pT > 1 GeV", fontsize=12)
    a1.legend(fontsize=9)
    for k in range(4):
        a2.plot(c * mpt_gen[k], rc["shape"][k], "o-", ms=3, lw=1, label=labs[k])
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlim(0.3, 300)
    a2.set_ylim(3e-3, 0.6)
    a2.set_xlabel(r"$R_L \cdot \langle p_T^{gen} \rangle$ [GeV]")
    a2.set_ylabel("unfolded charged EEC, unit area")
    a2.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / "eec_mod2011_unfolded_halfrise.png"))
    plt.close(fig)

    # 3) validation: closure and stress test per pT bin
    fig, axs = plt.subplots(2, 4, figsize=(18, 6.5), sharex=True, sharey="row")
    for k in range(4):
        ax = axs[0, k]
        ax.fill_between(c, 1 - rc["closure_err"][k], 1 + rc["closure_err"][k], color="0.8",
                        lw=0, label="half-sample stat.")
        for n, col in ((1, "C3"), (rc["n_iter"], "k")):
            ax.plot(c, rc["closure_by_n"][n][k], color=col, lw=1.3, label=f"n = {n}")
        ax.set_ylim(0.9, 1.1)
        ax.axhline(1, color="grey", ls=":")
        ax.set_title(labs[k], fontsize=12)
        if "stress_by_n" in rc:
            ax2 = axs[1, k]
            for n, col in ((1, "C3"), (2, "C0"), (rc["n_iter"], "k"), (8, "C4")):
                if n in rc["stress_by_n"]:
                    ax2.plot(c, rc["stress_by_n"][n][k], color=col, lw=1.3, label=f"n = {n}")
            ax2.axhline(1, color="grey", ls=":")
            ax2.set_ylim(0.9, 1.1)
            ax2.set_xlabel(r"$R_L$")
            ax2.set_xscale("log")
    axs[0, 0].set_ylabel("closure:\nunfolded / truth")
    axs[1, 0].set_ylabel("stress test:\nunfolded / true (P8 shape)")
    axs[0, 0].legend(fontsize=9)
    axs[1, 0].legend(fontsize=9)
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / "eec_mod2011_unfolding_validation.png"))
    plt.close(fig)

    # 4) relative uncertainties per source
    fig, axs = plt.subplots(1, 4, figsize=(18, 4.6), sharey=True)
    for k in range(4):
        ax = axs[k]
        s = rc["shape"][k]
        ax.plot(c, rc["stat"][k] / s, color="k", lw=2, label="data statistics")
        for i, (n, v) in enumerate(rc["syst"].items()):
            ax.plot(c, np.abs(v[k]) / s, lw=1.3, color=f"C{i}", label=n)
        ax.plot(c, rc["tot"][k] / s, color="k", ls="--", lw=1.5, label="total")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 0.3)
        ax.set_xlabel(r"$R_L$")
        ax.set_title(labs[k], fontsize=12)
    axs[0].set_ylabel("relative uncertainty")
    axs[-1].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / "eec_mod2011_unfolding_uncertainties.png"))
    plt.close(fig)

    # 5) all PFCs vs charged, unfolded
    fig, ax = plt.subplots(figsize=(7, 5))
    for k in (0, 3):
        ax.plot(c, rc["shape"][k], color=f"C{k}", label=f"charged, {labs[k]}")
        ax.plot(c, ra["shape"][k], color=f"C{k}", ls="--", label=f"all PFCs, {labs[k]}")
        ax.plot(c, ra["sim_gen"][k], color=f"C{k}", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(3e-3, 0.6)
    ax.set_xlabel(r"$R_L$")
    ax.set_ylabel("unfolded EEC, unit area (dotted: Pythia 6, all)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    plot_mod.save(fig, str(outdir / "eec_mod2011_unfolded_all_vs_charged.png"))
    plt.close(fig)
    print("      figures -> results/figures/eec_mod2011_unfold*.png")


if __name__ == "__main__":
    raise SystemExit(main())
