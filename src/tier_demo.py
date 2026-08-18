"""Data-tier demo — what can you actually DO with each tier?

Writes three ROOT files whose branch structure mimics real CMS Open Data tiers, then
attempts five analyses on each and reports which succeed. The point is to make the
information loss concrete: every tier stores the *same events*, just less of them.

    tier          branches stored                                  mimics
    ----          ---------------                                  ------
    NanoAOD       Jet_pt/eta/phi/mass                               CMS NanoAOD
    PFNano        + PFCand_* with PFCand_jetIdx (in-jet particles)  CMS PFNano / JMENano
    AOD           all event particles, no jets stored at all        CMS AOD (PFCandidates)

Note the AOD file stores *no jets*: at AOD level you run the clustering yourself, which
is exactly why it is the only tier where you may change R or the algorithm.

Usage:  python src/tier_demo.py --n-events 300 --seed 12345
Output: data/interim/tier_demo/{nanoaod,pfnano,aod}_like.root
        results/figures/tier_capability.png
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import awkward as ak
import numpy as np
import uproot
import vector

sys.path.insert(0, str(pathlib.Path(__file__).parent))
vector.register_awkward()

REPO = pathlib.Path(__file__).resolve().parents[1]
OUTDIR = REPO / "data" / "interim" / "tier_demo"
FIGDIR = REPO / "results" / "figures"

R_NOMINAL = 0.4          # the radius the "experiment" reconstructed jets with
JET_PT_MIN = 100.0       # GeV


# --------------------------------------------------------------------------------------
# 1. Truth sample — a toy stand-in for Pythia (kept dependency-free so this always runs)
# --------------------------------------------------------------------------------------
def make_truth_events(n_events: int, rng: np.random.Generator) -> ak.Array:
    """Two back-to-back hard jets plus soft background, as final-state particles.

    Constituent angular distribution ~ 1/theta and a steeply falling energy sharing,
    which is enough structure to give an EEC with a realistic shape.
    """
    events = []
    for _ in range(n_events):
        pt, eta, phi, mass = [], [], [], []

        jet_pt = 100.0 + rng.exponential(120.0)
        jet_eta = rng.uniform(-1.5, 1.5)
        jet_phi = rng.uniform(-np.pi, np.pi)

        for sign in (0, 1):  # two back-to-back jets
            axis_eta = jet_eta if sign == 0 else -jet_eta
            axis_phi = jet_phi if sign == 0 else (jet_phi + np.pi + rng.normal(0, 0.1))
            axis_phi = (axis_phi + np.pi) % (2 * np.pi) - np.pi
            this_pt = jet_pt * (1.0 if sign == 0 else rng.uniform(0.75, 1.0))

            n_const = rng.integers(18, 45)
            # steeply falling momentum fractions
            z = rng.dirichlet(np.full(n_const, 0.45))
            # angular spread ~ 1/theta out to the jet radius
            u = rng.uniform(0, 1, n_const)
            dr = 0.005 * np.exp(u * np.log(R_NOMINAL / 0.005))
            ang = rng.uniform(0, 2 * np.pi, n_const)

            pt.extend(this_pt * z)
            eta.extend(axis_eta + dr * np.cos(ang))
            phi.extend(axis_phi + dr * np.sin(ang))
            mass.extend(np.full(n_const, 0.139))  # pion mass

        # soft underlying event, spread over the whole acceptance
        n_soft = rng.integers(20, 60)
        pt.extend(rng.exponential(0.7, n_soft))
        eta.extend(rng.uniform(-3, 3, n_soft))
        phi.extend(rng.uniform(-np.pi, np.pi, n_soft))
        mass.extend(np.full(n_soft, 0.139))

        events.append({"pt": pt, "eta": eta, "phi": phi, "mass": mass})

    arr = ak.Array(events)
    return ak.zip(
        {"pt": arr.pt, "eta": arr.eta, "phi": arr.phi, "mass": arr.mass},
        with_name="Momentum4D",
    )


# --------------------------------------------------------------------------------------
# 2. Write the three tiers
# --------------------------------------------------------------------------------------
def cluster_nominal(particles: ak.Array):
    """Reconstruct jets the way the experiment would have, at the nominal radius."""
    import fastjet

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, R_NOMINAL)
    seq = fastjet.ClusterSequence(particles, jetdef)
    jets = seq.inclusive_jets(min_pt=JET_PT_MIN)
    const = seq.constituents(min_pt=JET_PT_MIN)
    keep = abs(jets.eta) < 2.0
    return jets[keep], const[keep]


def write_tiers(particles: ak.Array) -> dict[str, pathlib.Path]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    jets, const = cluster_nominal(particles)
    paths: dict[str, pathlib.Path] = {}

    # -- NanoAOD: jet kinematics only ---------------------------------------------------
    p = OUTDIR / "nanoaod_like.root"
    with uproot.recreate(p) as f:
        f["Events"] = {
            "Jet_pt": jets.pt,
            "Jet_eta": jets.eta,
            "Jet_phi": jets.phi,
            "Jet_mass": jets.mass,
        }
    paths["NanoAOD"] = p

    # -- PFNano: jets + the particles inside them, linked by jetIdx ---------------------
    # Real PFNano stores a flat per-event candidate list plus an index back to the jet.
    jet_idx3 = ak.broadcast_arrays(ak.local_index(jets, axis=1), const.pt)[0]
    jet_idx = ak.flatten(jet_idx3, axis=2)
    cflat = ak.flatten(const, axis=2)

    p = OUTDIR / "pfnano_like.root"
    with uproot.recreate(p) as f:
        f["Events"] = {
            "Jet_pt": jets.pt,
            "Jet_eta": jets.eta,
            "Jet_phi": jets.phi,
            "Jet_mass": jets.mass,
            "PFCand_pt": cflat.pt,
            "PFCand_eta": cflat.eta,
            "PFCand_phi": cflat.phi,
            "PFCand_mass": cflat.mass,
            "PFCand_jetIdx": ak.values_astype(jet_idx, np.int32),
        }
    paths["PFNano"] = p

    # -- AOD: every particle in the event, no jets stored -------------------------------
    p = OUTDIR / "aod_like.root"
    with uproot.recreate(p) as f:
        f["Events"] = {
            "PFCand_pt": particles.pt,
            "PFCand_eta": particles.eta,
            "PFCand_phi": particles.phi,
            "PFCand_mass": particles.mass,
        }
    paths["AOD"] = p

    return paths


# --------------------------------------------------------------------------------------
# 3. Read each tier back and probe what it supports
# --------------------------------------------------------------------------------------
def load(path: pathlib.Path) -> dict:
    """Read a tier file, returning whatever it happens to contain."""
    with uproot.open(path) as f:
        tree = f["Events"]
        branches = set(tree.keys())
        arrs = tree.arrays()

    out: dict = {"branches": sorted(branches), "jets": None, "constituents": None,
                 "particles": None}

    if "Jet_pt" in branches:
        out["jets"] = ak.zip(
            {"pt": arrs.Jet_pt, "eta": arrs.Jet_eta,
             "phi": arrs.Jet_phi, "mass": arrs.Jet_mass},
            with_name="Momentum4D",
        )

    if "PFCand_pt" in branches:
        cands = ak.zip(
            {"pt": arrs.PFCand_pt, "eta": arrs.PFCand_eta,
             "phi": arrs.PFCand_phi, "mass": arrs.PFCand_mass},
            with_name="Momentum4D",
        )
        if "PFCand_jetIdx" in branches:
            # regroup the flat candidate list back into per-jet lists via jetIdx
            order = ak.argsort(arrs.PFCand_jetIdx, axis=1)
            counts = ak.run_lengths(arrs.PFCand_jetIdx[order])
            out["constituents"] = ak.unflatten(
                cands[order], ak.flatten(counts), axis=1
            )
        else:
            out["particles"] = cands  # full event record — reclusterable

    return out


def probe(tier: str, data: dict) -> dict[str, tuple[bool, str]]:
    """Attempt five analyses; return {analysis: (ok, note)}."""
    from jeteec import eec as eec_mod

    res: dict[str, tuple[bool, str]] = {}
    jets, const, parts = data["jets"], data["constituents"], data["particles"]

    # (a) inclusive jet pT spectrum
    if jets is not None:
        res["Jet pT spectrum"] = (True, f"<pT> = {ak.mean(ak.flatten(jets.pt)):.1f} GeV")
    elif parts is not None:
        j, _ = cluster_nominal(parts)
        res["Jet pT spectrum"] = (
            True, f"<pT> = {ak.mean(ak.flatten(j.pt)):.1f} GeV (self-clustered)")
    else:
        res["Jet pT spectrum"] = (False, "no jets, no particles")

    # (b) dijet invariant mass
    def dijet(j):
        two = j[ak.num(j) >= 2]
        return (two[:, 0] + two[:, 1]).mass
    if jets is not None:
        res["Dijet mass"] = (True, f"<m_jj> = {ak.mean(dijet(jets)):.0f} GeV")
    elif parts is not None:
        j, _ = cluster_nominal(parts)
        res["Dijet mass"] = (True, f"<m_jj> = {ak.mean(dijet(j)):.0f} GeV (self-clustered)")
    else:
        res["Dijet mass"] = (False, "no jets")

    # (c) constituent multiplicity — needs particles inside jets
    if const is not None:
        res["Constituent multiplicity"] = (
            True, f"<n> = {ak.mean(ak.flatten(ak.num(const, axis=2))):.1f}")
    elif parts is not None:
        _, c = cluster_nominal(parts)
        res["Constituent multiplicity"] = (
            True, f"<n> = {ak.mean(ak.flatten(ak.num(c, axis=2))):.1f} (self-clustered)")
    else:
        res["Constituent multiplicity"] = (False, "jet 4-vectors only")

    # (d) EEC — the target observable
    c_use = const
    if c_use is None and parts is not None:
        _, c_use = cluster_nominal(parts)
    if c_use is not None:
        flat = ak.flatten(c_use, axis=1)
        norm = eec_mod.normalisation_check(flat)
        res["EEC"] = (True, f"{len(flat)} jets, norm = {norm:.4f}")
    else:
        res["EEC"] = (False, "needs constituents")

    # (e) recluster at a different radius / groom — needs the full event
    if parts is not None:
        j04, _ = cluster_nominal(parts)  # nominal, for comparison
        import fastjet
        seq = fastjet.ClusterSequence(
            parts, fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8))
        j08 = seq.inclusive_jets(min_pt=JET_PT_MIN)
        # a wider cone sweeps up more soft radiation, so <pT> rises
        res["Recluster at R=0.8"] = (
            True,
            f"<pT> = {ak.mean(ak.flatten(j08.pt)):.1f} GeV at R=0.8 vs "
            f"{ak.mean(ak.flatten(j04.pt)):.1f} GeV at R=0.4",
        )
    elif const is not None:
        res["Recluster at R=0.8"] = (
            False, "only in-jet particles kept; R=0.8 needs everything outside R=0.4 too")
    else:
        res["Recluster at R=0.8"] = (False, "no particles at all")

    return res


# --------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    print(f"Generating {args.n_events} toy events (seed={args.seed}) ...")
    particles = make_truth_events(args.n_events, rng)
    print(f"  {len(particles)} events, "
          f"<particles/event> = {ak.mean(ak.num(particles)):.0f}\n")

    print("Writing tier files ...")
    paths = write_tiers(particles)
    sizes = {t: p.stat().st_size for t, p in paths.items()}
    for tier, p in paths.items():
        per_evt = sizes[tier] / args.n_events
        print(f"  {tier:<9} {p.name:<20} {sizes[tier]/1024:8.1f} kB "
              f"({per_evt:7.0f} B/event)")
    print()

    analyses = ["Jet pT spectrum", "Dijet mass", "Constituent multiplicity",
                "EEC", "Recluster at R=0.8"]
    tiers = ["NanoAOD", "PFNano", "AOD"]
    results, loaded = {}, {}
    for tier in tiers:
        loaded[tier] = load(paths[tier])
        results[tier] = probe(tier, loaded[tier])

    # -- capability table ---------------------------------------------------------------
    print("Branches stored")
    for tier in tiers:
        print(f"  {tier:<9} {', '.join(loaded[tier]['branches'])}")
    print()

    w = max(len(a) for a in analyses) + 2
    print("Capability matrix")
    print("  " + "".ljust(w) + "".join(t.ljust(12) for t in tiers))
    for a in analyses:
        row = "  " + a.ljust(w)
        for tier in tiers:
            row += ("YES" if results[tier][a][0] else "NO").ljust(12)
        print(row)
    print()

    print("Detail")
    for tier in tiers:
        print(f"  [{tier}]")
        for a in analyses:
            ok, note = results[tier][a]
            print(f"     {'+' if ok else '-'} {a:<26} {note}")
    print()

    # -- figure: EEC from the tiers that can do it -------------------------------------
    from jeteec import eec as eec_mod
    from jeteec import plot as plot_mod
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_mod.use_style("ATLAS")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)

    for tier, style, lw in [("PFNano", "-", 3.5), ("AOD", "--", 2.0)]:
        d = loaded[tier]
        c = d["constituents"]
        if c is None:
            _, c = cluster_nominal(d["particles"])
        centres, values = eec_mod.eec(ak.flatten(c, axis=1))
        plot_mod.plot_eec(centres, values, label=tier, ax=ax1, ls=style, lw=lw)
    ax1.axvline(R_NOMINAL, color="grey", ls=":", lw=1.5)
    ax1.text(R_NOMINAL * 0.92, 5e-1, "R = 0.4", color="grey",
             fontsize=13, ha="right")
    ax1.text(1.2e-3, 3e-1, "NanoAOD: not computable",
             fontsize=13, color="#C44E52")
    ax1.set_title("EEC needs constituents\n(PFNano = AOD: consistency check)",
                  fontsize=15)
    ax1.legend(fontsize=13, loc="lower right")

    sizes_kb = [sizes[t] / 1024 for t in tiers]
    bars = ax2.bar(tiers, sizes_kb, color=["#4C72B0", "#DD8452", "#C44E52"])
    ax2.set_ylabel("file size [kB]", fontsize=14)
    ax2.set_title(f"Storage cost, same {args.n_events} events", fontsize=15)
    ax2.tick_params(labelsize=13)
    ax2.set_ylim(0, max(sizes_kb) * 1.18)
    for b, s in zip(bars, sizes_kb):
        ax2.text(b.get_x() + b.get_width() / 2, s, f"{s:.0f} kB",
                 ha="center", va="bottom", fontsize=13)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / "tier_capability.png"
    plot_mod.save(fig, str(out))
    print(f"Figure -> {out.relative_to(REPO)}")

    # -- headline conclusion ------------------------------------------------------------
    print("\nTakeaway")
    print("  NanoAOD : jet kinematics only        -> no EEC, no substructure")
    print("  PFNano  : + in-jet particles         -> EEC works; radius is frozen at R=0.4")
    print("  AOD     : every particle in the event-> full freedom, largest files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
