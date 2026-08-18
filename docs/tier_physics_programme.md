# What to do at each data tier

Project goal: **understand jet properties and the strong interaction through statistical
analysis of collider data.** This document maps that goal onto each CMS data tier — what is
physically reachable, what it costs, and what it cannot do.

A correction to the earlier working assumption is recorded in §0; read it first.

---

## 0. Important: what CMS Open Data actually ships

| Format | On the Open Data portal? | Readable with | Constituents |
|---|---|---|---|
| **NanoAOD** | ✅ 2016 onward, alongside MiniAOD | **uproot** (plain C++ types) | ✗ |
| **MiniAOD** | ✅ 2015 onward | **CMSSW required** | ✅ all, compressed |
| **AOD** | ✅ mainly Run 1 (2010–2012) | **CMSSW required** | ✅ full detail |
| **RECO** | ✗ not released | — | — |
| **PFNano** | ❌ **not a released format** | uproot, once you make it | ✅ in-jet |

**PFNano is not something you download.** It is a CMSSW framework
([cms-jet/PFNano](https://github.com/cms-jet/PFNano)) that you run *yourself* over MiniAOD to
produce NanoAOD-like files that additionally carry jet constituents. So the earlier phrasing
"PFNano is the sweet spot, no CMSSW needed" was wrong on the second half: reaching it requires
a CMSSW pass. It remains the right *destination* — just budget for the production step.

Practical consequence: **any EEC analysis on CMS Open Data requires CMSSW at least once.**
Either you analyse MiniAOD directly inside CMSSW, or you use CMSSW once to convert MiniAOD →
PFNano and then work in pure Python forever after. The second is almost always better.

Size scaling: NanoAOD is roughly **20× smaller than MiniAOD** and **200× smaller than AOD**.

---

## 1. NanoAOD — jet kinematics and summary substructure

**Read with:** uproot, no framework. **~1–2 kB/event.** Start here regardless of ambition.

### What it stores

Per-jet summary numbers, not particles: `Jet_pt/eta/phi/mass`, energy fractions
(`Jet_chHEF`, `Jet_neHEF`, `Jet_chEmEF`, ...), b-tag discriminants, and — importantly for us —
a **FatJet** collection with *pre-computed* substructure: `FatJet_msoftdrop`,
`FatJet_tau1..4`, and ParticleNet/DeepAK8 taggers.

> Check the branch list of your specific file; contents vary by NanoAOD version.

### Analyses that serve the goal

**1a. Inclusive jet pT spectrum in rapidity bins.**
The classic QCD measurement. The falling spectrum and its rapidity dependence test pQCD
directly and are sensitive to both αs and the gluon PDF at high x. Doable end-to-end in a
week with uproot + hist.

**1b. Dijet mass and dijet angular distributions.**
The χ = exp(|y₁−y₂|) distribution is a beautiful strong-interaction observable: Rutherford
scattering gives a flat χ distribution, and deviations measure the QCD matrix element. Also
the standard place to set limits on contact interactions.

**1c. Jet-multiplicity ratios (R₃₂ = N₃jet/N₂jet).**
Historically an αs extraction channel — each extra jet costs one power of αs. Reachable with
jet counting alone, so NanoAOD is enough.

**1d. Summary-level jet substructure with FatJets.**
τ₂₁ and τ₃₂ separate 1-, 2-, 3-prong jets; soft-drop mass gives the groomed mass peak. You can
study the QCD jet mass spectrum and quark/gluon differences without ever seeing a constituent.
This is *real* substructure physics — just not correlator physics.

**1e. Deep learning at jet-summary level.**
Train on the ~20 jet-level features (kinematics + energy fractions + tags) for quark/gluon
discrimination or jet-flavour classification. Modest performance ceiling, but the full ML
workflow — splits, calibration, systematics, overtraining checks — is identical to what you
will need later. Build the muscle here where iteration takes seconds.

### What it cannot do

**No EEC, no custom clustering, no grooming you did not inherit, no per-particle anything.**
The constituents are simply absent.

---

## 2. MiniAOD — the first tier where the real target opens up

**Read with:** CMSSW (Docker image). **~30–50 kB/event.** 2015 onward on the portal.

### What it stores

`packedPFCandidates`: **every particle-flow candidate in the event**, with compressed
kinematics, reduced precision, and limited PID — plus charge and vertex association. Also full
jet collections, MET, leptons, and trigger information.

The key point is *every* candidate, not just in-jet ones. That is strictly more than PFNano
gives you.

### Analyses that serve the goal

**2a. Two-point energy correlator, E2C — the project's target observable.**

$$\text{EEC}(R_L) = \sum_{i,j} \frac{p_{T,i}\,p_{T,j}}{p_{T,\text{jet}}^2}\,\delta(R_L - \Delta R_{ij})$$

The distribution has three regions with distinct physics: a free-hadron region at small angle,
a **scaling region** governed by perturbative QCD, and a turnover near R. The transition
between them images the confinement transition directly.

**2b. Three-point correlator E3C, and the ratio E3C/E2C → extract αs.**
This is the sharpest form of "understand the strong interaction" available here. CMS did
exactly this and obtained αs(mZ) = 0.1229 +0.0040 −0.0050 — the most precise αs from jet
substructure to date ([arXiv:2402.13864](https://arxiv.org/abs/2402.13864); data on
[HEPData](https://www.hepdata.net/record/ins2760466)). Their setup: anti-kT R=0.4, jet pT from
97 GeV to 1.8 TeV. **You can reproduce this on Open Data.** The ratio is powerful because
normalisation, jet energy scale, and much of the non-perturbative physics cancel.

**2c. Track-based ("charged-only") correlators.**
packedPFCandidates carry charge, so you can build EECs from charged particles only. Tracks have
far better angular resolution than calorimeter deposits, which matters enormously at small R_L
where the interesting physics lives. Cost: a non-perturbative track-function correction.

**2d. Reclustering at any R and any algorithm; grooming.**
Because all candidates are present, you are free to run anti-kT R=0.8, C/A for Lund planes,
soft-drop with your own z_cut and β, trimming, and so on. PFNano cannot do this — it kept only
what was inside R=0.4.

**2e. R-dependence and pT-dependence of the EEC scaling region.**
Scanning jet pT moves the scaling window; that running *is* the running of αs made visible.

**2f. Deep learning on constituents.**
Particle Flow Networks / Energy Flow Networks, ParticleNet-style graph nets, transformers over
constituents. This is where jet ML actually lives, and it needs exactly what MiniAOD provides.

**2g. Produce your own PFNano here.** One CMSSW pass, then pure Python thereafter.

### What it cannot do

Reduced numerical precision and pT thresholds on stored candidates limit very-low-pT and
very-small-angle work. Full track parameters, hits, and detailed reconstruction internals are
gone — so you inherit CMS's reconstruction and cannot study or re-derive it.

---

## 3. AOD — full reconstruction detail

**Read with:** CMSSW. **~200 kB/event.** Mostly Run 1 (2010–2012) on the portal.

### What it stores

Full `PFCandidate` objects with uncompressed kinematics, complete track collections with fit
parameters and quality flags, vertices, calorimeter clusters, and the full reconstruction
provenance.

### Analyses that serve the goal

**3a. Small-angle EEC with controlled tracking systematics.**
The most interesting EEC region is small R_L, precisely where tracking resolution and
efficiency dominate the uncertainty. AOD lets you apply your own track quality cuts, study
efficiency vs. angular separation, and quantify two-track resolution effects. If you want the
free-hadron region to be *believable*, this is where that work happens.

**3b. Jet energy corrections and their uncertainty, from scratch.**
Derive JEC yourself, or at least study how the correction propagates into the observable. At
MiniAOD you take JEC on faith.

**3c. Particle-flow composition studies.**
How the EEC differs for charged vs neutral hadrons vs photons, with full PID rather than
MiniAOD's compressed version.

**3d. Run 1 vs Run 2 comparison — collision energy dependence.**
Since Run 1 (7/8 TeV) is where AOD lives and Run 2 (13 TeV) is MiniAOD/NanoAOD, comparing the
EEC across √s is a genuine physics result about how the scaling region evolves.

### What it cannot do / costs

Large images, large files, slow. Realistically 10–100× the CPU and storage of the MiniAOD path
for a small improvement in most observables. **Only go here when a specific systematic forces
you to.**

---

## 4. RECO — detector-level

**~3 MB/event. Not released as Open Data**, and CMS itself barely uses it for analysis.

Contents are hits, clusters, and full reconstruction intermediates. The physics it enables is
*detector* physics: calibration, alignment, reconstruction-algorithm development, tracking
efficiency from first principles.

**For this project: not applicable, and not obtainable.** Listed only for completeness of the
chain RAW → RECO → AOD → MiniAOD → NanoAOD.

---

## 5. Recommended progression

```
  Step 0   Pythia truth          validate_eec_pythia.py       — days
           No data. Prove the EEC code is right where you control everything.

  Step 1   NanoAOD (uproot)      inclusive jet pT, dijet χ    — 2–4 weeks
           Real data, no framework. Learn selection, triggers, JEC, unfolding
           basics, and get a QCD measurement finished end-to-end.

  Step 2   MiniAOD → PFNano      one CMSSW pass in Docker     — 1–2 weeks
           The single unavoidable framework step. Output: uproot-readable
           files with constituents.

  Step 3   PFNano (uproot)       E2C, then E3C/E2C → αs       — the main project
           Pure Python from here. This is where the goal is actually met.

  Step 4   AOD                   only if tracking systematics — optional
           block a claim you want to make at small R_L.
```

The reason to respect this order: each step's failures are cheap to diagnose only if the
previous step is known-good. Debugging an EEC anomaly is very different when you already know
your clustering, selection, and corrections reproduce a published jet spectrum.

---

## 6. Tier → capability summary

| | NanoAOD | PFNano | MiniAOD | AOD | RECO |
|---|:---:|:---:|:---:|:---:|:---:|
| Jet pT / dijet spectra | ✅ | ✅ | ✅ | ✅ | ✅ |
| FatJet τ-ratios, softdrop mass | ✅ | ✅ | ✅ | ✅ | ✅ |
| Constituent multiplicity | ✗ | ✅ | ✅ | ✅ | ✅ |
| **EEC / E2C / E3C** | ✗ | ✅ | ✅ | ✅ | ✅ |
| Recluster at new R, custom grooming | ✗ | ✗ | ✅ | ✅ | ✅ |
| Charged-only correlators | ✗ | ~ | ✅ | ✅ | ✅ |
| Own track-quality cuts, tracking systematics | ✗ | ✗ | ✗ | ✅ | ✅ |
| Derive JEC yourself | ✗ | ✗ | ✗ | ✅ | ✅ |
| Detector calibration / alignment | ✗ | ✗ | ✗ | ~ | ✅ |
| **Needs CMSSW** | no | to produce | yes | yes | yes |
| **Size / event** | 1–2 kB | ~10–20 kB | 30–50 kB | ~200 kB | ~3 MB |
| **On Open Data portal** | ✅ 2016+ | ✗ make it | ✅ 2015+ | ✅ Run 1 | ✗ |

`~` = partially, with caveats.

---

## 7. Epistemic notes

- Sizes are the standard quoted figures and vary with dataset, era, and event content; treat
  them as order-of-magnitude.
- NanoAOD branch content **differs between versions** (v9 vs v12 etc.). Verify against the
  actual file rather than this table.
- The αs value quoted is CMS's published result, listed as a reproduction target — not
  something this project has verified.
- ATLAS/ALICE/LHCb use different tier names and different Open Data formats. The *logic* here
  transfers; the branch names do not.

---

## Sources

- [CMS Open Data — about CMS data formats](https://opendata.cern.ch/docs/about-cms)
- [Getting started with CMS NanoAOD Open Data](https://opendata.cern.ch/docs/cms-getting-started-nanoaod) · [MiniAOD](https://opendata.cern.ch/docs/cms-getting-started-miniaod)
- [First CMS open data from LHC Run 2](https://cms.cern/news/first-cms-open-data-lhc-run-2-released)
- [cms-jet/PFNano](https://github.com/cms-jet/PFNano)
- Petrucciani, Rizzi, Vuosalo, *Mini-AOD: A New Analysis Data Format for CMS*, [arXiv:1702.04685](https://arxiv.org/abs/1702.04685)
- CMS, *Measurement of energy correlators inside jets and determination of αS(mZ)*, [arXiv:2402.13864](https://arxiv.org/abs/2402.13864) · [HEPData](https://www.hepdata.net/record/ins2760466)
- Larkoski, Moult, Neill, *Energy correlation functions for jet substructure*, [arXiv:1305.0007](https://arxiv.org/abs/1305.0007)
