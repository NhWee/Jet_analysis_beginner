# Jet Analysis — Setup

**Goal:** understand jet properties and the strong interaction through statistical analysis of
collider data. Concretely: **jet reconstruction** + **Energy-Energy Correlators (EEC)** on LHC
Open Data, with a deep-learning component. Research notes (incl. failure reports) live in Notion;
setup and code are documented here for GitHub.

> 📌 **처음 오셨다면 [`docs/WORKFLOW.md`](docs/WORKFLOW.md) 부터 읽어주세요.**
> 저장소 구조, 재현 방법, 커밋·기록 규칙, 지금까지의 판단 근거, 알려진 함정이 정리돼 있습니다.

---

## 1. Architecture — two layers

Docker and Conda are **complementary, not competing**:

```
┌─────────────────────────────────────────────────────────────┐
│  Layer 1 — DATA ACCESS  (Docker, reproducible)               │
│  Experiment's official image → download / skim Open Data     │
│  Needed as a framework ONLY for CMS AOD (CMSSW). Flat        │
│  ntuples (NanoAOD / ATLAS ntuple) skip this entirely.        │
│  Output: flat ROOT / parquet with particle-level 4-vectors   │
└───────────────────────────┬─────────────────────────────────┘
                            │  flat ntuples (constituents)
                            ▼
┌─────────────────────────────────────────────────────────────┐
│  Layer 2 — ANALYSIS + ML  (Conda env `jet-eec`)             │
│  uproot → awkward → fastjet (cluster) → EEC → hist → plot    │
│  PyTorch for deep learning.  This env is experiment-agnostic.│
└─────────────────────────────────────────────────────────────┘
```

The **same** `environment.yml` / `Dockerfile` in this repo cover all four LHC experiments — the
only thing that changes per experiment is how you fetch the data in Layer 1.

---

## 2. Required programs (checklist)

### System-level (install once, outside the env)

| Program | Why |
|---|---|
| **Miniforge / Mambaforge** (conda + mamba) | environment manager — your chosen approach; handles compiled HEP packages cleanly |
| **git** | version control → GitHub documentation goal |
| **Docker** | reproducible container + data-access layer for Open Data |
| **nvidia-container-toolkit** *(only for GPU)* | lets Docker use the GPU for deep learning |

### Inside the `jet-eec` conda env (see `environment.yml`)

| Package | Role |
|---|---|
| `numpy` `scipy` `pandas` `matplotlib` | core scientific stack |
| `numba` | JIT-accelerate the O(N²) EEC constituent-pair loop |
| `uproot` | read ROOT / NanoAOD files **without** installing ROOT |
| `awkward` | jagged arrays for variable-length jet constituents |
| `vector` | Lorentz 4-vectors (pt, η, φ, E) |
| `fastjet` (scikit-hep) | **jet clustering** (anti-kₜ, kₜ, C/A) + built-in **energy correlators**, awkward interface |
| `hist` `boost-histogram` `mplhep` | histogramming + publication-style plots |
| `particle` `hepunits` | PDG particle properties / units |
| `coffea` | columnar analysis framework — scales the pipeline across many files |
| `xrootd` | stream Open Data over `root://` without downloading everything |
| `energyflow` | EFPs, EFN/PFN networks, EMD, ready-made jet datasets — jets ∩ ML ∩ correlators |
| `pytorch` | deep learning (CPU by default; CUDA variant for GPU) |
| `scikit-learn` | baselines, preprocessing, metrics |
| `pythia8` | event generator — build & validate the EEC pipeline on **truth** particles first |
| `jupyterlab` `ipykernel` | notebooks |

**Optional:** `root` (conda-forge) if you want PyROOT / RDataFrame — installs cleanly, unlike pip.

---

## 3. Setup

### Option A — Conda (local)

```bash
# 1. install Miniforge (provides conda + mamba), then:
mamba env create -f environment.yml
conda activate jet-eec

# 2. sanity check
python -c "import uproot, awkward, vector, fastjet, energyflow, torch; print('ok')"

# 3. freeze exact versions for reproducibility (commit this file)
conda env export --no-builds > environment.lock.yml
```

**GPU:** the env ships CPU PyTorch. For CUDA, after creating the env:

```bash
mamba install -n jet-eec pytorch pytorch-cuda=12.4 -c pytorch -c nvidia
```

### Option B — Docker (reproducible)

```bash
docker build -t jet-eec .
docker run -it --rm -p 8888:8888 -v "$PWD":/work jet-eec      # open http://localhost:8888
# GPU: add  --gpus all  (needs nvidia-container-toolkit + a CUDA PyTorch build)
```

---

## 4. Data access — what actually exists

| Source | Public Open Data? | Notes |
|---|---|---|
| **CMS** | ✅ Yes | Largest release. NanoAOD (light, uproot-readable) → AOD (needs CMSSW Docker). |
| **ATLAS** | ✅ Yes | 13 TeV flat ntuples, lightest to start; uproot / RDataFrame. Official Jupyter Docker image. |
| **ALICE** | ✅ Yes | On the CERN Open Data portal. |
| **LHCb** | ✅ Yes | On the CERN Open Data portal (more limited but growing). |
| **RHIC** (STAR / sPHENIX) | ❌ No public portal | Collaboration access only. sPHENIX (2023–2026) is jet-focused and physically ideal for EEC, but not openly downloadable. |
| **Fermilab Tevatron** (CDF / D0) | ❌ No public portal | Data *preserved* (~9 PB each, on CVMFS) for the collaborations — not open data. |

**Practical path:** the realistic playground is the four LHC experiments via
[opendata.cern.ch](https://opendata.cern.ch/). Start with **ATLAS 13 TeV ntuples** or **CMS
NanoAOD** (light, no framework), reserve CMSSW-Docker for CMS AOD only if you need it.

> ⚠️ **EEC needs jet *constituents*** (particle-level 4-vectors inside each jet). Plain NanoAOD
> stores jets but not always all constituents — use **PFNano**-style releases or AOD (CMS), or the
> track/topocluster collections (ATLAS). This is the one place the data tier matters, so develop
> the pipeline on Pythia MC first (full truth control), then map it onto the richest Open Data tier.

---

## 5. Jet + EEC pipeline (outline)

```python
import fastjet, awkward as ak, vector
vector.register_awkward()

# 1. read constituents (uproot) -> awkward array of particle 4-vectors per event
# 2. cluster jets
jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
cluster = fastjet.ClusterSequence(particles, jetdef)
jets = cluster.inclusive_jets(min_pt=100)          # GeV
constituents = cluster.constituents()

# 3. EEC: for each jet, sum over constituent pairs (i,j)
#    weight  w_ij = (pt_i * pt_j) / pt_jet**2   binned in ΔR_ij (or angle)
#    fastjet also exposes energy-correlator helpers directly.
```

For the ML component, `energyflow` provides Energy/Particle Flow Networks (permutation-invariant
over constituents) and EMD — a natural fit alongside EEC and easy to wrap in PyTorch.

---

## 6. Repository layout

```
Jet_analysis/
├── environment.yml          # conda env definition
├── Dockerfile               # reproducible container (analysis + ML layer)
├── data/                    # NOT committed — see data/README.md sample register
│   ├── raw/                 #   as downloaded, never edited
│   ├── interim/             #   skimmed ntuples with constituents
│   └── processed/           #   analysis-ready parquet
├── src/
│   ├── jeteec/
│   │   ├── jets.py          #   clustering wrapper (fastjet, awkward)
│   │   ├── eec.py           #   EEC computation + normalisation check
│   │   └── plot.py          #   mplhep publication styling
│   ├── validate_eec_pythia.py  # run this FIRST — validates pipeline on truth
│   ├── toy_generator.py     #   dependency-free stand-in for Pythia
│   ├── tier_demo.py         #   what each data tier can/cannot do (see docs/)
│   └── tier_analyses.py     #   one real analysis per tier (see docs/)
├── notebooks/               # exploration (strip outputs before committing)
├── results/figures|tables/  # outputs, gitignored
├── models/                  # checkpoints, gitignored
└── docs/                    # notes mirrored to Notion
```

### First run

```bash
conda activate jet-eec
python src/validate_eec_pythia.py --n-events 2000 --seed 12345
```

This generates Pythia dijets, clusters anti-kT R=0.4 jets, and checks that the EEC weights
sum to 1.0 per jet. **It must PASS before any Open Data is touched** — a failure there is a
pipeline bug, which is far cheaper to find on truth particles than on detector data.

### Which data tier do you need?

```bash
python src/tier_demo.py --n-events 400 --seed 12345
```

Writes the same events at NanoAOD / PFNano / AOD structure and reports which analyses each
supports. Short answer: **EEC needs constituents, so NanoAOD is out and PFNano is the sweet
spot.** Mechanics in [`docs/data_tiers.md`](docs/data_tiers.md); what each tier lets you
*measure*, and in what order to attack them, in
[`docs/tier_physics_programme.md`](docs/tier_physics_programme.md).

### One real analysis per tier

```bash
python src/tier_analyses.py --n-events 2000 --seed 12345
```

Runs the analysis each tier actually permits — NanoAOD: jet spectrum, dijet χ, R₃₂ ·
PFNano: fragmentation, girth/pTD, EEC · MiniAOD: radius scan, trimming, charged-only ·
AOD: tracking systematics. Results and figures:
[`docs/tier_analyses_results.md`](docs/tier_analyses_results.md).

> ⚠️ **fastjet gotcha found here:** passing a record with extra fields (e.g. `charge`) to
> `ClusterSequence` silently yields wrong jets — off by up to 168 GeV in pT in our test, with
> no error raised. Always strip to `pt, eta, phi, mass` first (`jeteec.jets.four`).

---

## 7. Reproducibility notes (for GitHub / Notion)

- Commit `environment.lock.yml` (exact versions), not just `environment.yml`.
- Record the **dataset DOI / record ID** from the CERN Open Data portal for every sample used.
- Set and log random seeds (`numpy`, `torch`) in every training run.
- Keep a `failures.md` (or the Notion failure log) — mismatched constituent collections and
  jet-area / pileup subtraction are the usual early pitfalls.

---

## Sources

- [CMS Open Data — Docker guide](https://cms-opendata-guide.web.cern.ch/tools/docker/) · [Running CMS analysis with Docker](https://opendata.cern.ch/docs/cms-guide-docker)
- [ATLAS Open Data — Jupyter/notebooks](https://opendata.atlas.cern/docs/notebooks/intro) · [notebooks Docker image](https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata)
- [scikit-hep/fastjet](https://github.com/scikit-hep/fastjet) · [fastjet docs (awkward interface)](https://fastjet.readthedocs.io/en/latest/Awkward.html)
- [Data preservation at the Fermilab Tevatron (CDF/D0)](https://arxiv.org/abs/1701.07773)
