# Data

Data files are **not committed** (see `.gitignore`). Instead, every sample used must be recorded
here with its persistent identifier so the analysis is reproducible from this repo alone.

## Layout

| Dir | Contents |
|---|---|
| `raw/` | As downloaded from the source. Never edited. |
| `interim/` | Skimmed / flattened ntuples (constituent 4-vectors), Pythia output. |
| `processed/` | Analysis-ready arrays (parquet) fed to clustering / EEC / training. |

## Sample register

Record every dataset here. Template:

| ID | Experiment | Dataset / record | DOI or record URL | Tier | Constituents? | Used in |
|----|-----------|------------------|-------------------|------|---------------|---------|
| P8-001 | Pythia8 8.318 (`pythia8mc`) | HardQCD:all, pp 13 TeV, pTHatMin=100 GeV, 10k events | n/a (generated, seed=12345) | truth, hadron level | yes | EEC validation, [`docs/pythia_validation_results.md`](../docs/pythia_validation_results.md) |
| P8-002 | Pythia8 8.318 (`pythia8mc`) | as P8-001 with `HadronLevel:all = off` | n/a (generated, seed=12345) | truth, parton level | yes | hadronisation comparison |
| P8-003 | Pythia8 8.318 (`pythia8mc`) | HardQCD:all, pp 7 TeV, pTHatMin=150 GeV, bias2Selection pow 4 ref 150, 40k weighted events; `--preset mod2011` | n/a (generated, seed=12345) | truth, hadron level | yes | MOD-matched baseline, [`docs/pythia_validation_results.md`](../docs/pythia_validation_results.md) §4 |
| P8-004 | Pythia8 8.318 (`pythia8mc`) | as P8-003 with `HadronLevel:all = off` | n/a (generated, seed=12345) | truth, parton level | yes | MOD-matched hadronisation comparison |
| CMS-001 *(candidate, not downloaded)* | CMS | 2011A Jet primary dataset, Jet300, MOD HDF5 (`CMS2011AJets`, `cms`) | [10.5281/zenodo.3340205](https://doi.org/10.5281/zenodo.3340205) | reprocessed AOD | yes (PFCs) | planned: first real-data EEC |

**Reminder:** EEC requires jet *constituents*. Before registering a sample, confirm it stores
particle-level objects inside jets (PFCandidates / tracks / topoclusters), not just jet kinematics.
