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
| gen-pythia-01 | Pythia 8.3 | HardQCD:all, pp 13 TeV, pTHatMin = 100 GeV, 2000 events | n/a (generated, seed=12345) | truth | yes | `src/validate_eec_pythia.py` |
| gen-toy-01 | toy generator | `src/make_sample.py`, 200 events -> `raw/sample_events.root` | n/a (generated, seed=12345) | particle-level | yes | `notebooks/01_what_is_a_jet.ipynb` |
|  |  |  |  |  |  |  |

Generated samples are registered too: with the seed they are as reproducible as a DOI.

**Reminder:** EEC requires jet *constituents*. Before registering a sample, confirm it stores
particle-level objects inside jets (PFCandidates / tracks / topoclusters), not just jet kinematics.
