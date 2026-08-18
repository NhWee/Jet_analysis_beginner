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
| _example_ | Pythia8 | dijet pp 13 TeV, pThat>100 | n/a (generated, seed=12345) | truth | yes | EEC validation |
|  |  |  |  |  |  |  |

**Reminder:** EEC requires jet *constituents*. Before registering a sample, confirm it stores
particle-level objects inside jets (PFCandidates / tracks / topoclusters), not just jet kinematics.
