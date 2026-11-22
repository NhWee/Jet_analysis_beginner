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
| CMS-001 | CMS | 2011A Jet primary dataset, Jet300, MOD HDF5 (`CMS2011AJets`, `cms`), read with `energyflow` 1.4.0. **Used: 3 of 18 files**, `CMS_Jet300_pT375-infGeV_{0,1,2}_compressed.h5` (md5-checked), stored in `raw/mod/datasets/CMS2011AJets/CMS_Jet300_pT375-infGeV/`. 284059 jets after pT > 375 GeV, \|η\| < 1.9 | [10.5281/zenodo.3340205](https://doi.org/10.5281/zenodo.3340205), CC BY 4.0 | reprocessed AOD (detector level) | yes (PFCs) | first real-data EEC, [`docs/mod_data_results.md`](../docs/mod_data_results.md) |
| CMS-002 | CMS (MC) | 2011A MOD SIM: Pythia 6.4.25 Z2 + CMS GEANT4 simulation, reco PFCs plus matched gen-level particles (`CMS2011AJets`, `sim`), read with `energyflow` 1.4.0. **Used: 27 of 355 files**: SIM170 1/1, SIM300 6/24, SIM470 6/73, SIM600 5/78, SIM800 3/79, SIM1000 2/40, SIM1400 2/40, SIM1800 2/20 (files `_0` … `_{n-1}`, md5-checked, 6.0 GB), in `raw/mod/datasets/CMS2011AJets/SIM*_Jet300_pT375-infGeV/`. 2587771 jets after the data selection, weighted n_eff 742139 | SIM170 [10.5281/zenodo.3341500](https://doi.org/10.5281/zenodo.3341500), SIM300 [.3341498](https://doi.org/10.5281/zenodo.3341498), SIM470 [.3341419](https://doi.org/10.5281/zenodo.3341419), SIM600 [.3364139](https://doi.org/10.5281/zenodo.3364139), SIM800 [.3341413](https://doi.org/10.5281/zenodo.3341413), SIM1000 [.3341502](https://doi.org/10.5281/zenodo.3341502), SIM1400 [.3341770](https://doi.org/10.5281/zenodo.3341770), SIM1800 [.3341772](https://doi.org/10.5281/zenodo.3341772); CC BY 4.0 | detector simulation + matched truth | yes (PFCs and gens) | detector check, [`docs/mod_sim_results.md`](../docs/mod_sim_results.md) |

**Reminder:** EEC requires jet *constituents*. Before registering a sample, confirm it stores
particle-level objects inside jets (PFCandidates / tracks / topoclusters), not just jet kinematics.
