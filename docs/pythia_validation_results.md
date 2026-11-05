# Pythia validation — results

Step 1 of the README Quick start, run on Pythia 8 generator truth. This is the first result in
the repository that comes from a real parton shower and hadronisation model, not the toy
generator.

```bash
python src/validate_eec_pythia.py --n-events 10000 --seed 12345 --compare-hadronisation
```

Pythia 8.318 (PyPI `pythia8mc`), `HardQCD:all = on`, pp at 13 TeV, `PhaseSpace:pTHatMin = 100`
GeV. Final-state visible particles with |η| < 4 and pT > 0.2 GeV. Anti-kT R = 0.4, jet pT > 100
GeV, |η| < 2. Runtime about 2 minutes on one CPU core.

> **Generator truth, not data.** There is no detector, no pileup and no unfolding. These numbers
> show that the pipeline produces the expected physics on a known input. They are not a
> measurement.

---

## 1. Pipeline check

| quantity | result |
|---|---|
| events | 10000 |
| jets | 9129 |
| ⟨pT⟩ | 139.5 GeV |
| ⟨n_const⟩ | 26.5 |
| EEC normalisation | **1.000000** ✅ |

**Verdict: PASS.** Step 1 of the Quick start has now actually been passed.

---

## 2. EEC in jet-pT slices

![EEC in pT slices](../results/figures/eec_pythia_pt_slices.png)

| jet pT [GeV] | jets | EEC peak R_L | peak R_L·⟨pT⟩ [GeV] | median rel. stat. error |
|---|---|---|---|---|
| 100–150 | 6788 | 0.102 | 12.3 | 2.8% |
| 150–200 | 1610 | 0.059 | 10.0 | 5.2% |
| 200–300 | 631 | 0.051 | 11.9 | 7.4% |
| 300–500 | 91 | 0.026 | 9.2 | 15.4% |

Errors are the standard error of the mean over jets (see `jeteec/eec.py`).

**The small-angle side scales with 1/pT.** In R_L (left panel) the curves separate: higher-pT
jets reach further into small angles. In R_L·⟨pT⟩ (right panel), the rising small-angle side of
all four slices falls on one curve up to about 3 GeV. The peak moves by a factor of 4 in R_L but
stays at 9–12 GeV in R_L·pT. So the transition is set by a fixed momentum scale and not by an
angle, which is what a hadronisation (confinement) transition should look like. The large-angle
edge, by contrast, is fixed at R_L ≈ R = 0.4 in the left panel and spreads out in the right one,
because it comes from the jet radius rather than from QCD dynamics.

**Caveats.**
- The peak position is read off from 50 log bins over three decades, so it is quantised in steps
  of about 15%.
- The 300–500 GeV slice has 91 jets. Use more events or a higher `pTHatMin` before quoting it.
- The absolute R_L·pT value of the peak depends on how "peak" is defined and binned. Only its
  constancy across slices is a robust statement here; don't compare 9–12 GeV directly with a
  number in the literature.

---

## 3. Hadronisation on vs off

![Hadronisation comparison](../results/figures/eec_pythia_hadronisation.png)

Same seed, same cuts, `HadronLevel:all = off`. Jets are clustered from the post-shower partons.

| | hadron level | parton level |
|---|---|---|
| jets | 9129 | 9306 |
| ⟨n_const⟩ | 26.5 | 7.9 |
| EEC peak R_L | 0.102 | 0.269 |
| binned EEC area | 0.876 | 0.621 |

| shape ratio parton/hadron (unit area in range) | |
|---|---|
| 0.1 < R_L < 0.4 (perturbative) | 1.19 |
| R_L < 0.01 (hadronisation) | **0.28** |

**The small-angle part of the EEC comes from hadronisation.** At large angles the parton-level
and hadron-level shapes agree to about 20%. Below R_L = 0.01 the parton level carries only about
a quarter of the hadron-level weight. Hadronisation moves energy flow into small angles.

**A pitfall found while doing this.** The first comparison showed the parton level about 30%
*below* hadron level everywhere, including the perturbative region. That is not physics. Every
EEC sums to 1 *including* the self-pairs at ΔR = 0, which fall outside the bins. A parton-level
jet has about 8 hard constituents, so its self-pair weight Σzᵢ² is much larger and less weight is
left for the visible bins (binned area 0.62 vs 0.88). The comparison must use shapes normalised
to a common range. **The same trap applies to any EEC comparison between samples with different
constituent multiplicity**: charged-only vs all particles, different pT thresholds, or data vs MC
at different tiers.
