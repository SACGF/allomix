# Issue #8 — LOD alongside MAE, and depth × markers sensitivity curves

Tracking: https://github.com/SACGF/allomix/issues/8

## Goal

Two complementary changes:

1. **Formally define LOD** in the paper following CLSI EP17-A2 (which is what clinical labs validate against, and is what the AlloSeq / Devyser / ScisGo numbers cited in our Table 2 effectively report).
2. **Characterise LOD as a function of (panel size, sequencing depth, donor-host relatedness)** so a reader with a panel of *X* markers at *Y* mean depth can read their expected LOD off a curve. Two relatedness facets only: unrelated and full sibling (the best and worst clinically relevant cases).

The end product is a figure that clinicians can look at and immediately situate their own assay against. If it works, this figure is a strong candidate to be promoted into the main paper (current Fig 4, replacing the MAE-only relatedness panel).

## LOD operational definition

We follow **CLSI EP17-A2** (the standard cited in Pierson-Perry 2012 and used by Blouin 2024 for ScisGo and Qama 2026 for Devyser). Three quantities:

- **LoB (Limit of Blank).** The 95th percentile of estimated donor fraction across replicates of a *blank* sample (true donor fraction = 0). Operationally:
  - `LoB = quantile_95(est_frac | true_frac = 0)`  computed over N replicates of pure-host admixture at the (relatedness, depth, n_markers) cell.
  - Conceptually this is "how nonzero can the estimator look on a pure-host sample, by chance?"
- **LoD (Limit of Detection).** The lowest true donor fraction at which ≥95% of replicates return `est_frac > LoB`. Operationally:
  - For each candidate true fraction *f* > 0, compute detection rate = fraction of replicates with `est_frac > LoB`.
  - Fit a logistic / probit curve `P(detected | f)` vs `log10(f)` across the dilution series.
  - LoD = `f` where the fitted curve = 0.95.
- **LoQ (Limit of Quantitation).** Out of scope for this issue — would require defining an acceptable bias/CV target (e.g. ≤25% CV) and is best treated as a separate piece of work once we have wetlab data.

Two notes on choices:

- We use the **EP17 LoB + ≥95% detection** definition rather than the simpler alternative ("CI lower bound > 0 in ≥95% of replicates") because it matches what the commercial tools in Table 2 actually report. The CI-based definition is included as a secondary column in the raw CSV so we can compare if reviewers ask.
- `est_frac` is bounded ≥0 by the optimiser, so LoB is well-defined (not symmetric like an absorbance/signal LoB would be).

These definitions go into `paper/methods.md` as a new subsection ("Limit of detection") and are referenced from `paper/results.md` wherever LOD is reported.

## Sweep design

Grid:

| Axis | Values | N |
|------|--------|---|
| Relatedness | unrelated, sibling | 2 |
| Mean depth (×) | 100, 250, 500, 1000, 2000 | 5 |
| Panel size (markers) | 25, 50, 75, 100, 200, 400 | 6 |
| True donor fraction | 0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05 | 7 |
| Replicates per cell | 60 | 60 |

Total estimator runs: 2 × 5 × 6 × 7 × 60 = **25,200**.

Rationale:

- **Panel sizes** chosen so our lab's 76-SNP rhAmpSeq panel sits naturally on the curve (the 75 and 100 points bracket it). The 400-marker point covers AlloSeq HCT (~200 markers, but they're well-chosen high-MAF) and gives readers a sense of where the curve flattens.
- **Depths** cover the regime of interest. 100× is the panel cutoff QC warns about; 1000× is our empirical operating depth; 2000× is what one would push to for low-fraction detection.
- **True fractions** chosen to bracket the expected LOD region (≤2%). The blank (f=0) anchors LoB. The probit fit needs at least 4-5 fractions spanning ~10% detection to ~99% detection; this gives 6 non-zero points.
- **Replicates** at 60 are enough to estimate a 95% detection point with a logistic fit to acceptable SE (rough Wald SE for a 95% binomial is ~3% absolute at N=60, propagating to a usable LoD CI). Smaller pilot runs (N=10) of the existing relatedness script suggest the cell-level cost is small enough that 60 is comfortable.

Cost sanity check: per replicate ≈ generate genotypes + blend + classify + estimate ≈ ~100–300 ms (most of the cost is `estimate_single_donor_bb` doing the grid + profile-likelihood CI). 25k runs × 0.2 s ≈ 80 minutes single-threaded, well under an hour with `-j 8`.

Pilot plan: first run a 10-replicate version of the grid (~4k runs, 10–15 min) to sanity-check that the LoB/LoD curves look monotonic in depth and panel size before committing to the full 60-rep run.

## Reuse from existing code

The new script should reuse, not duplicate:

| Need | Reuse |
|------|-------|
| Generate synthetic host/donor pair at given relatedness | `simulate.generate_related_genotypes(n_markers, relatedness, rng, maf_range=(0.2, 0.5))` (`src/allomix/simulate.py:508`) — same function the existing relatedness script uses. |
| Write host/donor as VCF | `simulate.write_genotype_vcf` (`simulate.py:640`) |
| Generate chimeric sample at given depth/fraction | `simulate.blend_vcfs(...)` (`simulate.py:689`) with `realistic_biases=True, error_rate=0.01, locus_dropout_rate=0.016, depth_cv=0.43` — matches the existing relatedness and depth scripts' noise model. |
| Write chimeric VCF | `simulate.write_vcf` (`simulate.py:898`) |
| Classify markers | `genotype.classify_markers` |
| Estimate fraction + CI | `chimerism.estimate_single_donor_bb` (returns `donor_fraction`, `donor_fraction_ci`, `n_informative`) |

What's new (must be written):

- Outer sweep loop over the grid.
- Caching: host/donor genotypes depend only on (relatedness, n_markers, replicate seed), not on depth or fraction. Generate once per `(relatedness, n_markers, rep)`, reuse across all depth × fraction cells. Write the host/donor VCFs once per `(relatedness, n_markers, rep)` to a temp dir.
- LoB/LoD computation per cell: `compute_lob(blanks)`, `fit_lod(detection_rates_by_fraction)`.
- Logistic fit on `log10(f)`: 4-parameter logistic is overkill; standard probit / 2-parameter logistic (`P = 1 / (1 + exp(-(a + b·log10(f))))`) is sufficient. Solve for `log10(f95) = (logit(0.95) - a) / b`.
- Plotting code for the new figure.

No new src/ code needed in `allomix/` itself — this is all paper-scripts territory. If the LoB/LoD helpers turn out to be useful for the CLI's QC output later, we can promote them to `src/allomix/qc.py` then, not now.

## Output CSVs

Following the existing facts/ convention (one row per atomic observation; vibepaper interpolates from named summary CSVs).

### `output/facts/lod_grid_raw.csv` — per-replicate raw data

One row per `(relatedness, depth, n_markers, true_frac, replicate)`:

```
relatedness, depth, n_markers, true_frac, rep, seed, est_frac, ci_lo, ci_hi, n_informative
```

This is the durable artefact. If we want to change the figure (different curve colours, log/linear axes, add a panel) we re-render from this CSV — no recompute.

### `output/facts/lod_summary.csv` — per-cell LOD

One row per `(relatedness, depth, n_markers)`:

```
relatedness, depth, n_markers,
  lob_pct,           # 95th percentile of est_frac at f=0, in %
  lod_pct,           # logistic-fit 95% detection point, in %
  lod_pct_ci_lo,     # bootstrap CI on LoD (200 resamples of replicates)
  lod_pct_ci_hi,
  mean_informative,
  median_n_informative
```

Bootstrap CIs on LoD are nice-to-have (clinical reviewers like seeing uncertainty on a sensitivity claim). 200 resamples × cheap logistic fit is fast.

### `output/facts/lod_headline.csv` — paper-text-ready snapshot

One row, headline numbers for the Results paragraph and Table 2:

```
unrelated_lod_1000x_76markers_pct,
sibling_lod_1000x_76markers_pct,
unrelated_lod_500x_200markers_pct,
...
```

Whichever cells the Results prose actually quotes. Easier than indexing into the bigger CSV from vibepaper templates.

## Figure design

**Title (working):** "Limit of detection as a function of panel size and sequencing depth."

**Layout:** 1 row × 2 columns (facets). Left = unrelated, right = full sibling. Shared y-axis (log10 LoD %), shared x-axis (log10 panel size). One coloured curve per depth (5 curves: 100×, 250×, 500×, 1000×, 2000×). Markers at each data point, shaded band for bootstrap LoD CI.

**Annotations:** dashed horizontal lines at clinically meaningful thresholds (0.1%, 0.5%, 1%) labelled on the right of the right facet. Vertical reference line at 76 markers (our panel) on both facets so readers can immediately see "at my panel size, at my depth, I'd expect X% LoD."

**Axes:**
- x: panel size, log scale, ticks at the actual sweep values (25, 50, 75, 100, 200, 400)
- y: LoD (%), log scale, range maybe 0.05% to 10%
- legend: depth (×)

**Why two facets, not four:** the cousin / half-sibling cases sit between unrelated and sibling and clutter the figure. They stay in the existing Fig 4 (relatedness MAE) as the spread reference. Sibling is the clinically worst case; unrelated is what most patients are. Two facets is enough to make the point.

**Promotion question:** I think this *should* go in the main paper. Current Fig 4 has three panels: (a) informative markers vs relatedness, (b) MAE vs relatedness, (c) truth-vs-estimated scatter. Options:

- **Option A (preferred):** New main-paper figure (becomes Fig 5; current Fig 5 multidonor → Fig 6, etc.). Keep current Fig 4 intact. Two big additions to the paper, but each tells a distinct story.
- **Option B:** Replace current Fig 4 panel (b) (MAE vs relatedness) with the LoD curves. Saves a figure slot, but loses the head-to-head MAE bar that's currently the easiest read of "sibling vs unrelated."
- **Option C:** Demote the existing Fig 4 to supplementary, promote the new one to main.

Recommend A. The MAE-by-relatedness panel is the cleanest evidence that the estimator behaves well across the engrafted range even when informative markers are scarce; LoD curves are the clearest evidence that the *low-fraction* behaviour is competitive with commercial tools. Different audiences, both worth keeping.

## File changes

### New

- `paper/scripts/run_lod_validation.py` — the sweep. Outputs `output/facts/lod_grid_raw.csv` and `output/facts/lod_summary.csv`. Plotting kept in a separate script (below) so we can re-render the figure without re-computing.
- `paper/scripts/plot_lod_curves.py` — reads `lod_summary.csv` and produces `output/facts/fig5_lod_curves.png`. Cheap to re-run.
- `tests/test_lod_validation.py` — small unit test on the LoB / logistic-fit helpers. Synthetic input where the answer is known.

### Modified

- `paper/Snakefile` — add `lod_validation` rule (depends on nothing in `tests/test_data/`, generates `LOD_FACTS`) and `lod_plot` rule (depends on `lod_summary.csv`, outputs `fig5_lod_curves.png`). Wire into `ALL_FACTS` / `ALL_FIGS`. Split rules so the cheap plot can re-run without rerunning the sweep.
- `paper/methods.md` — add "Limit of detection" subsection with the EP17-A2 definitions above. Cite CLSI EP17-A2 and Pierson-Perry 2012 in `paper/references.bib`.
- `paper/results.md`:
  - Add a new subsection "Limit of detection" after the depth section, referencing the new figure and quoting LoD at our panel size / clinically representative cells.
  - Update Table 2 (`results.md:39`) to replace "~0.6% MAE (in silico)" in the **allomix** row with both MAE and LoD numbers, e.g. "0.6% MAE / 0.X% LoD (100 markers, 1000×, unrelated)". Keep the qualifier explicit because LoD depends on panel + depth.
  - The existing line about "limit of detection has not been formally characterized following AMP guidelines" (`discussion.md:38`) becomes obsolete — remove or rewrite.
- `paper/discussion.md` — short paragraph contextualising our LoD against vendor numbers (AlloSeq 0.3%, Devyser 0.05–0.1%, ScisGo 0.5%). Honest framing: the comparison is in silico vs analytical, so it sets the ceiling on what our estimator can achieve given perfect VCF input, not the floor of what a wetlab run will deliver.
- `paper/references.bib` — add CLSI EP17-A2 and Pierson-Perry 2012.

### Not changed

- `src/allomix/` — no source changes. Estimator and simulator already have everything needed.
- Existing scripts (`run_relatedness_validation.py`, `run_depth_validation.py`) — left alone. The new script does not subsume them; the MAE-by-relatedness story still belongs in the existing Fig 4.

## Snakemake integration

```python
LOD_FACTS = [
    f"{FACTS_DIR}/lod_grid_raw.csv",
    f"{FACTS_DIR}/lod_summary.csv",
    f"{FACTS_DIR}/lod_headline.csv",
]
LOD_FIGS = [f"{FACTS_DIR}/fig5_lod_curves.png"]

rule lod_validation:
    output:
        LOD_FACTS,
    shell:
        "python paper/scripts/run_lod_validation.py"

rule lod_plot:
    input:
        f"{FACTS_DIR}/lod_summary.csv",
    output:
        LOD_FIGS,
    shell:
        "python paper/scripts/plot_lod_curves.py"
```

Wire `LOD_FACTS + LOD_FIGS` into `ALL_FACTS` / `ALL_FIGS`. The sweep is the long-running rule; splitting plot out keeps figure-iteration fast.

Consideration: shard the sweep by relatedness for parallelism (two rules `lod_validation_unrelated`, `lod_validation_sibling` writing to per-relatedness CSVs, plus a merge rule). Not necessary for first cut — single-script ~80 min is fine.

## Open questions / things to confirm before implementing

1. **CI on LoD.** Bootstrap over replicates (200 resamples) is straightforward. Worth doing, or overkill for a first cut? Reviewers in clinical journals tend to like seeing it. Lean toward including.
2. **MAF range for synthetic markers.** Currently `(0.2, 0.5)`. This is favourable; real panels are typically `(0.1, 0.5)` with the tail dominating. Worth a sensitivity sub-analysis? I'd say: keep `(0.2, 0.5)` for the headline curves (matches the existing relatedness script for consistency), and document the assumption in methods. A supplementary sweep at `(0.1, 0.5)` could follow if reviewers ask.
3. **Profile-likelihood CI behaviour at very low f.** The estimator's CI lower bound is bounded at 0. At f=0.001 with 25 markers, we expect a lot of cells where the MLE is 0 and the CI is `(0, hi)`. This is fine for EP17-style detection (`est_frac > LoB`), but it means the "CI lower > 0" alternative definition will fail at very small panel × shallow depth. Worth flagging in the methods text so we don't get a reviewer comment about it.
4. **Sibling MAF assumption.** `generate_related_genotypes("sibling", ...)` uses Mendelian segregation given parental allele frequencies. For HLA-matched sibling donors specifically there's additional non-random IBD around the HLA locus — irrelevant for sample-ID panels which avoid HLA, but worth a methods sentence.
5. **Number of replicates.** 60 vs 100. 60 is fine for the headline LoD point; 100 tightens the bootstrap CI. ~30% more compute. I'd start at 60 and only re-run at 100 if the LoD CIs look too wide to make claims.
6. **Whether to also report LoB.** Yes — it's cheap, and reviewers will ask. One column in `lod_summary.csv` and one sentence in Results.

## Followup: align all validation figures on bias-corrected estimates

A pilot run of this sweep with bias correction *off* (matching the existing depth and relatedness scripts) produced LoB values of 1–3% across most cells, dominated by the simulator's heavy-tailed per-marker bias mixture. With perfect bias known the LoB collapses and the LoD becomes meaningfully better than vendor specs, which is the headline number this paper should be reporting.

Most published assay characterizations report best-case calibrated performance up front. The depth, relatedness, and multidonor validations currently default to no bias correction, with a separate "Effect of Per-Marker Bias Correction" subsection at one cell (2000×, 100 markers) to demonstrate the gain. This is a defensible structure but undersells the tool's achievable performance.

Followup work (separate to issue #8): rerun `run_depth_validation.py`, `run_relatedness_validation.py`, and `run_multidonor_validation.py` with bias correction on, so Tables 1, 3, 4 and Figures 1–4 all sit on the same calibrated baseline. The current bias-correction subsection would be reframed as the *uncalibrated* baseline for comparison. Numbers in those tables and figures will change; expect Table 1 MAE to drop, Table 3 sibling MAE to drop, and Table 4 per-donor MAE to drop.

For issue #8 itself, the LoD sweep is run with bias correction *on* (perfect per-marker biases provided to the estimator) so the LoD headline numbers reflect achievable analytical sensitivity given panel calibration. A short Methods sentence calls out the deliberate inconsistency with the rest of the paper until the followup is done.

## Validation that the implementation is working

Before committing to the full 25k-run sweep:

- **Sanity 1:** at `depth=2000, n_markers=400, unrelated`, LoD should be well below 0.5% (this is the easiest cell — lots of informative markers, deep coverage, max signal).
- **Sanity 2:** at `depth=100, n_markers=25, sibling`, LoD should be much worse, probably several percent (few informative markers, shallow depth, low signal).
- **Sanity 3:** LoD curves monotonic: increasing depth at fixed n_markers should not make LoD worse; increasing n_markers at fixed depth should not make LoD worse (within noise).
- **Sanity 4:** LoB should be small (well under 1%) across all cells. If LoB blows out at small panel sizes that's a real finding, but worth eyeballing before promoting the figure.
- **Sanity 5:** at the (depth=1000, n_markers=100, unrelated) cell, the LoD should land roughly in the 0.1–0.5% range — anything wildly different (e.g., 5%) suggests a bug in the LoB/LoD computation, not the estimator.

If the pilot run (10 reps) clears these sanity checks, scale to 60 reps.
