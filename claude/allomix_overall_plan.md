# Donor Chimerism Tool — Overall Build Plan

This is the master plan for building a general-purpose, panel-agnostic NGS chimerism monitoring tool. Each step is designed to be fed into a prompt as a self-contained task.

---

## Step 1: Project Setup and README ✅ COMPLETE

- Project name: **allomix** (registered on PyPI, v0.0.1 published)
- Directory structure: `src/allomix/`, `tests/`, `docs/`, `scripts/`, `data/`
- `pyproject.toml` with dependencies (cyvcf2, numpy, scipy), dev deps (pytest, ruff), CLI entry point
- `README.md` — general-purpose tool description, workflow, input/output contract, comparison table
- `CLAUDE.md` — project context + dev conventions
- CLI stub with `monitor` and `timeline` subcommands

---

## Step 2: BAM vs VCF — Determine Primary Input Format ✅ COMPLETE

**Decision: VCF as primary input. No BAM support in v1.**

Full analysis in `claude/step2_bam_vs_vcf_decision.md`. Key points:

- VCF AD fields provide everything needed (ref count, alt count, depth)
- Joint calling ensures hom-ref samples at variant sites still have 2-element AD — confirmed with `data/joint_called_example.vcf`
- Pipeline is flexible and can be adjusted to produce what allomix needs
- Minimum required FORMAT fields: GT, AD, DP
- Includes an audit script for running on /tau to verify existing VCFs

### Example Data Available

- `data/idt_rhampseq_sid_example.vcf` — single sample, 15 markers (de-identified coordinates)
- `data/joint_called_example.vcf` — 114 samples, 9 markers (de-identified coordinates) confirms joint-calling provides ref+alt AD at all sites

---

## Step 3: Test Data Generation — Synthetic Chimeric Files ✅ COMPLETE

Built in `src/allomix/simulate.py` + `scripts/generate_test_data.py` + `tests/test_simulate.py` (64 tests).

- VCF blending: takes host + donor VCFs, mixture fraction, target depth → synthetic chimeric VCF with binomial-sampled allele counts
- **Capture bias simulation**: `marker_bias_sd` parameter adds per-marker capture/amplification bias drawn from N(0, sd). 0.0 = ideal, 0.02 = realistic for capture panels. Each marker gets a fixed bias that shifts its observed VAF relative to truth.
- Plain-text VCF parsing (no cyvcf2 dependency) so simulation code is lightweight
- CLI scripts: `generate_test_data.py` (chimerism series), `generate_timeline_data.py` (engraftment + relapse scenario)
- Supports custom fractions, depth, random seed, and `--bias-sd`

---

## Step 4: Clone and Examine Reference Open-Source Projects ✅ COMPLETE

Full analysis in `claude/step4_reference_tool_analysis.md`. (now claude/historical/step4_reference_tool_analysis.md) 9 repos examined.

**Key findings:**

| Tool | Key Takeaway |
|------|-------------|
| **Demixtify** | MLE likelihood framework (Crysup & Woerner 2022 Formula 5) is our starting point. Known genotypes simplify dramatically. AGPL — reimplement math, cite paper. |
| **Chimerism-Bias** | Per-marker amplification bias correction is essential; ~30% error reduction. Port the correction formulas. |
| **Chimerism-FABCASE/nMarkers** | 76 markers is more than sufficient even for sibling pairs. MIT. |
| **All-FIT** | Grid search + outlier removal pattern. MIT. |
| **EuroForMix** | Multi-contributor architecture reference. LGPL. |
| **Conpair** | Brent refinement after grid search; log-space arithmetic. Non-commercial. |
| **somalier** | Extract-then-analyze architecture; site selection. MIT. |

**License approach:** Implement math independently (published science, cite Crysup & Woerner 2022), do not copy AGPL/non-commercial code.

---

## Step 5: Detailed Implementation Plan ✅ COMPLETE

Full plan in `claude/step5_implementation_plan.md`. Defines:

- 6 modules: genotype, chimerism, bias, qc, report, cli
- MLE algorithm: Demixtify Formula 5 with known-genotype simplification
- Grid search (1001 points) + Brent refinement + profile likelihood CI
- 9-phase implementation order
- Full test plan

---

## Step 6: Implement Core Algorithm ✅ COMPLETE (single-donor)

### Implemented modules:

| Module | File | Purpose |
|--------|------|---------|
| genotype | `src/allomix/genotype.py` | VCF parsing (cyvcf2), marker joining by (chrom,pos,ref,alt), Vynck type classification (6 types), depth/GQ/PASS filtering |
| chimerism | `src/allomix/chimerism.py` | MLE estimation: per-marker log-likelihood (Crysup & Woerner Formula 5 with known genotypes), 1001-point grid search, Brent refinement, profile likelihood 95% CI, per-marker residuals, 3-SD outlier flagging |
| qc | `src/allomix/qc.py` | QC assessment: marker counts, depth stats, GOF chi-squared, CI width, pass/fail with warnings |
| report | `src/allomix/report.py` | TSV output (summary + verbose per-marker detail), JSON, timeline format |
| cli | `src/allomix/cli.py` | `allomix monitor` and `allomix timeline` wired end-to-end with all options |

### Test coverage: 261 tests passing

- 29 genotype tests (parsing real VCFs, classification, filtering)
- 55 chimerism tests (MLE math, estimation accuracy at multiple fractions, CI coverage, edge cases)
- 48 multi-donor tests (unit + integration + CLI)
- 21 bias tests
- 12 QC tests (pass/fail conditions, warnings)
- 17 report tests (TSV/JSON format, timeline)
- 64 simulate tests (blending logic, round-trips)
- 15 integration tests (full pipeline: synthetic VCF → genotype → chimerism → qc → report → CLI)

### What works now:

```bash
# Single-donor chimerism from VCFs
allomix monitor --host host.vcf --donor donor.vcf --sample admix.vcf -o results.tsv

# JSON output
allomix monitor --host host.vcf --donor donor.vcf --sample admix.vcf --format json

# Timeline across timepoints
allomix timeline --host host.vcf --donor donor.vcf --sample d30.vcf --sample d60.vcf -o timeline.json

# Verbose per-marker detail
allomix monitor --host host.vcf --donor donor.vcf --sample admix.vcf --verbose

# Estimate per-marker bias from training samples
allomix estimate-bias --vcfs *.vcf.gz -o bias_table.tsv

# Monitor with bias correction
allomix monitor --host host.vcf --donor donor.vcf --sample admix.vcf --bias-table bias.tsv
```

---

## Step 7: Implement Multi-Donor Support ✅ COMPLETE

**Goal:** Extend chimerism estimation to host + 2 donors.

### Implementation (complete, 2026-04-08)

Detailed plan: `claude/multi_donor_plan.md`

**simulate.py** — `generate_sibling_trio_genotypes()` (Mendelian segregation from shared parents), `_mendelian_child()`, `expected_vaf_multi()`, `blend_from_genotype_dicts()`

**genotype.py** — `InformativeMarker` gained `marker_types: list[int | None]` and `informative_for: list[bool]` fields. `classify_markers()` now includes markers informative for ANY donor (was: first donor only).

**chimerism.py** — `MultiDonorResult` dataclass, `expected_weight_multi()`, `total_log_likelihood_multi()`, `estimate_multi_donor()` (triangular grid search at 101 steps → Nelder-Mead → profile likelihood CIs with chi2 df=1 per donor), `_per_marker_results_multi()`.

**qc.py** — `QCReport.per_donor_n_informative`, per-donor CI width and informativity warnings.

**report.py** — `_write_tsv_multi()`, multi-donor branches in `to_json()` and `timeline_json()`.

**cli.py** — `_run_single_sample()` auto-detects: 1 donor → `estimate_single_donor()`, 2+ → `estimate_multi_donor()`.

**Test data** — `tests/test_data/multidonor/`: 3-brother sibling VCFs (100 markers, 61 informative for any donor, 46/41 per donor) + 22 chimeric VCFs at a grid of (f1, f2) points. Generated by `scripts/generate_multidonor_test_data.py`.

**Tests** — `tests/test_multidonor.py`: 48 tests (unit + integration + CLI). 261 total tests pass, zero regressions.

**Validation** — Estimation accuracy on sibling donors at 1000x: pure host <1%, balanced 25/25 → 24.3/26.0%, asymmetric 30/10 correctly distinguished, pure donor1 >98%.

### Paper updates ✅ COMPLETE

All paper sections updated: methods (multi-donor extension subsection), results (sibling donor validation), discussion (moved from limitation to capability), abstract (multi-donor mention), README. Validation script (`paper/scripts/run_multidonor_validation.py`) and figure (`paper/scripts/generate_multidonor_figure.py`, `paper/figures/fig_multidonor.png`) complete.

---

## Step 8: Implement Bias Correction ✅ COMPLETE

Implemented in `src/allomix/bias.py`.

- `estimate_bias()` — estimate bias from training VCFs: bias = median(VAF_het - 0.5) per marker
- `load_bias_table()` / `save_bias_table()` — TSV I/O for bias tables
- Bias correction integrated into MLE: adjusts expected reference allele weight per marker
- CLI: `allomix estimate-bias --vcfs *.vcf.gz -o bias_table.tsv`
- CLI: `allomix monitor --bias-table bias.tsv ...`
- Validation shows bias correction reduces MAE ~15% and max error ~25% at 2000x depth with realistic biases

---

## Step 9: In-Silico Validation ✅ COMPLETE

### Simulation framework

The simulator (`src/allomix/simulate.py`) models four sources of measurement noise, all calibrated from empirical panel characterisation (210 joint-called VCFs, 18,047 samples, 76-SNP rhAmpSeq panel — results in `paper/empirical_results/`):

1. **Per-marker amplification bias** — heavy-tailed Gaussian mixture: 95% from N(0, 0.012), 5% from N(0, 0.08), yielding overall SD ~0.018 matching the empirical distribution (median |bias| 0.005, 95th pct 0.041, max 0.10)
2. **Non-uniform depth across markers** — per-marker depths drawn from log-normal distribution matching empirical CV=0.43 (mean 1,732x, range 285–2,789x)
3. **Sequencing errors** — symmetric error model at ε=0.01
4. **Locus dropout** — 1.6% per-marker dropout rate matching empirical no-call rate

### Validation scripts (in `paper/scripts/`)

| Script | What it tests |
|--------|--------------|
| `run_depth_validation.py` | Accuracy across 5 depths (50x–1000x), 12 donor fractions |
| `run_relatedness_validation.py` | Accuracy across 4 relatedness levels (unrelated to sibling), 10 replicates each |
| `compare_bias_correction.py` | Side-by-side comparison with/without bias correction |
| `generate_paper_facts.py` | Generates all facts CSVs for the paper |

### Key results

**Depth validation** (100 markers, 80 informative, realistic noise):

| Depth | MAE (%) | RMSE (%) | Max Error (%) |
|:---:|:---:|:---:|:---:|
| 50x | 0.97 | 1.09 | 2.04 |
| 100x | 0.69 | 0.79 | 1.30 |
| 200x | 0.66 | 0.73 | 1.10 |
| 500x | 0.62 | 0.66 | 1.03 |
| 1,000x | 0.61 | 0.68 | 0.96 |

**Relatedness validation** (100 markers, 500x, 10 replicates):

| Relatedness | Mean Informative | MAE (%) |
|:---|:---:|:---:|
| Unrelated | 58 | 0.99 |
| 1st cousin | 55 | 1.04 |
| Half-sibling | 52 | 0.98 |
| Full sibling | 38 | 1.16 |

All MAE values sub-2% (clinically acceptable). Even sibling donors maintain sufficient informative markers (min 27, well above the minimum 3 required).

**CI coverage** is 25–58% (below nominal 95%). This is expected: the binomial likelihood does not model the systematic biases and non-uniform depth. The paper discusses this and notes approaches for improvement (bias correction, beta-binomial likelihoods, empirical recalibration).

---

## Step 10: VariantGrid Integration 🔲 TODO

- JSON output schema agreed with VG team
- VG stores donor/host genotypes, exports as VCF
- VG ingests allomix JSON results per patient
- VG renders timeline chart from chimerism results across timepoints
- Exact API integration TBD

---

## Step 11: Real Sample Validation ✅ COMPLETE

- Joint-called VCFs for the idt_rhampseq_sid panel produced on /tau; joint-called VCF available locally at `output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz`
- Batch runner `scripts/run_xls_batch.py` drives allomix across the patient list in `output/Chimerism project patient list.xlsx`, using the bias table from `output/bias_training/bias_table.tsv` and appending clinical reference columns (`Donor`, `Chimerism result TP2`) to `batch.tsv`
- Validation runs captured in `output/validation_run_new_bias2/` (post-BB / error-adjusted GoF fix). All 7 PASS samples now produce non-trivial `gof_pval` (previously all 0.0000).
- Concordance assessment vs. clinical sorted-cell chimerism (CD45 / CD3 / CD13) is not a direct apples-to-apples comparison: allomix reports bulk DNA chimerism, which is a cell-type-weighted average and tracks CD13 myeloid more closely than CD45 in samples with strong lineage disparity (e.g. 20_MO RCAR: allomix 40.79%, CD45 46.78%, CD3 93.19%, CD13 30.58%).
- **Re-run pending against the new two-phase pipeline (2026-05-29).** The validation above used the all-GATK joint-calling pipeline, which we have since discovered strips minority ALT reads from `FORMAT/AD` at hom-ref calls (confirmed empirically: 0 ALT reads across ~9M reads at hom-ref calls in joint-called VCFs). The new `pipeline/Snakefile` runs GATK only on HOST/DONOR and forced `bcftools mpileup` on ADMIX, preserving raw AD. Real-data results should be regenerated against this pipeline before any downstream analysis is treated as final. Blocks the paper methods/discussion rewrite in Step 18.
- Next phase (out of this step): controlled dilution series for quantitative accuracy validation.

---

## Step 12: Per-Marker Likelihood Context Refactor 🔲 TODO

Pure-refactor pre-step for Steps 14, 15, 16, and 17. Each of those steps independently
proposes adding a new optional kwarg to `total_log_likelihood_bb`,
`total_log_likelihood_multi_bb`, `estimate_single_donor_bb`, `estimate_multi_donor`,
and `_profile_likelihood_cis_multi`, then threading it through every nested closure.
By the time all four land, those signatures will have grown four new optional
parameters each, plus the existing `marker_biases`, and every call site has to
forward all of them.

This step does that plumbing once. Introduce a `PerMarkerContext` dataclass that
the estimators build once per call from whatever inputs are configured (today:
`error_rate` + optional `marker_biases`); the aggregators take a single
`ctx: list[PerMarkerContext]` aligned with the markers list. Each downstream step
then mutates one field on the context rather than threading another kwarg through
every closure.

No CLI changes, no behavioural changes. Success criterion is the existing 261-test
suite passing unchanged plus a numerical-regression spot-check on the multi-donor
fixture and the April-24 validation batch.

Detailed plan: `claude/12_marker_context_refactor_plan.md`.

---

## Step 13: Beta-Binomial Goodness-of-Fit ✅ COMPLETE

The MLE already uses a beta-binomial likelihood (fits both `f` and overdispersion `ρ`), but the gof chi-squared in `qc.py` standardised residuals by binomial variance, so `gof_pval` was ~0 on every real-data sample even when the fit was fine. Fix landed:

- `ρ` now plumbed through `ChimerismResult` and `MultiDonorResult` (default `float("inf")` for backward compatibility with old fixtures).
- `_compute_gof_pval` in `qc.py` uses beta-binomial variance `p(1-p)(n+ρ)/(n(ρ+1))`, df corrected to `n_markers - n_fitted_params` (single-donor: 2, multi-donor: k+1).
- Follow-up fix discovered during validation: at f near 0 or 1, the raw `expected_vaf` saturates at 0 or 1, making the variance floor (clamped at `1-1e-6`) collapse. A typical ~1% sequencing-error residual against that tiny floor produced spurious chi-sq blow-ups for 100%-donor samples. Introduced `_error_adjusted_p_alt()` in `qc.py` using the same 4-state error model as `log_likelihood_marker_bb`, driven by `result.error_rate`. The variance floor now reflects the actual error rate at saturated markers.
- Real-data result: all 7 PASS samples in the idt_haem validation batch now produce sensible `gof_pval` (0.46–1.00 range), vs all 0.0000 previously.

Detailed plan: `claude/beta_binomial_plan.md`.

---

## Step 14: Empirical Per-Site Error Rates 🔲 TODO

Replace the global `--error-rate 0.01` constant with empirically measured rates from the bias-training cohort. At hom-ref sites measure observed ALT-read rate, at hom-alt sites measure observed REF-read rate. Per-site (preferred) or per-sample fallback. Removes a chunk of modelling slack and reduces reliance on a hand-tuned constant. Output a per-site error rate table from `estimate-bias`-style tooling so it can be loaded alongside the bias table.

Detailed plan: [`claude/14_empirical_error_rates_plan.md`](14_empirical_error_rates_plan.md).

---

## Step 15: Per-Site Dropout Rate 🔲 TODO

The bias-training cohort already gives us per-site no-call rates. Integrate a per-site dropout probability into the likelihood so flaky sites are automatically downweighted rather than treated as fully informative when they happen to call. Estimate alongside the bias and error-rate tables; load as an optional input to `monitor`.

Detailed plan: [`claude/per_site_dropout_plan.md`](per_site_dropout_plan.md).

---

## Step 16: GQ-Weighted Marker Contributions 🔲 TODO

Currently `--min-gq 20` is a hard pass/fail. Replace with a per-marker weight on the likelihood contribution so borderline-confidence genotypes (e.g. GQ 20-30) stay informative but downweighted. Likely small gains relative to Steps 14 and 15, but cheap to add once the per-marker likelihood is already being modified.

Detailed plan: [`claude/16_gq_weighted_markers_plan.md`](16_gq_weighted_markers_plan.md).

---

## Step 17: Per-Base Quality-Aware Likelihood 🔲 TODO (skeptical, may not ship)

**Status note:** we are not convinced this is worth doing. It is the most invasive
of the remaining algorithm steps (new upstream `bcftools mpileup -a FORMAT/QS`
pipeline rule, new `bq.py` module, simulator extension to emit per-read BQs, new
`MarkerResult` field, new CLI flag, real-data revalidation) for what is likely a
small accuracy gain on the current panel where most reads sit at Q30+. Step 14
(empirical per-site error rates) attacks the same source of slack — over-reliance
on the global `--error-rate` constant — at a fraction of the engineering cost,
and the mean-phred → mean-error approximation in the BQ plan (Jensen's inequality;
see plan §"Mean-phred vs mean-error approximation") means BQ-aware is itself only
a partial fix. Decision gate: revisit only if Steps 14, 15, 16 land and CIs on
real data are still wider than we want.

The current MLE uses a flat sequencing error rate (default 1%) for all reads at all
markers. A more accurate model would weight each read's contribution to the likelihood
by its base quality (BQ/QUAL), following the approach used by Conpair (Bergmann et al.):

- Parse per-read base qualities from the VCF (if available) or BAM pileup
- Replace the flat error rate `e` in the likelihood with a per-read error probability
  derived from the Phred quality score: `e_i = 10^(-Q_i/10)`
- Per-marker likelihood becomes a product over individual reads rather than a binomial
  with aggregate counts
- This gives more weight to high-quality reads and down-weights low-quality bases,
  improving accuracy at low depths and near the detection limit
- Requires either BQ annotation in VCF FORMAT fields, or falling back to BAM access
- Profile likelihood CIs should improve as the model better captures per-read uncertainty

Implementation notes:
- Add optional `--bq-aware` flag to CLI (default off for backwards compatibility)
- If BQ data unavailable, fall back to current flat error rate
- Benchmark accuracy improvement vs computational cost on real data

Sequenced last because it requires an upstream pipeline change (`bcftools mpileup -a
FORMAT/QS` + `bcftools annotate`) and partially overlaps Step 14 in motivation
(both reduce dependence on the global `--error-rate` constant). Defer until 14, 15,
16 are in and we have seen whether they close the gap on real-data CIs.

Detailed plan: `claude/17_bq_aware_plan.md`.

---

## Step 18: Publication 🟡 IN PROGRESS

- Method paper describing the approach — framework set up with vibepaper
- Target journal: Journal of Molecular Diagnostics (Technical Advance)
- Paper sections in `paper/`, analysis scripts in `paper/scripts/`
- Cite: Crysup & Woerner 2022 (Demixtify MLE framework), Vynck et al. (bias correction)
- In silico validation complete (depth series, relatedness, bias correction, multi-donor)
- Multi-donor paper updates complete (methods, results, discussion, abstract, figures)
- Simulation calibrated from empirical panel characterisation (210 VCFs, 18,047 samples)
- Supplementary figures (S1-S6) for simulation model validation complete (`paper/scripts/generate_supp_synthetic.py`, Snakefile rules, supplementary text)
- Validation with real samples (Step 11) still needed
- Open-source tool release (MIT license, PyPI)

### Remaining paper tasks

- [ ] Add bias stability figure (`fig_bias_stability.png`) to `paper/results.md` near the "Effect of Per-Marker Bias Correction" section. This validates the fixed-bias-per-marker assumption (r = correlation between |median_bias| and within-marker SD). Caption template in `supp_synthetic.csv` facts.
- [ ] Decide: should the ablation study (Figure S4) also include a "no overdispersion" baseline (standard binomial vs beta-binomial)?
- [ ] **Rewrite joint-calling references for the two-phase pipeline.** Specifically: `methods.md:11` (the "joint calling preserves admix AD" claim is wrong — GATK HaplotypeCaller -ERC GVCF strips minority ALT reads at hom-ref blocks; new pipeline uses GATK for HOST/DONOR + `bcftools mpileup` for ADMIX) and `discussion.md:29` (same claim restated). `supplementary.md:5` is fine as-is (bias estimation uses het sites, which weren't affected). Gated on the Step 11 re-run landing so methods text and real-data results can be updated in one pass.

---

## Step 19: Intronic Shoulder Marker Evaluation 🔲 TODO

Our Haem capture panel's depth extends past the exon boundaries into the flanking introns (reads sequencing off the ends of captured fragments), forming a declining-depth shoulder. These intronic positions are a potential source of extra informative markers. Introns are under weaker purifying selection than exons, so their site frequency spectrum is shifted toward common (near-0.5 MAF) variants, which are exactly the high-heterozygosity markers most likely to distinguish host from donor. More informative markers means a tighter chimerism estimate and lower LOD.

The open question is whether they carry allele-specific bias that would distort VAF (and therefore the chimerism estimate):

- **Capture (hybridization) bias: expected to be negligible.** The intronic SNP sits outside the probe footprint (the probe is over the exon), so the polymorphic base does not affect duplex stability and both haplotypes are pulled down equally. The only capture effect is the depth drop, which is allele-symmetric and already handled by the beta-binomial weighting plus depth/dropout QC.
- **Read-end mapping bias: the real risk to check.** By construction these SNPs sit near read ends (reads sequence from the exon out into the intron). Alt reads carry a mismatch and are more likely to be soft-clipped or MAPQ-penalised there, dropping alt reads from AD and skewing observed VAF toward reference. This is allele-asymmetric and is NOT caught by the depth filter, because depth is allele-blind: a marker can have healthy depth and still carry a quiet reference skew.

**Analysis to run** (summary-stats script against /tau, no coordinates or patient IDs in output): for intronic-shoulder markers vs exon-core markers, report per-marker median het VAF (should center on 0.5), depth distribution, and dropout rate, binned by distance into the intron (intron offset).

- If het VAF only departs from 0.5 once depth has already fallen below the QC threshold, the depth filter alone suffices and the introns can be harvested directly.
- If het VAF drifts off 0.5 while depth is still healthy, the read-end mapping bias is real in the intermediate band and an explicit allele-balance filter (orthogonal to depth) is needed. Where bias is moderate and stable it folds into the existing per-marker bias term (`bias.py`); where large, filter the marker out.

Optionally rank candidate intronic markers by population MAF / expected heterozygosity (gnomAD AF) to prioritise the most informative ones.

Came out of a design discussion on 2026-05-27.

---

## Step 20: Host-Presence Detection at Donor-Homozygous Markers 🔲 TODO

A dedicated detection test for "is the host present at all?", separate from the fraction MLE, aimed at low-level host re-occurrence (relapse) post-HSCT. It uses only the markers where the donor is homozygous and the host carries the donor-absent allele: there the donor-absent allele sits at the sequencing-error background in a pure-donor sample, so its read counts give a one-sided count test (and an LRT yielding a host-fraction estimate) against that background, combined across markers. Same reads the MLE already sees, but reframed as detection and freed from the single shared overdispersion `ρ` and the global error rate, both of which blunt the MLE at very low fractions. Reported alongside the MLE (route A); a unified two-component likelihood is a follow-up (route B).

Depends on Step 14 (empirical per-site error rates) for the per-site background that sets the achievable detection limit; without it the test falls back to the global `--error-rate` and the limit is error-floor-bound. Soft dependency on Step 12 (per-marker context refactor) for route B only.

Build the validation controls first: we do not yet have extremely-low-fraction synthetic data (the issue #8 LoD sweep stops at 0.1%). The plan adds a control-generation + calibration step before any CLI work, generating low-fraction positive controls and error-only negative controls (EP17 LoB/LoD applied to the detection statistic) and checking the test is calibrated. Quality scores (Step 17) are not needed: the detector uses AD counts plus per-site error rates, and the controls declare their error rate rather than deriving it from per-read qualities.

Came out of a design discussion on 2026-05-28.

Detailed plan: [`claude/20_host_presence_detection_plan.md`](20_host_presence_detection_plan.md).

## Step 21: Calibrate Simulator Overdispersion for Realistic LoD 🟡 IN PROGRESS

Discovered 2026-05-28 while reconciling the in-silico LoD against real run3 patient LoDs (~0.5–1%) vs the paper's headline in-silico LoD (~0.13–0.32%). The simulator drew reads from a pure binomial, so the in-silico LoD reflects near-binomial sampling and is optimistic by ~3–5x. The per-marker beta-binomial variance approaches `p(1-p)/(ρ+1)` as depth grows, so the effective depth caps near `ρ+1` reads and the LoD saturates; overdispersion, not depth, is the dominant LoD control at clinical coverage.

Done:
- `simulate.sample_allele_counts` / `blend_vcfs` / `blend_from_genotype_dicts` now take a `rho` arg (default `inf` = binomial, unchanged). Tests in `tests/test_simulate.py`.
- New paper artefacts: `paper/scripts/plot_lod_saturation.py` (LoD vs depth, reconciles sim vs real) and `paper/scripts/run_overdispersion_lod.py` (LoD vs ρ). Wired into `paper/Snakefile` (rules `lod_saturation`, `overdispersion_lod`); figures added as Supplementary S7/S8 with `overdispersion_lod_headline.csv` facts; discussion + methods updated.
- `scripts/diagnose_sample.py` prints each real sample's fitted `rho` (per-sample, authoritative).

TODO:
- [ ] Calibrate `rho` from real per-sample fits (`diagnose_sample.py` on run3 VCFs), then re-run `lod_validation` with that `rho` so the **headline** LoD reflects real overdispersion rather than the binomial best case. This re-runs the expensive `lod_validation` job (warn before triggering).
- [ ] The simulator applies a single global `rho` to every marker/allele uniformly, including the near-zero donor-absent allele where overdispersion is not physical (it is a het/intermediate-marker amplification phenomenon). A marker-type-aware (or allele-aware) overdispersion model is needed before `rho` can be used to validate host-presence detection (Step 20) — otherwise turning on a global `rho` miscalibrates the presence-test null. See the note added to `claude/20_host_presence_detection_plan.md`.
- [ ] Decide whether the headline-LoD wording in `discussion.md` should switch from the binomial number to the overdispersion-calibrated number once the above lands.

---

## Step 22: Decide Whether to Fully Switch to the Pileup / Two-VCF Model ✅ COMPLETE

Committed to pileup-only on 2026-05-29. Verification gate ran against `output/validation_run6/batch.tsv` (wide-BED discovery + 71 force-called SID panel sites, additive Snakefile semantics); n_informative and host-presence magnitudes matched the prior run4 numbers, p-values preserved. Migration landed: `--vcf` removed from `monitor` / `timeline` (still present on `estimate-bias` / `estimate-errors`), tests/test_integration.py + tests/test_multidonor.py CLI tests rewritten as panel/admix pairs (joint VCF passed twice for synthetic fixtures), `_resolve_vcf_inputs` deleted, README/CLAUDE.md/doc/joint_calling.md updated. Full test suite: 320 pass.



Pipeline and CLI now both support the two-VCF model (panel VCF for host/donor `GT` from GATK; admix VCF for `AD` from forced `bcftools mpileup`). The original single-joint-VCF model is still supported as a back-compat path: `allomix monitor --vcf <single>`, synthetic test data (`tests/test_data/`, `scripts/generate_*.py`), the in-silico validation harness, and the paper figures all still live on the single-VCF model.

This step is the decision and the follow-through: do we keep dual support, or commit to pileup-only and delete the single-VCF branches?

Arguments for fully switching:
- The two-phase pipeline is now the only sane way to produce admix `AD` (see the empirical 0-ALT-reads-at-hom-ref result in Step 11 and `doc/joint_calling.md`).
- Dual support is maintenance drag: every new CLI arg, every estimator change, every paper figure has to consider both modes.
- Synthetic data generated under the single-VCF assumption can hide AD-stripping bugs that real data exposes.

Arguments against:
- The single-VCF mode is exactly how the in-silico simulator currently emits data — synthetic ground truth is convenient because we control both the panel GTs and the admix AD in one file.
- Existing paper validation runs (Steps 8, 9, 13) were done against single-VCF synthetic data; tearing this out invalidates a chunk of reproducibility.
- All 21 integration tests in `tests/test_integration.py` currently drive the single-VCF path. A switch means rebuilding the fixtures as two-file panel/admix pairs.

If we commit to pileup-only, the work to do:
- [ ] Drop `--vcf` from `monitor` / `timeline`; make `--panel-vcf` + `--admix-vcf` the only mode. Update `_resolve_vcf_inputs` accordingly.
- [ ] Update synthetic data generation (`scripts/generate_test_data.py`, `scripts/generate_timeline_data.py`, `scripts/generate_multidonor_test_data.py`, `src/allomix/simulate.py`) to emit a panel VCF + admix VCF pair rather than one joint VCF. Decide whether the simulator should model `bcftools mpileup`'s actual behaviour (raw AD, no GVCF rounding) or keep the current binomial draw and just route the GTs/ADs into two files.
- [ ] Rebuild `tests/test_data/` fixtures as panel+admix pairs; update all 21 single-VCF integration tests and the multi-donor / LOD / detection fixtures.
- [ ] Re-run all paper validation scripts (`paper/scripts/run_*.py`) against pipeline-style synthetic data and update facts CSVs, figures, methods text, and any results numbers that move.
- [ ] Remove the back-compat branch in `_resolve_vcf_inputs` and the `test_monitor_two_vcf_mode` parity test (becomes the only mode).
- [ ] Sweep for residual single-VCF assumptions in `doc/`, `claude/`, `README.md`.

Gated on: Step 11 re-run landing (need to see that the two-VCF results agree with or improve on the old single-VCF batch before pulling the floor out from under the existing validation).

Decision required before any of the above ships — this is the explicit "do we even want to do this" step.

---

## Step 23: Widen Force-Output Panel to Recover Marginal Markers ✅ COMPLETE

Landed 2026-05-29 (see [`2026-05-29_wider_panel_validation_notes.md`](2026-05-29_wider_panel_validation_notes.md) for the day's full write-up, before/after numbers, file inventory, and open items).

What was actually needed (the wider panel alone wasn't enough):

1. gnomAD v4.1-derived panel build: `scripts/build_force_output_panel.sh` + `pipeline/gnomad_refseq_to_hg38_chrs.tsv`. Output panel at `output/union_sid_haem_gnomad_af05.vcf.gz` (258 sites) recommended; af01 alternative kept.
2. Pipeline fix: `bcftools call -A` in `pipeline/Snakefile:401` so the admix VCF preserves the panel ALT at hom-ref sites (without this ~48 informative SNPs/patient were lost in the join).
3. Pipeline fix: `-e 'ALT="."'` in `panel_tsv` rule so force-output REF-only rows don't reach `bcftools call -A` (where they trigger malformed PL).
4. allomix fix: skip indels in `parse_vcf` (pileup can't count indel reads the way GATK's local reassembly does).
5. allomix fix: GT/AD consistency check on host/donor (drops GATK miscalls where the called het has AD VAF < 0.35 or > 0.65, etc.).

run9 verification (`output/validation_run9/batch.tsv`): n_informative up across the board (BHOA 103→144, NDAD 91→146, GBRI 85→128, PCAH 95→132, BCOL 98→150 etc.), donor% matches flow truth on every sample, host-presence detection ("any detected") preserved with stronger p-values, full test suite 320/320.

### Open items rolled out of this step

- **Post-doc follow-up (user owns):** share run9 marker-count gains and "any detection" results with the post-doc. Notes file has the table to send. LNAN host-presence p=0.16 is borderline — discuss treatment.
- **Possible Step 24 — wider-panel rho calibration.** NDAD, BHOA, PCAH show `gof_pval = 0.0000` in run9 → QC=REVIEW. The chimerism fraction is fine (matches flow); the issue is residual variance exceeding beta-binomial expectation because the wider panel has more diverse per-marker behaviour. Refit `rho` via `scripts/diagnose_sample.py` on a run9 sample and check whether a recalibrated global `rho` resolves it; if not, this folds into the per-marker-type overdispersion work already noted under Step 21.
- **Comparison plot refresh.** `output/run1_vs_run2_vs_run3_vs_run6.png` should become run1-vs-run2-vs-run3-vs-run9 once the post-doc signs off.
- **Plan tidying** (cosmetic): Step 22 section above this one still has its original "if we commit to pileup-only, the work to do:" subsection embedded under the ✅ COMPLETE header; remove next time someone is in the file.

---

## Step 24: Calibrate Wider-Panel Overdispersion (rho) for run9 REVIEW Samples 🔲 TODO

Formalises the item earmarked under Step 23. NDAD, BHOA, PCAH come up `gof_pval = 0.0000` in run9 (QC=REVIEW): the chimerism fraction matches flow, but residual per-marker variance exceeds beta-binomial expectation. Step 27 partly explains why: some "outlier" markers carry real localised recipient signal from clonal LOH, not noise, so a single global `rho` may not be the right fix. Options: refit `rho` on a run9 sample via `scripts/diagnose_sample.py`; move to per-marker-type / per-locus overdispersion (Step 21 TODO); or mask the genuine clonal-imbalance markers (Step 27) before fitting. Decide after Step 27's timepoint test clarifies which markers are real signal vs noise.

## Step 25: Host-Presence Visualisation and Per-Marker Diagnostics ✅ COMPLETE (2026-06-01)

A suite of standalone diagnostic / validation plots for the host-presence detector (Step 20), plus the per-marker tooling that drove the Step 27 discovery. All scripts are standalone (not in the package), matching `scripts/plot_chimerism_comparison.py`.

- `scripts/plot_presence_lod_curve.py` — simulated detection probability vs spiked level, x-axis in donor % (log-spaced by distance from 100%), two panels (binomial vs beta-binomial overdispersion) showing the rho penalty. Reads `output/facts/presence_lod_*.csv`.
- `scripts/plot_host_presence.py` — run9 cohort presence forest. Donor % on the Y axis (house style, log-from-100, full donor at top), samples across X, exact presence p-values listed along the bottom coloured by call. A p-value bar was rejected on purpose (its magnitude saturates; only the alpha crossing matters for the call).
- `scripts/plot_host_presence_per_marker.py` — per donor-homozygous marker, the dose-normalised implied host fraction (background-subtracted) vs marker rank, with the pooled MLE line and a count of markers whose CI excludes the pooled value ("off the line" = overdispersion / a few drivers). Recomputed MLE matches `batch.tsv` exactly.
- `scripts/host_presence_manhattan.py` — genomic-position ("Manhattan") view of per-marker implied host fraction, chromosome-banded, with a per-chromosome mean line and nearest-gene labels on Bonferroni-upregulated markers. Local-only (x-axis shows coordinates). Uses `output/refseq109_genes.bed` (protein-coding genes extracted once from `/data/cdot_data/downloads/Homo_sapiens_GRCh38_RefSeq_109.gff3.gz`).
- `scripts/host_presence_markers_vcf.py` — per-sample VCF of donor-homozygous markers with INFO annotations (HOSTY/DP/DOSE/RAWVAF/IMPLIEDF/POOLEDF/FOLD/PUP) and an UPREG flag, sorted by position for VEP / driver-panel intersection. Local-only (coordinates).
- Comparison-plot integration: `scripts/plot_chimerism_comparison.py` now rings the primary-run points green/grey by host-presence call and lists the presence p along the bottom (reads `host_present_p` from `batch.tsv`). `output/run2_run3_run9_presence.png` is the run9 example (partially addresses the Step 23 "comparison plot refresh" item, with run2/run3/run9 rather than run1).

Convention settled this session: donor % wherever a measured chimerism value is plotted (Graph 1 X-axis, Graph 2 Y-axis, the comparison-plot Y-axis); host fraction only in the per-marker diagnostics where the spread of the small signal is the point.

### Open: LoD overlay full re-run 🟡
Presence LoD was overlaid on the panel-size LoD figure (`fig5_lod_curves.png`). Code done: `paper/scripts/run_presence_lod.py` gained `--host-fractions` (so the LoD resolves across the low-depth / small-panel cells); `paper/scripts/plot_lod_curves.py` gained `--presence-summary`; `paper/scripts/run_lod_validation.py` gained `--error-rate` (was hardcoded 0.01). A quick 20-blank sweep produced `output/facts/presence_lod_bypanel_summary.csv` and `output/fig5_lod_curves_with_presence.png` as proof of concept. NOT a full re-run.
- [ ] **Decision before the full re-run:** the chimerism LoD curves are at 1% simulated error; presence is far more error-sensitive (its LoD collapses to ~14% at 1% error). For a fair overlay both sweeps must run at the same, realistic per-site error. Pick the error rate (ideally the empirical panel value from Step 14), then re-run both `run_lod_validation.py` and `run_presence_lod.py` at that rate and regenerate. Expensive job; warn before triggering.

## Step 26: Sex-Chromosome Handling ✅ COMPLETE (2026-06-01)

Sex and mitochondrial contigs (X/Y/M) are unreliable for chimerism in sex-mismatched transplants (recipient/donor allele dosage on chrX/chrY is wrong).
- `genotype.classify_markers` gained `use_sex_chroms: bool = False` (default excludes X/Y/M) and reports `n_sex_chrom_excluded`. Helper `genotype.is_sex_chrom`.
- CLI: `--use-sex-chroms` on `monitor` / `timeline` (default off); prints the dropped count per sample to stderr.
- Reportable results default to excluding sex chroms; the diagnostic scripts (Step 25) force `use_sex_chroms=True` so chrX stays visible for investigation, just not in reportable output.
- Cost in run9: 5 / 6 / 7 chrX informative markers (no chrY) for NDAD / BHOA / PCAH.
- Re-enable per run once recipient and donor sex are confirmed to match. Sex is being added to the project xls files (user owns).

## Step 27: Per-Marker Clonal Allelic Imbalance (TP53/17p LOH) as an Incidental Relapse Signal 🟡 HYPOTHESIS (2026-06-01)

The most interesting find of the session. The host-presence detector flagged markers (via the Step 25 diagnostics) where the recipient-distinguishing allele reads several-fold higher than the genome-wide recipient fraction. The Bonferroni-upregulated markers were TP53 (chr17:7676483) in all three REVIEW samples (NDAD, BHOA, PCAH) and BCOR (chrX) in PCAH.

**It is real recipient signal, not a mapping artifact.** Proven by a forced `samtools mpileup` (mapq>=20, baseq>=20, `--excl-flags 1796`, matching the pipeline's admix pileup) on the two TRUE donor BAMs (now at `output/bams/`) at chr17:7676483: the host-free donors show ~0% of the recipient allele (QUDO donor 0.00%, PCAH donor 0.03% at >3000x), while patient admix samples show 2-3%. GATK AD on the donors also reads ~0% but cannot see a pileup artifact (local reassembly discards the mismapped reads), so the pileup is the definitive control. Note: the earlier "pure donor" samples (GBRI/PNOL) were TP1-proxy donors, NOT verified host-free, and must not be used as controls.

**Mechanism hypothesis:** the residual recipient cells are the relapsing TP53-aberrant clone carrying 17p CN-LOH (acquired UPD) or amplification, so the retained recipient allele is locally over-represented relative to the genome-wide recipient fraction. Loss of the wild-type allele (function down) and increased reads of the retained allele (signal up) are the same event. CN-LOH alone gives ~2x; the observed ~2-5x (after dose correction) implies copy gain on top, consistent with complex-karyotype 17p amplification. Recurrence across three unrelated patients fits a relapse-surveillance cohort enriched for TP53/17p-aberrant disease, not an artifact.

**Literature search (deep-research, 2026-06-01): the specific idea appears novel / underexplored.**
- Confounder is old and robust: LOH/CNV at marker loci distorts chimerism; STR loci with gain/loss are routinely excluded (Vietz/Lin PMC3219907; ASH 2018 Blood 132:5135; UK NEQAS Br J Haematol 2014, 10.1111/bjh.13073).
- Feature version exists but narrowly: STR gain/loss loci as a relapse marker (Vietz/Lin, n=4 case series); and NGS HLA-loss assays (Devyser HLA Loss, PMC12259423; HLA-CLN, 10.1177/09636897221102902; concept from Vago Blood 2009;115:3158) that detect CN-LOH/UPD of 6p by per-marker imbalance over the HLA region BY DESIGN.
- TP53/17p biology firmly established but in diagnostic/MDS papers, not chimerism (Jasek Leukemia 2009, PMID 19759556; Sugimoto Leuk Res 2016, PMID 26851439).
- Mixed chimerism + MRD integration uses genome-wide-average chimerism plus a separate MRD readout, never per-marker LOH (Kanaan Blood Adv 2023; FIGARO/Pearce Blood Adv 2023; Lee/Jo J Clin Med 2019, 121 SNPs averaged).
- Mainstream NGS-SNP chimerism treats chromosomal deletion as a false-negative to dilute out with more markers (AlloSeq HCT review PMC10380370; Kim/Stahl J Mol Diagn 2024, PMID 38307253, no mention of CNV/LOH/TP53). The opposite of exploiting it.
- The gap: nobody applies incidental per-marker clonal CN-LOH detection at non-HLA driver loci (TP53/17p) within a general SNP/NGS chimerism panel as a relapse feature. The HLA-loss assays are the precedent/roadmap (they prove the concept works in NGS), not competition. allomix already exposes per-marker chimerism genome-wide (Step 25), so generalising it is the novel step.
- Direction nuance: the literature focuses on the deletion / reduced-signal (false-negative) direction; the over-representation (CN-LOH up-signal) direction we observe is even less characterised.

**Open tests / next steps (in order of decisiveness):**
- [ ] **Timepoint test (the discriminator).** Pull the TP53 recipient-allele VAF across serial timepoints for relapsers; the MRD signature is a rise over time, ideally leading the genome-wide chimerism. Data is local. This makes or breaks the hypothesis.
- [ ] **Segmental check.** Real 17p LOH lifts a run of adjacent 17p markers, not a lone SNP. Read the per-chromosome mean line on the Manhattan (Step 25). A lone spike forces a focal / structural explanation instead.
- [ ] **Ground truth.** Cross-check the patients' molecular / cytogenetic reports for TP53 mutation, 17p loss, or complex karyotype.
- [ ] **BCOR/chrX.** Run the same true-donor pileup on the BCOR chrX sites (chrX has the sex-mismatch dosage caveat, Step 26).
- [ ] **Per-marker null.** Any rigorous version must separate true clonal imbalance from per-marker technical bias / dropout (STR uses per-locus baselines; the HLA assays compare region-vs-rest). This is the same calibrated per-marker null as Step 14 (per-site error table) and Step 24 (overdispersion), so the machinery that cleans up quantification is what would make the relapse feature rigorous.
- [ ] Proper PubMed / Scholar plus patent sweep before claiming novelty (the search above is web-based; "appears novel" is not "proven novel"). Two sub-claims were refuted in verification and must not be cited: a "12% vs 2%" HLA-loss discordance figure, and a "HLA-CLN is a general-LOH precedent" framing.

Implication for the rest of the tool: these markers are not noise to blacklist. They may detect a real residual clone (potentially useful), but they remain a quantification confound that feeds the Step 24 overdispersion / REVIEW flags. Handling is a policy choice once the timepoint test resolves whether the signal is what we think it is.
