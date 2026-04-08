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
- `data/joint_called_example.vcf` — 114 samples, 9 markers, confirms joint-calling provides ref+alt AD at all sites

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

Full analysis in `claude/step4_reference_tool_analysis.md`. 9 repos examined.

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

### Test coverage: 192 tests passing

- 29 genotype tests (parsing real VCFs, classification, filtering)
- 55 chimerism tests (MLE math, estimation accuracy at multiple fractions, CI coverage, edge cases)
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
```

---

## Step 7: Implement Multi-Donor Support 🔲 TODO

**Goal:** Extend chimerism estimation to host + 2 donors.

- 2D grid search over (f_d1, f_d2) with constraint f_d1 + f_d2 ≤ 1
- Nelder-Mead refinement (scipy)
- Profile likelihood CIs per donor (chi-square 2df)
- A marker is informative for donor_i if host and donor_i genotypes differ
- `InformativeMarker.donor_gts` already stores multiple donor genotypes
- Update CLI: `--donor d1.vcf --donor d2.vcf` already accepted
- Update report: output per-donor fractions

---

## Step 8: Implement Bias Correction 🔲 TODO

**Goal:** Per-marker amplification bias correction (Vynck et al.).

- `src/allomix/bias.py` — estimate bias from training VCFs, apply correction
- Bias = median(VAF_het - 0.5) per marker across a training set
- Analytic correction formulas for each Vynck marker type
- CLI: `allomix estimate-bias --vcfs *.vcf.gz -o bias_table.tsv`
- CLI: `allomix monitor --bias-table bias.tsv ...`
- Optional — tool works without bias correction, just less accurate

---

## Step 9: In-Silico Validation 🟡 PARTIALLY COMPLETE

### Done

Test data infrastructure and initial validation are complete:

- `scripts/make_synthetic_genotypes.py` — generates 100-SNP synthetic host + donor VCFs (80 informative markers)
- `scripts/generate_test_data.py` — blends host + donor at specified fractions, outputs `host_X_donor_Y.vcf` naming. Supports `--bias-sd` for simulating per-marker capture/amplification bias.
- `scripts/generate_timeline_data.py` — generates 7-timepoint engraftment + relapse scenario (`day030_donor_95.vcf` through `day300_donor_40.vcf`)
- `scripts/run_validation.py` — reads truth table, runs allomix on each sample, computes accuracy metrics (bias, MAE, RMSE, CI coverage), produces validation_results.tsv, validation_summary.tsv, and 3 plots (scatter, residuals, CI coverage)
- `tests/test_data/` — 11 chimeric VCFs at 0–100% in 10% steps + timeline data (7 timepoints)
- Capture bias simulation added to `simulate.py` via `marker_bias_sd` parameter (0.0 = ideal, 0.02 = realistic)

### Validation Results (100 markers, 2000x depth, bias_sd=0.02)

| Metric | Value |
|--------|-------|
| Mean signed error | +0.08% |
| Mean absolute error | 0.30% |
| RMSE | 0.37% |
| Max absolute error | 0.59% |
| CI coverage rate | **55%** (target: 95%) |
| Mean CI width | 0.51% |

**Accuracy is excellent. CI coverage is the known gap.** The 95% profile likelihood CIs only cover the truth ~55% of the time because the model assumes pure binomial sampling variance but does not account for:
1. Per-marker capture/amplification bias (systematic VAF shifts, ~2% SD per marker)
2. Overdispersion beyond binomial (real sequencing data is slightly overdispersed)

**Expected fix:** Step 8 (Vynck bias correction) should correct the systematic per-marker shifts. Additionally, a variance inflation factor or empirical CI calibration will be needed to achieve proper coverage.

### Remaining

- Run at finer fractions near the sensitivity limit (0.1%, 0.5%, 1%, 2%, 5%) to characterise low-fraction performance
- Generate multiple donor-host pairs (different genotype distributions, related pairs with fewer informative markers)
- Vary depth (500x, 1000x, 2000x, 5000x) to measure depth effect on CI width and accuracy
- Vary marker count (10, 20, 50, 100) to measure panel size effect
- Implement CI calibration/inflation after Step 8 bias correction
- Validate on real data (Step 11)

---

## Step 10: VariantGrid Integration 🔲 TODO

- JSON output schema agreed with VG team
- VG stores donor/host genotypes, exports as VCF
- VG ingests allomix JSON results per patient
- VG renders timeline chart from chimerism results across timepoints
- Exact API integration TBD

---

## Step 11: Real Sample Validation 🔲 TODO

- Run allomix on real post-HSCT samples from /tau
- Compare with STR-based chimerism results (current method)
- Assess concordance, sensitivity at low fractions
- Requires running the /tau audit script first (in step2 document)

---

## Step 12: Per-Base Quality-Aware Likelihood 🔲 TODO

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

---

## Step 13: Publication 🔲 IN PROGRESS

- Method paper describing the approach — framework set up with vibepaper
- Target journal: Journal of Molecular Diagnostics (Technical Advance)
- Paper sections in `paper/`, analysis scripts in `paper/scripts/`
- Cite: Crysup & Woerner 2022 (Demixtify MLE framework), Vynck et al. (bias correction)
- In silico validation complete (depth series, relatedness, bias correction)
- Validation with real samples (Steps 9 + 11) still needed
- Open-source tool release (MIT license, PyPI)
