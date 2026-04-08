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

## Step 13: Publication 🟡 IN PROGRESS

- Method paper describing the approach — framework set up with vibepaper
- Target journal: Journal of Molecular Diagnostics (Technical Advance)
- Paper sections in `paper/`, analysis scripts in `paper/scripts/`
- Cite: Crysup & Woerner 2022 (Demixtify MLE framework), Vynck et al. (bias correction)
- In silico validation complete (depth series, relatedness, bias correction, multi-donor)
- Multi-donor paper updates complete (methods, results, discussion, abstract, figures)
- Simulation calibrated from empirical panel characterisation (210 VCFs, 18,047 samples)
- Validation with real samples (Step 11) still needed
- Open-source tool release (MIT license, PyPI)
