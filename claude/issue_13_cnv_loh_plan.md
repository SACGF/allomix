# Issue #13 — Simulation: CNV / LoH in host genome

Tracking: https://github.com/SACGF/allomix/issues/13

## Motivation

The HSCT recipient ("host") is almost always a haematological malignancy
patient (AML, MDS, ALL, ...). The residual or relapsing host clone routinely
carries somatic copy-number aberrations that the chimerism model does not
account for:

- **Copy-neutral LoH (CN-LoH / acquired uniparental disomy).** The clone
  retains two copies of one germline homolog. A host heterozygous marker
  becomes an effective homozygote with no change in copy number. Common in AML
  (e.g. UPD13q with FLT3-ITD, UPD11p).
- **Deletions / monosomy** (del5q, del7q, −7) and **gains / trisomy** (+8).
  These also change the total DNA mass the locus contributes.

Every VAF calculation in `simulate.py` and the likelihood in `chimerism.py`
assumes all contributing genomes are diploid:

```
expected_vaf = ((1 - f) * host_dose + f * donor_dose) / 2
```

CN-LoH breaks the allele-balance assumption; deletions/gains additionally break
the DNA-mass assumption, so the local mixing fraction at an aberrant locus is
no longer the genome-wide `f`. With a sparse genome-wide SNP panel only a few
markers fall inside any given event, so the failure mode is a handful of
high-leverage rogue markers per sample rather than a global shift.

This issue adds CN-LoH (and the general copy-number machinery) to the simulator
and measures the impact on chimerism estimates.

## Scope of this pass

- CN-LoH, deletion (CN1), and gain (CN3). CN-LoH is copy-neutral (allele-balance
  effect, het markers only); deletion and gain also change locus DNA mass, so
  they bias even at homozygous markers and are applied to any genotype.
- The host genotype VCF stays clean germline, so marker classification is
  unaffected. The aberration is applied only to the admixture sample, modelling
  the realistic two-phase pipeline where host GT comes from a clean reference
  and the residual host clone in the admix sample is the disease. (Corrupting
  the host reference GT itself, the "GT called from diseased sample" case, is a
  separate follow-up.)

## Model

`simulate.py` additions:

- `HostAberration(cn, alt_copies, clonal_fraction)` — the clone's state at one
  marker, as a mixture of normal diploid host (fraction `1 - clonal_fraction`)
  and aberrant clone (`clonal_fraction`).
- `cn_weighted_vaf(host_gt, donor_gts, donor_fractions, host_aberration)` —
  copy-number-weighted mixture VAF:

  ```
  VAF = sum_i frac_i * cn_i * alt_frac_i  /  sum_i frac_i * cn_i
  ```

  Reduces exactly to `expected_vaf_multi` when there is no aberration.
- `assign_cnv_aberrations(markers, fraction_affected, clonal_fraction, rng, kind)`
  — builds the clone by mutating one random germline homolog: `cnloh` duplicates
  the retained homolog (cn=2, het markers only), `deletion` drops one homolog
  (cn=1, any genotype), `gain` duplicates one homolog (cn=3, any genotype).
  `assign_cnloh_aberrations` is a thin `kind="cnloh"` wrapper.
- `blend_vcfs(..., host_aberrations=...)` — per-shared-marker aberration list
  (aligned like `fixed_biases`); affected markers use `cn_weighted_vaf`.

Worked example (host het, hom-ref donor, donor fraction f = 0.2, pure clone
retaining ALT): the diploid model expects VAF = 0.5·(1−f) = 0.4, but the CN-LoH
clone is hom-ALT so the locus reads VAF = (1−f) = 0.8. That +0.4 shift at the
affected markers is the bias the estimator has to absorb.

## Sweep design

`paper/scripts/run_cnv_loh_validation.py`. The metric is the **limit of
detection (LoD)** of the minor component, the quantity the rest of the paper
reports (CLSI EP17-A2, as in `run_lod_validation.py`), plotted on a log donor-%
axis (0.3, 0.5, 1, 2, 5, 10, 20 %) like `plot_lod_curves.py`. The recipient
clone always carries the aberration. Both low-fraction detection directions are
swept (one sweep, `--modes`):

- **`host` (relapse early-warning, the primary use):** detect the recipient
  relapse clone (the minor component, carrying the aberration) against a clean
  donor background. The blank (true host = 0) is pure donor.
- **`donor` (mixed chimerism / substantial recipient):** detect the donor (the
  minor component) against a recipient CN-LoH background. The blank (true donor
  = 0) is a pure host carrying the aberration.

For each direction: LoB = mean + 1.645·SD of the estimated minor fraction over
blanks; LoD = lowest true minor fraction at which ≥95% of replicates exceed LoB
(logistic fit of P(detected) vs log10(f)); a LoD past the 20% probed ceiling is
"above range" (undetectable). Genotypes and capture biases are fixed per
(relatedness, replicate); standard and robust (`--robust auto`) are computed per
cell in one pass.

| Axis | Values |
|------|--------|
| Detection direction | host (relapse), donor (mixed chimerism) |
| Relatedness | unrelated, sibling |
| Aberration kind | cnloh, deletion, gain |
| Burden (fraction of eligible markers) | 0.0 (baseline), 0.1, 0.25, 0.5 |
| Clonal fraction | 1.0 (pure clone) |
| True minor fraction | 0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2 |
| Replicates | 30 |
| Markers / depth | 100 / 1000x |

(For `cnloh` only het markers are eligible; for `deletion`/`gain` every marker
is, so they affect ~2x as many markers at matched burden.) Fixed noise model
matches the other scripts (error 0.01, locus dropout 0.016, depth CV 0.43,
realistic heavy-tailed biases).

Outputs: `output/facts/cnv_loh_{raw,summary,headline}.csv` (with a `mode` column)
and `output/facts/fig_cnv_loh.png` (2 rows = directions, 3 cols = kinds; log
donor-% LoD axis; std solid / robust dashed). Wired into `paper/Snakefile`
(`cnv_loh_validation`, `cnv_loh_plot`).

## Results

30 replicates, 100 markers, 1000x, pure clone. The two detection directions
behave **very differently**, which is why both are reported. `>20%` = LoD above
the probed ceiling (component undetectable). std → robust refit.

### Relapse detection (host minor; the early-warning use)

Baseline relapse LoD **0.49%** (unrelated) / **0.78%** (sibling). Across every
kind and burden it stays put (LoD ~0.2-1.6%); robust refit makes no material
difference. Recipient CN-LoH/CNV does **not** degrade relapse detection.

| Rel | Kind | Burden | relapse LoD std → robust |
|-----|------|-------:|--------------------------|
| unrelated | deletion | 0.50 | 0.68% → 0.80% |
| sibling | cnloh | 0.50 | 1.25% → 1.59% |
| sibling | deletion | 0.50 | 1.45% → 1.51% |

Why: the aberration rides the *minor* component, the blank (true host = 0) is
pure donor so the LoB is clean, and the strongest detector markers are
host-homozygous (host 1/1 vs donor 0/0), which CN-LoH (a het-only effect) does
not touch. A germline-referenced relapse stays detectable despite the clone's
aberrations. (Deletion in siblings is the only hint of degradation, ~2x.)

### Donor detection (host major + CN-LoH; mixed chimerism)

Baseline donor LoD **0.52%** (unrelated) / **1.14%** (sibling). The recipient's
CN-LoH/CNV background badly inflates it:

| Rel | Kind | Burden | donor LoD std → robust |
|-----|------|-------:|------------------------|
| unrelated | cnloh    | 0.10 | 18.4% → 9.5% |
| unrelated | cnloh    | 0.25 | >20% → >20% |
| unrelated | deletion | 0.10 | >20% → 6.6% |
| unrelated | deletion | 0.25 | >20% → >20% |
| unrelated | gain     | 0.25 | 1.8% → 0.54% |
| unrelated | gain     | 0.50 | 4.4% → 0.69% |
| sibling | cnloh    | 0.10 | 17.3% → >20% |
| sibling | deletion | 0.10 | 17.3% → 16.3% |
| sibling | gain     | 0.50 | 4.2% → 1.4% |

1. **CN-LoH/deletion wreck the donor LoD**: 0.5% baseline → ~17-18% at 10%
   burden, undetectable (>20%) by 25%. Deletion is worst (allele skew + DNA mass
   + hits every genotype); gain is mildest (LoD only to a few %).
2. **Robust refit recovers the gain case strongly** (4.4% → 0.7%) and rescues
   low-burden deletion in unrelated donors (>20% → 6.6%); it cannot rescue
   CN-LoH/deletion once burden ≥25% (no clean marker majority) — those stay
   undetectable and are flagged REVIEW, not trusted.
3. **Neutral on clean data**: baseline LoD essentially unchanged by robust.

Figure: `output/facts/fig_cnv_loh.png` (top row relapse LoD, bottom row donor
LoD, 3 kinds, log donor-% axis, std solid / robust dashed).

**Clinical read:** routine relapse early-warning is robust to recipient
CN-LoH/CNV, but measuring the *donor* fraction in a sample with a substantial
CN-LoH-bearing recipient (mixed chimerism) can be badly biased — there the
robust refit and the REVIEW flag matter. A complementary point-estimate (MAE)
analysis in the donor-dominant regime is consistent: deletion/CN-LoH push the
donor fraction up (residual host under-called), gain pushes it down.

## Key finding: rogue markers are not rejected

The estimator computes per-marker residuals at the fitted MLE and flags those
beyond 3 SD (`chimerism.py:_compute_per_marker_results`), but it **never refits
without them** — `n_markers_used` only counts survivors. When a meaningful
fraction of markers is biased the same way, those markers inflate the residual
SD and shift the mean, so almost none exceed 3 SD. The sweep confirms this
across all three kinds: at the highest burden ~20-50 markers are biased but the
flag catches ~0.1-0.4 of them on average, and the point estimate is pulled
accordingly. The outlier path is cosmetic here, not protective.

## How to fix / detect it

Three families of approach, in rough order of effort. The simulator now provides
the ground-truth generator to validate any of them.

### 1. Robust refit (IMPLEMENTED, the default mitigation)

Iterate: fit, compute per-marker residuals, drop those beyond `k` robust SDs
using **median/MAD** (not mean/SD, so a cluster of contaminating markers cannot
set the scale the way it defeats the current 3-SD flag), refit on the survivors
(f, rho, CI, detection limits, GoF), repeat until the set stabilises. Shipped as
`estimate_single_donor_bb` / `estimate_multi_donor` `robust=` and CLI
`--robust {off,auto,force}` / `--robust-k` (`_robust_refit`; k=3.5, ≤5 iters).

This matters only in the **donor-detection (mixed chimerism)** regime; relapse
detection needs no help. The sweep records standard *and* robust per cell in one
pass; `cnv_loh_summary.csv` carries `lod_std`/`lod_robust` and `plot_cnv_loh.py`
overlays them. The donor-LoD recovery (Results table above), short version:

- **Gain:** strong recovery, LoD back near baseline (unrelated 25% burden
  1.8% → 0.54%, 50% burden 4.4% → 0.69%).
- **Deletion / CN-LoH at low burden, unrelated:** partial recovery (deletion 10%
  >20% → 6.6%; CN-LoH 10% 18.4% → 9.5%).
- **CN-LoH / deletion at ≥25% burden, and siblings:** little or no recovery
  (no clean marker majority; the marker floor binds for siblings).

**Neutral on clean data (why it is safe as the default).** At burden 0 the gate
never engages, so the LoD is unchanged within replicate noise (baseline donor
LoD ~0.5% unrelated / ~1.1% sibling, std ≈ robust). Turning it on costs nothing
when there is nothing to fix, which keeps the already-validated cohort stable.

**Limit at high burden.** Once burden ≥25% for CN-LoH/deletion the donor is
undetectable (LoD >20%) with or without the refit, because the aberrant markers
are no longer a minority and median/MAD has no clean reference set. This is
exactly why a high drop fraction is flagged REVIEW rather than trusted: the
regime needs orthogonal CNV information (section 3), not more trimming.

### Packaging: default-on, gated, and flagged (IMPLEMENTED)

- **Default on (`--robust auto`), not opt-in.** The users who most need it
  (recipient with CNV/LoH) would not know to set a flag. The 3-SD code already
  computed `included`/`n_markers_used` but never refit; this completes it for the
  single- and multi-donor paths (threaded via `analysis.analyse_sample`).
- **Gated so clean runs are unchanged.** `auto` engages only when the first
  trimming pass finds more outliers than `_robust_trigger(n) = max(3,
  ceil(0.03 n))`; below that the result is returned untouched (byte-identical),
  keeping clean and already-validated samples stable. `force` drops the gate
  (trigger 1) and the marker floor (down to `ROBUST_HARD_MIN`) for experiments.
- **Floor so sparse panels are protected.** `auto` never trims below
  `ROBUST_MIN_MARKERS` (15) survivors.
- **Always flagged.** Dropped markers are recorded (`included=False`,
  `n_robust_excluded`, `robust_drop_fraction`) with a QC warning, and a drop
  fraction above `ROBUST_REVIEW_FRACTION` (0.15) promotes QC status to REVIEW
  (alongside the existing GoF trigger in `qc.assess_quality`).

Net: behaves like "always run", but the gate makes it "detect and kick in", and
it always flags. **Next:** re-run the real-data validation batch with
`--robust auto` to confirm clean-data neutrality on actual samples (expected:
near-identical numbers, a few new REVIEW flags on genuinely aberrant samples).

### 2. Robust loss in the MLE

Fold the down-weighting into the likelihood itself (Huber / trimmed / t-distributed
per-marker loss) so a few wild markers cannot dominate, instead of the hard
drop-and-refit. More invasive (touches `chimerism.py` closures, ties into the
Step 12 per-marker-context refactor) but smoother and avoids the hard cutoff.

### 3. Detect the aberration directly (best, complements 1-2)

The aberrant markers carry their own signature, so flag them rather than only
absorbing them:

- **Allele-balance signature.** A host-het marker whose admix VAF sits near 0 or
  1, where the donor genotype cannot explain it, is a CN-LoH/deletion candidate.
  Surface as a per-marker QC flag (overlaps issue #14, marker similarity QC).
- **Use the host genotyping sample's own depth.** If the host reference VCF (or a
  paired tumour) is sequenced, a per-marker copy-number/B-allele-frequency
  profile localises deletions/gains directly; those markers can be excluded or
  CN-corrected via `cn_weighted_vaf` (the math already supports CN ≠ 2).
- **Spatial clustering.** Real CNV/LoH events are contiguous genomic blocks. On a
  dense panel, several adjacent informative markers sharing a consistent VAF
  shift is a much stronger signal than any single marker; a positional scan would
  catch arm-level events the per-marker flag misses.
- **Goodness-of-fit gate.** Heavy aberration burden inflates the beta-binomial
  residual / lowers `gof_pval`; route those samples to manual REVIEW (consistent
  with how Steps 24/27/28 treat artifact-driven misfit).

### Recommended sequence

Ship robust refit (1) as the cheap default with a high-drop-fraction warning,
add the per-marker allele-balance QC flag (3) so users see *which* markers are
suspect, and reserve the host-depth CN profile and robust-loss work for when
paired tumour/host data or a denser panel makes them worthwhile.
