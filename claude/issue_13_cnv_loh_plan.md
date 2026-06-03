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

`paper/scripts/run_cnv_loh_validation.py`. Genotypes and per-marker capture
biases are fixed per (relatedness, replicate) and reused across every cell, so
within a replicate only the aberration and the sequencing draw change. The
no-aberration baseline is run once per replicate and shared across kinds.

| Axis | Values |
|------|--------|
| Relatedness | unrelated, sibling |
| Aberration kind | cnloh, deletion, gain |
| Burden (fraction of eligible markers) | 0.0 (baseline), 0.1, 0.25, 0.5 |
| Clonal fraction | 0.5, 1.0 |
| True donor fraction | 0.2, 0.5, 0.8, 0.9, 0.95, 0.99 |
| Replicates | 20 |
| Markers / depth | 100 / 1000x |

(For `cnloh` only het markers are eligible; for `deletion`/`gain` every marker
is, which is why they affect roughly twice as many markers at matched burden.)

Fixed noise model matches the other validation scripts (error 0.01, locus
dropout 0.016, depth CV 0.43, realistic heavy-tailed biases).

Metrics per cell: MAE, signed bias, RMSE, 95% CI coverage, mean markers
affected, and **mean markers flagged by the estimator's 3-SD outlier rule**.

Outputs: `output/facts/cnv_loh_{raw,summary,headline}.csv` and
`output/facts/fig_cnv_loh.png`. Wired into `paper/Snakefile`
(`cnv_loh_validation`, `cnv_loh_plot`).

## Results

20 replicates, 100 markers, 1000x, averaged over the six true donor fractions.
Pure clone (clonal_fraction = 1.0). Baseline (no aberration): MAE 0.0016
(unrelated) / 0.0023 (sibling), coverage 0.97 / 0.96.

| Relatedness | Kind | Burden | MAE | Signed bias | 95% CI cov | Affected | Flagged |
|-------------|------|-------:|----:|------------:|-----------:|---------:|--------:|
| unrelated | cnloh    | 0.10 | 0.0114 | +0.010 | 0.89 | 3.7 | 1.4 |
| unrelated | cnloh    | 0.25 | 0.0225 | +0.020 | 0.79 | 11.4 | 1.3 |
| unrelated | cnloh    | 0.50 | 0.0409 | +0.039 | 0.62 | 19.9 | 0.4 |
| unrelated | deletion | 0.10 | 0.0113 | +0.011 | 0.77 | 8.7 | 1.2 |
| unrelated | deletion | 0.25 | 0.0391 | +0.039 | 0.22 | 25.4 | 0.7 |
| unrelated | deletion | 0.50 | 0.0640 | +0.064 | 0.11 | 50.8 | 0.3 |
| unrelated | gain     | 0.10 | 0.0047 | -0.004 | 0.82 | 9.7 | 1.3 |
| unrelated | gain     | 0.25 | 0.0109 | -0.011 | 0.47 | 23.8 | 1.0 |
| unrelated | gain     | 0.50 | 0.0210 | -0.021 | 0.20 | 46.7 | 0.4 |
| sibling | cnloh    | 0.10 | 0.0143 | +0.007 | 0.87 | 4.0 | 0.8 |
| sibling | cnloh    | 0.25 | 0.0276 | +0.023 | 0.80 | 11.8 | 0.5 |
| sibling | cnloh    | 0.50 | 0.0401 | +0.036 | 0.78 | 20.8 | 0.2 |
| sibling | deletion | 0.10 | 0.0127 | +0.012 | 0.88 | 9.2 | 0.7 |
| sibling | deletion | 0.25 | 0.0410 | +0.041 | 0.50 | 24.7 | 0.4 |
| sibling | deletion | 0.50 | 0.0697 | +0.070 | 0.22 | 48.5 | 0.1 |
| sibling | gain     | 0.10 | 0.0048 | -0.004 | 0.93 | 9.1 | 0.6 |
| sibling | gain     | 0.25 | 0.0087 | -0.008 | 0.78 | 23.6 | 0.6 |
| sibling | gain     | 0.50 | 0.0239 | -0.024 | 0.37 | 51.0 | 0.2 |

MAE inflation at the highest burden vs baseline: deletion **30-40x**, CN-LoH
**18-26x**, gain **10-13x**. Figure: `output/facts/fig_cnv_loh.png`.

Five things to note:

1. **Even a modest burden matters.** 10% burden (a handful of markers on a
   100-marker panel) already raises MAE 3-7x and drops coverage. Real arm-level
   events on a sparse genome-wide panel land in exactly this few-markers regime.
2. **Deletion is the worst kind.** It both skews allele balance and halves the
   locus DNA mass, and it hits every genotype (not just hets), so it affects
   ~2x as many markers and collapses coverage hardest (0.11-0.22 at burden 0.5).
3. **The bias direction depends on the kind.** Deletion and CN-LoH bias *upward*
   (donor overestimated -> residual host **underestimated**, the dangerous
   under-call direction for relapse). Gain biases *downward* (donor
   underestimated -> host **overestimated**, a false-alarm direction). A deletion
   removes host DNA so the locus looks more donor; a gain adds host DNA so it
   looks less donor.
4. **Gain is the mildest** because tripling host copies is a smaller relative
   DNA-mass change (x1.5) than halving it (x0.5), and the allele skew is gentler.
5. **Siblings are generally hit at least as hard** as unrelated donors, because
   they have fewer informative markers so each rogue marker carries more
   leverage.

Partial clone (clonal_fraction = 0.5) shows the same pattern at roughly half the
magnitude, as expected from the mixture model.

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

The validation sweep records standard *and* robust per cell in one pass (no
separate runs), so `cnv_loh_summary.csv` carries `mae`/`mae_robust`/
`ci_coverage`/`ci_coverage_robust`, and `plot_cnv_loh.py` overlays them (solid =
standard, dashed = robust). Pure clone, avg over true fractions, 20 reps:

| Rel | Kind | Burden | MAE std → robust | CI cov std → robust | dropped |
|-----|------|-------:|------------------|---------------------|--------:|
| unrelated | cnloh    | 0.10 | 0.0114 → 0.0048 | 0.89 → 0.95 | 2.4 |
| unrelated | cnloh    | 0.25 | 0.0225 → 0.0122 | 0.79 → 0.84 | 4.1 |
| unrelated | deletion | 0.10 | 0.0113 → 0.0064 | 0.77 → 0.85 | 2.2 |
| unrelated | deletion | 0.25 | 0.0391 → 0.0297 | 0.22 → 0.31 | 2.4 |
| unrelated | gain     | 0.25 | 0.0109 → 0.0085 | 0.47 → 0.50 | 4.0 |
| sibling | cnloh    | 0.25 | 0.0276 → 0.0180 | 0.80 → 0.81 | 2.2 |
| sibling | deletion | 0.25 | 0.0410 → 0.0344 | 0.50 → 0.52 | 0.7 |

At low-to-moderate burden the refit cuts MAE roughly 1.5-2.4x and recovers
coverage, most strongly for CN-LoH and deletion in unrelated donors. It helps
less for siblings (fewer informative markers, so the marker floor binds and the
engage-gate trips less often).

**Neutral on clean data (why it is safe as the default).** At burden 0 the gate
never engages, so the estimate is byte-identical: MAE 0.0016 → 0.0016 (unrelated),
0.0023 → 0.0023 (sibling). Turning it on costs nothing when there is nothing to
fix, which is what keeps the already-validated cohort stable.

**Limit at extreme burden.** At 50% burden it barely helps (deletion 0.064 →
0.063) and gain can even tick worse, because the aberrant markers are no longer a
minority so median/MAD has no clean reference set. This is exactly why the high
drop fraction is flagged REVIEW rather than trusted: the regime needs orthogonal
CNV information (section 3), not more trimming.

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
