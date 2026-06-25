# Step 30 design note: per-marker contamination correction, gated by a per-flowcell table

Status: design, backed by the prototype in `output/figure_review/` (see
`srp434573_figure_review_findings.md`). Not yet implemented in the package.

## Problem

On co-pooled panels (index hopping on a patterned flowcell), donor-homozygous markers
carry extra host-allele reads from co-pooled genomes that happen to carry the host
(donor-absent) allele. The current contamination handling is a per-sample scalar, which
does not localize this: it taxes every donor-absent marker by the average. Consequences
measured on SRP434573 (real reads):
- At true-0%-host samples the MLE shows a positive floor (up to 0.165%) on the
  contaminated mixtures, statistically significant (0 outside the 95% CI), while the
  presence test correctly calls not-detected.
- At low dilutions the same contamination inflates the estimate upward at every level
  (e.g. M3->F2: 1% read as 1.18%, 2.5% as 3.05%, 10% as 11.1%).
- The contamination is donor-hom-marker-specific and scales with co-pooled carrier dose
  (M3->F2 host-allele background: 0.29% at 1 carrier rising to 0.85% at 5), and is missed
  by the consensus-hom contamination floor (which reads ~0 there).

## Mechanism (what Step 30 does)

Predict the per-marker contamination on donor-hom markers from the co-pooled carrier dose
and subtract it from the host-allele count before the MLE. The host signal is identical at
every donor-hom marker (independent of carrier count); contamination scales with carrier
count. So per marker:

    host_allele_reads_corrected = host_allele_reads - slope * n_carriers * depth

where `slope` is the per-carrier contamination rate and `n_carriers` is the number of
co-pooled individuals carrying the host allele at that site. Only the dose-dependent part
is subtracted; the flat error floor stays handled by the existing per-site error model
(error table #23), so it is not double-counted.

## Prototype evidence (output/figure_review/step30_prototype.py, step30_sweep.py)

Pool-level slope (one slope per pool), applied to every timepoint, across all 10 mixtures:
- Zero-host floor removed in 10/10 mixtures.
- Dilution-series accuracy: 5 HELP, 2 neutral, 3 HARM. The 3 HARM cases (F2->M1, F2->M2,
  M1->M2) are clean/low-contamination mixtures with noise-level slopes (0.013-0.018%/carrier)
  where subtracting a spurious slope removes a little real signal.

So Step 30 helps contaminated runs and always kills the false floor, but UNCONDITIONAL
application harms clean runs. It needs a gate.

## The gate problem, and why a per-flowcell table is the answer

A hardcoded slope threshold (the SRP434573 data suggests ~0.03-0.04%/carrier separates
HELP from HARM) is tuned to these samples and would not transfer to another panel, depth,
or flowcell. That is the wrong kind of knob.

Contamination/index hopping is a PER-FLOWCELL property: it depends on which samples were
co-pooled in that run and how much cross-talk that flowcell produced. So the gate and the
correction magnitude should both be measured per run, from the run's own data, not
hardcoded. Concretely, build a per-flowcell contamination table:

1. Inputs: the flowcell's joint-called germline genotypes for all co-pooled samples (the
   same joint-call output that already feeds the #23 error table) plus the admix AD.
2. At CONSENSUS-HOMOZYGOUS sites (no contributor carries the minor allele, so the minor
   allele is pure background), pool across the flowcell and fit the minor-allele fraction
   as a function of co-pooled carrier count. This gives, per run:
   - c(0): the sequencing-error floor (no-carrier sites).
   - slope: the per-carrier contamination rate for this flowcell.
   - a SIGNIFICANCE for the slope (is the dose-response real given this run's depth/noise?).
3. Gate = the consensus-hom slope's significance on THIS flowcell. Clean run -> slope not
   significant -> no correction (no harm). Contaminated run -> significant -> correct.
   (Validated: consensus-hom slope predicts the informative-marker slope at r=0.92, so the
   gate transfers reliably even though the exact magnitude does not -- see open question 1.)
4. Magnitude: do NOT apply the consensus slope 1:1 to informative markers (it over-predicts
   on clean, under-predicts on dirty). Once gated in, estimate the per-marker correction
   slope ON the informative donor-hom markers themselves (host-allele frac ~ a + b*n_carriers,
   weighted), pooled across the patient's serial timepoints to beat per-sample noise. Subtract
   b * n_carriers * depth (dose part only) from each donor-hom host-allele count.

This replaces the hardcoded magnitude with a per-run empirical estimate plus a statistical
significance test. The only remaining knob is the significance level (e.g. p<0.05) and a
minimum effect size worth correcting, both principled rather than tuned to one dataset. A
clean flowcell self-selects out; a dirty one self-calibrates its own correction. Division of
labour: consensus-hom = the gate and the run-level contamination level; informative-marker
regression = the per-marker correction magnitude.

This is a natural extension of the existing per-patient error table (#23): same joint-call
inputs, but adds the per-carrier contamination-dose model across the whole flowcell, and a
go/no-go significance gate.

## Open questions / validation before implementing

1. Marker-class transfer -- TESTED (`output/figure_review/step30_slope_transfer.py`).
   Per mixture, the per-carrier contamination slope from consensus-hom sites vs from the
   informative donor-hom markers (at the true-0% endpoint), carrier-COUNT dose:
   Pearson r = 0.92, median info/cons ratio = 1.00 (unbiased on average). Allele-copies
   dose was worse (r = 0.86), so carrier COUNT is the dose to use. CAVEAT: per-mixture
   ratios span 0.35-1.75 and the deviation is systematic, not just noise -- on CLEAN
   mixtures consensus OVER-predicts the informative slope (ratio ~0.4 -> over-correction),
   on CONTAMINATED mixtures it UNDER-predicts (ratio ~1.5-1.75 -> under-correction). This
   is exactly the pattern that made the naive consensus-hom subtraction over-correct clean
   mixtures in the prototype. Interpretation:
   - The GATE transfers well: r=0.92 means a significant consensus-hom slope reliably flags
     real informative-marker contamination, so the per-flowcell significance gate is sound.
   - The MAGNITUDE does not transfer 1:1: do NOT apply the raw consensus slope to the
     informative markers. Use consensus-hom for the GATE + run-level contamination level,
     and calibrate the actual per-marker correction magnitude on the informative markers
     themselves (within-sample regression, pooled across the patient's serial timepoints to
     beat per-sample noise; deployment-valid since it needs no endpoint).
   The FLOOR also differs by marker class (informative markers carry a ~0.15% mapping floor
   consensus-hom sites lack); subtracting only the dose part (slope * n) and leaving the
   floor to the error model handles that.
2. Deployment without co-pooled genotypes. In a cohort the carrier counts are known
   directly. In deployment the co-pooled genotypes may not be visible; the dose-response is
   still estimable from the run's own consensus-hom sites, and serial monitoring (>=3
   timepoints/patient) gives extra data to stabilize the per-patient slope.
3. Donor-het markers. The prototype corrects only donor-hom types 0/1. Decide whether
   donor-het markers need a contamination term too (smaller effect; near-0.5 balance).
4. Full validation: run the gated correction on all mixtures + the semi-synthetic ladder +
   genuinely clean high-depth controls, and confirm net improvement with zero harm on clean
   runs, before turning it on by default. Keep it behind a flag initially.

## Relationship to the cheap interim fix

The presence-test gate (report "not detected" when the host-presence test is negative)
already hides the cosmetic S12 floor at the report layer with no estimator change. Step 30
is the deeper fix that corrects the MLE number itself (the upward inflation at low
dilutions, not just the zero-host display). Ship the presence gate now if the floor display
is the immediate concern; pursue Step 30 (gated, per-flowcell) for sub-1% accuracy on
co-pooled panels.

## Implementation status (issue #30, behind a flag, off by default)

Implemented in `src/allomix/marker_contamination.py` and wired through the estimator:

- `ContaminationCorrection` (carrier counts + slope + gate decision) is carried on
  `PanelCalibration.contamination_correction` and applied in
  `estimate_single_donor_bb` before the MLE. When the field is None or the table gated
  out, `apply_contamination_correction` returns the marker list unchanged, so the default
  path is byte-identical (guarded by a unit test, like the #33 opt-out).
- `estimate_contamination_table` builds the table: carrier counts from the cohort
  genotypes, the per-flowcell gate from the consensus-hom dose response, and the
  per-patient magnitude from the informative donor-hom slope. Both slopes are fit by a
  shared-slope / per-timepoint-intercept weighted regression (`_grouped_slope_fit`), which
  pools a patient's serial timepoints without letting their differing host levels leak into
  the slope. This per-timepoint centering matters: pooling raw fractions under one intercept
  inflated the magnitude on moderately-contaminated pairs and caused over-correction.
- CLI: `allomix build-contamination-table` writes the table; `--contamination-table` plus
  `--contamination-correction` (BooleanOptionalAction, default OFF) apply it.

No-harm sweep reproduced with the package (`output/figure_review/step30_package_sweep.py`,
the prototype as oracle; SRP434573, all 10 mixtures):

- The contaminated mixtures are fixed: M3->F2 zero-host floor 0.165%->0.047% with dilution
  MAE 0.449%->0.394%; F3->F2 floor 0.101%->0.042%. 6/10 mixtures HELP, all 10 floors driven
  toward 0.
- The required clean pairs are preserved: F2->M1 and F2->M2 are exact no-ops (their
  informative slope is <=0 and clamps to 0); M1->M2 is within +/-0.008pp (noise level) with
  an improved floor. The clean pairs are preserved by the per-patient MAGNITUDE clamping to
  ~0, not by the gate: SRP434573 is one co-pooled (contaminated) flowcell, so the
  per-flowcell consensus gate correctly fires for every pair on it. A genuinely clean
  flowcell gives a flat consensus slope and gates out (unit-tested).

Open / deferred: donor-het markers are still not corrected (smaller, near-balanced effect);
multi-donor estimation does not apply the correction yet; the gate's `min_slope` effect-size
knob defaults to 0 (significance-only) since the magnitude self-clamps the low-contamination
pairs on this data.
