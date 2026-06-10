# Further improvements from the SRP434573 run

What the real-data run on SRP434573 (issue #16) suggests we could still improve in
detection and QC. SRP434573 is a HiSeq 3000 MIP panel (~1050 sites, single-strand
pear-merged reads) of known two- and three-person mixtures of 7 individuals,
titrated from 10% down to 0.5% minor fraction, with each contributor reused across
many runs. The minor (titrated) contributor is mapped to HOST, so the dilution
series is a relapse / declining-chimerism series.

The items below are grounded in numbers from that run (cited inline). Confidence is
marked: **[data]** when the observation is directly measured here, **[likely]** when
it is a well-supported inference, **[speculative]** when it needs its own check.

All measurements come from `allomix monitor` on the per-patient VCFs in
`output/genotypes/SRP434573/` and from `scripts/probe_*_srp434573.py`.

## Status (2026-06)

**Wave 1 (low-end MLE accuracy) is implemented and validated in silico.**

- **One-sided robust trim (Observation 1): DONE.** `chimerism.ROBUST_ONE_SIDED`.
  A marker whose residual deviates toward host presence is no longer trimmed from
  the fit. In-silico validation against the collapse this fixes (ground truth
  known, `scripts/validate_onesided_trim.py`): the symmetric trim collapsed the
  host estimate to ~0 in 47-80% of 0.5% replicates in noisy / low-depth regimes;
  one-sided collapses in 0% and recovers the true fraction (e.g. 0.51% at a true
  0.5%, tracking the presence estimate). MAE over non-blank fractions fell 30-42%
  in those regimes. At deployment depth (>1000x) it is a measured no-op
  (MAE 0.044% -> 0.043%, identical blank). **[data]**
- **Pre-trim goodness-of-fit (Observation 5): DONE.** `qc.goodness_of_fit_pval_pretrim`.
  GoF is now computed on the full marker set as well as the post-trim set, the
  REVIEW gate uses the worse of the two, and a distinct warning fires when a trim
  is masking a poor full-set fit.
- **Presence-estimate backstop (Observation 1, first improvement): deferred.**
  The trim fix removes the MLE collapse this was meant to catch (the rescued MLE
  now matches the presence estimate), so silently swapping the headline number is
  redundant and carries reporting risk. Revisit only if real SRP434573 / wetlab
  data shows residual collapse the trim does not catch.

**One carry-forward from wave 1:** the one-sided trim raises the blank estimate by
~0.05-0.07% in the noisiest low-depth regimes (it keeps host-direction noise that
the symmetric cut removed). This is small and absent at deployment depth, but it is
a second noise term that the LoD floor should account for, which feeds directly
into Observation 2. **[data]**

**Next: Observation 2** (contamination floor into the LoD and the host-presence
background).

**Cohort / multi-sample phase.** Observations 3, 4, and the pooled-bias-table part
of Observation 7 all block on a multi-sample entry point that `monitor` (single
patient) and `estimate-bias` (single pair) do not provide. They are consolidated
as a forward plan in `cohort_phase_plan.md`.

## Observation 1 (the big one): the robust-refit MLE collapses toward zero below ~1.5% host, while the presence-test estimate stays accurate

**Status: fixed by the one-sided trim (see Status above). The evidence that motivated it is kept below.**

Per-sample, two representative mixtures (known = true host %, mle = `100 - donor_pct`,
presF = host-presence `f_host_mle`, drop = robust markers excluded, gof = goodness of fit):

mix_M3_into_F3 (per-marker contamination floor ~0.18-0.31%):

| known % | mle % | presF % | robust drop | gof | qc |
|--:|--:|--:|--:|--:|--|
| 1.0 | 0.68 | 1.11 | 200 (35%) | 2e-4 | REVIEW |
| 1.25 | 1.23 | 1.44 | 160 (28%) | 7e-4 | REVIEW |
| 2.5 | 2.58 | 3.04 | 91 (16%) | 0.32 | REVIEW |
| 5.0 | 4.44 | 4.62 | 46 (8%) | 0.63 | REVIEW |
| 10.0 | 10.98 | 11.08 | 28 (5%) | 0.54 | PASS |

mix_F2_into_M1 (cleaner library, contamination floor ~0.06-0.09%):

| known % | mle % | presF % | robust drop | gof | qc |
|--:|--:|--:|--:|--:|--|
| 0.5 | 0.07 | 0.13 | 179 (31%) | 0.0 | REVIEW |
| 1.0 | 0.00 | 0.74 | 292 (50%) | 4e-4 | REVIEW |
| 1.25 | 0.00 | 0.81 | 334 (57%) | 1.0 | REVIEW |
| 2.5 | 1.40 | 1.88 | 106 (18%) | 0.26 | REVIEW |
| 5.0 | 3.66 | 4.30 | 71 (12%) | 0.96 | PASS |
| 10.0 | 8.41 | 8.70 | 0 | 1.0 | PASS |

The diagnosis (confirmed by the fix): the robust refit dropped a runaway fraction of
markers as the host fraction fell, and the markers it trimmed were the host-carrying
ones (at low host fraction their VAF deviation from the donor-dominated fit is small
and the symmetric median/MAD cut reads it as an outlier), so trimming removed the
very signal we want and the estimate collapsed toward the donor-only solution. The
one-sided trim protects residuals that deviate toward host presence, which removes
the collapse without touching clean-data or deployment-depth behaviour. **[data]**

### Still open under Observation 1

- **Recalibrate the robust-drop REVIEW so it discriminates at the low end.** With
  the one-sided trim the drop fractions are much smaller, so this is less acute, but
  a bare >15% drop count still does not separate "expected low-fraction trimming"
  from "real CNV/LoH/genotyping problem". Tie the threshold to the host fraction, or
  report the dominant *reason* for the drops (direction of trimmed residuals,
  clustering on a chromosome arm) instead of a count. Re-measure the post-fix drop
  distribution before changing the threshold. **[speculative]**

## Observation 2 (NEXT): a ~0.2% co-pooled contamination floor competes directly with sub-1% detection, and it is independent of the host-presence background

The contamination probe established a per-marker floor of ~0.2% (median; up to ~1.5%
at p95) from co-loaded genomes, dose-proportional to how many co-pooled samples carry
the allele, flat across host fraction, and varying by library (M3-F3 ~0.18-0.31%,
F2-M1 ~0.06-0.09%) **[data]**. The host-presence detector tests donor-absent reads
against a *sequencing-error* background only; it does not know about this
contamination floor, which lands on exactly the same donor-absent alleles wherever a
co-pooled genome carries the host's allele **[likely]**.

### Improvements

- **Feed the contamination estimate into the host-presence background and the
  reported LoD.** The presence test and `contamination.py` are currently independent.
  The contamination floor is a second noise term that should (a) raise the
  per-marker background the presence test compares against, and (b) floor the
  reported `lob_fraction` / `lod_fraction`, which today are computed from sequencing
  error and Fisher information alone. At a 0.2% floor, a 0.5% host call is within
  noise and the LoD should say so. The one-sided trim's small blank inflation
  (~0.05-0.07% at low depth, see Status) is a third noise term to fold into the same
  floor. This is the natural next step now that all three quantities exist.
  **[data-backed]**
- **Per-marker contamination flag for the presence test**, not just a per-sample
  scalar: a donor-homozygous presence marker where a co-pooled genome carries the
  donor-absent allele is the one whose count is inflated. With the cohort genotypes
  we showed this is identifiable; in deployment we cannot see co-pooled genotypes,
  but the dose-response means the floor can still be apportioned. **[speculative]**

## Observation 3 -> cohort plan

Batch / run-level contamination QC (group admix samples by `##allomixRunUnit` and
flag a whole flowcell lane when its samples share an elevated floor). Needs the
multi-sample entry point. See `cohort_phase_plan.md`. **[likely]**

## Observation 4 -> cohort plan

Cohort-recurrence bad-site detection (a panel-level blacklist for loci that are
systematically inconsistent across the repeated cohort). Needs the multi-sample
entry point. See `cohort_phase_plan.md`. **[likely]**

## Observation 5: goodness-of-fit and robust-drop disagree at the low end

**Status: addressed by pre-trim GoF (see Status above).**

GoF was poor for M3-F3 at 1-1.25% (p ~2e-4 to 7e-4) but 1.0 for F2-M1 1_79 despite a
57% robust drop **[data]**, because GoF was computed on the post-trim marker set, so
a fit that trimmed away its problems reported a clean GoF. GoF is now also computed
pre-trim and the REVIEW gate uses the worse of the two, so a fit cannot look clean by
virtue of having discarded the inconvenient markers.

## Observation 6: the contamination flag thresholds are calibrated on one dataset

`CONTAMINATION_WARN_FRACTION` (0.2%) and `CONTAMINATION_REVIEW_FRACTION` (1%), and the
empirical p10 floor, were set from SRP434573 alone. The p10-floor estimator assumes a
useful fraction of no-carrier sites and a particular allele-frequency spectrum. This
gates Observation 2: any LoD floor built on the contamination estimate inherits these
thresholds.

### Improvement

- **Characterize the contamination flag on clean controls and other panels.** Measure
  the false-positive rate on genuinely uncontaminated high-depth samples, and confirm
  the empirical floor behaves on panels whose marker allele frequencies differ from
  this one. This is validation work, not a code change, but it gates trusting the
  thresholds in deployment, and it should be done alongside Observation 2. **[speculative]**

## Observation 7: per-marker amplification bias is real and reproducible but is not the limiting factor here; correcting it (even correctly) does not help the estimate

The additive bias model overcorrected at extreme expected VAF and was replaced with a
multiplicative (logit-space) correction (issue #20). On `mix_M3_into_F3` the old
additive form, applied from a panel bias table, collapsed goodness-of-fit from 0.54 to
0.0 at the 10% point; the logit form leaves it essentially intact (0.54 -> 0.13)
**[data]**. So the model is now correct. The surprise is that, once correct, bias
correction still does not move the result on this data.

Measured on `mix_M3_into_F3` (`allomix monitor` baseline vs `--bias-table`, with a
pooled caller-consistent both-het table built by `estimate-bias --both-het`):

| | MAE (host %) | gof @ 1% / 1.25% | REVIEWs |
|--|--:|--:|--:|
| no correction | 0.392 | 0.00 / 0.00 | 4/5 |
| pooled per-marker table (66% informative coverage) | 0.396 | 0.00 / 0.00 | 4/5 |

Three measured reasons it does not help **[data]**:

1. **The global mean bias is ~0.** Signed mean per-marker bias is +0.30% (both-het,
   mpileup) and -0.08% (panel-het, GATK); medians -0.30% / 0.00%. There is no
   systematic REF/ALT skew on this panel to subtract.
2. **That spread is already absorbed by overdispersion.** A ~6% per-marker VAF spread
   is extra-binomial variance, which the jointly fit beta-binomial `rho` already turns
   into wider CIs.
3. **The markers you can measure are not the markers you need.** Per-marker bias is
   only measurable where a marker is heterozygous; the fully-informative markers at low
   host fraction (host hom-ref / donor hom-alt) are homozygous in both contributors, so
   their bias cannot be measured from that pair at all. Only a table pooled across a
   cohort reaches them.

This says the levers for the low-end problems were the robust refit (now fixed) and
the contamination floor (Observation 2), **not** per-marker bias.

### Improvements

- **Deprioritize per-marker bias correction as a low-fraction accuracy lever.** It is
  worth keeping (small interior gain, harmless now that it is logit-space), but the
  evidence here is that effort is better spent on Observation 2. **[data-backed]**
- **Warn on a caller mismatch between the panel and admix VCFs.** A bias table
  estimated from GATK-called panel het sites and applied to `bcftools mpileup` admix
  data makes results worse, because per-marker bias is caller-specific (issue #11). A
  cheap guard: detect the likely source of each VCF (`DP4`/`I16` present = mpileup,
  absent with GATK annotations = GATK) and warn on a mismatch. This is a small,
  self-contained change, independent of the cohort work. **[likely]**
- **Give `estimate-bias --both-het` a cohort entry point** -> cohort plan. A pair's
  both-het markers only help *other* pairs, so the pooled table needs the multi-sample
  entry point. See `cohort_phase_plan.md`. **[likely]**

## Lower-priority / needs-its-own-investigation

- **Input-quality QC predicting poor fits.** The v1 library (F2-M1) fit worse than v2
  (M3-F3) at matched fractions. A pre-estimate input-quality check (depth uniformity
  across the panel, per-site dropout) might predict which samples will fit poorly and
  warrant manual review before the estimate is even attempted. **[speculative]**
- **chrX recovery for sex-matched pairs.** The panel carries ~27 chrX amplicons that
  we default-drop. For a host/donor pair of the same sex they are usable informative
  markers; inferring sex from the data and keeping chrX when it matches would add
  markers, which matters most at low fractions where every informative marker counts.
  **[speculative]**
