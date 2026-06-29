# SRP434573 figure-review follow-ups (from the #33 paper-review pass, 2026-06-24)

Three questions Dave raised reviewing the rebuilt two-rho paper figures. Recorded so
they are not lost. Q1 and Q3 have overnight experiments (see
`output/figure_review/` and the findings doc `srp434573_figure_review_findings.md`).

Context: facts/figures were rebuilt under the per-marker-type (two-rho) overdispersion
default. See [[issue33-two-rho-validation]] and the plan note "Per-marker-type
overdispersion adopted as default".

## Q1. Figure 4: the 0.5% known-host points read ~0.28%, about half nominal. Is the realised mix fraction actually below the nominal label, or is the estimator biased low?

Observation: the 3 real 0.5% timepoints (F2->M1, F2->M2, M1->M2) estimate at
MLE 0.26-0.29% and presence 0.278-0.292% (the two independent readouts agree at ~0.28%),
about half the 0.5% label. Both are still detected (presence p~0).

Claim to verify with data, not assertion: "realised mix fraction likely below its
nominal label." Need a model-free, marker-level measurement of the realised host
fraction and a decomposition vs the contamination/error floor.

How (implemented in `paper`-independent script, reuses the allomix API):
- Load per mixture: `genotype.parse_vcf(panel, sample=host/donor)` for host+donor GTs,
  `parse_vcf(admix, sample=<titration>)` for admix AD; `genotype.classify_markers(...)`.
- Realised host fraction = pooled host-allele read fraction at DONOR-HOMOZYGOUS
  fully-informative markers (marker_type 0: host hom-ref/donor hom-alt -> host allele =
  REF -> host reads = admix_ad_ref; marker_type 1: host hom-alt/donor hom-ref -> host
  reads = admix_ad_alt). Pool: sum(host reads)/sum(dp). No model, no overdispersion.
- Floor = minor-allele fraction at CONSENSUS-HOMOZYGOUS markers (host and donor hom for
  the same allele; the minor allele cannot come from either) = sequencing-error +
  co-pooled contamination background. (`contamination.py` / `_select_consensus_hom`.)
- Compare per mixture x titration: nominal vs raw-realised vs (raw - floor). Fit slope
  of raw-realised vs nominal across the ladder. If raw-realised tracks nominal with
  slope ~1 + offset ~floor -> estimator/floor story; if raw-realised itself sits below
  nominal -> the mixing is genuinely below label.
- Decision: whichever it is, state it from the numbers in results.md (and reconsider the
  "scattered between 0.26 and 0.29%" wording, which is a tight low-biased cluster, not
  scatter).

## Q2. Where are the BAM (real-reads) mixture figures? (S7 is purely in silico.)

Answer (documentation, no experiment needed): S7 (`fig_lod_saturation.png`) and the main
LoD sweep (Fig 1 / `fig5_lod_curves.png`) are PURELY IN SILICO (simulator draws reads
from a binomial; `run_lod_validation.py`). The real-reads / public-data-BAM-derived
results are:
- Figure 4 (`fig_srp434573.png`): the real SRP434573 two-person dilution series + the
  real three-person mix (real reads through the full genotyping path).
- Figure 5 (`fig_subsample_lod_grid.png`): the real-data LoD, sub-sampled reads+markers
  from the high-depth SRP434573 mixture BAMs (`run_subsample_lod.py`). This is the
  BAM-derived LoD counterpart of the in-silico Fig 1/S7.
- The semi-synthetic sub-0.5% ladder (`fig_srp434573_synthetic.png`): real BAMs
  sub-sampled and remixed at known fractions.
Possible paper action: a one-line cross-reference in the S7 caption pointing to Fig 5 as
the real-reads counterpart, so the in-silico vs real-reads split is explicit.

## Q3. The MLE returns a small positive floor (0.01-0.17%) at true-zero-host samples (S12 pure-donor endpoints), with a wide CI, while the presence test correctly says not-detected. This is ugly. Can we fix it, and if not, write it up.

Observation: at the 10 pure-donor endpoints (true host = 0%), MLE host = 0.00-0.17%
(F3->F2/F2 0.10, M3->F2/F2 0.17, F2->F1/F1 0.04, M3->F3/F3 0.02, others ~0), presence
test p=1.0 (not detected, markers all on the "0" row). The MLE is a magnitude estimator
bounded at >=0 with no "absent" state, so at zero host it returns a small positive
best-fit value.

How (overnight experiment):
- Profile the host-fraction log-likelihood near 0 for each endpoint (fine grid 0-1%)
  using the single-donor BB likelihood (`likelihood.total_log_likelihood_bb` /
  `chimerism.estimate_single_donor_bb`). Is there a genuine interior peak > 0
  (contamination leaking onto the host allele) or is the curve flat/monotone to 0
  (boundary artifact)? Determines whether the floor is signal or estimator boundary.
- Test contamination-floor subtraction: does passing the in-data contamination estimate
  (error-table / consensus-hom floor) pull the endpoint MLE toward 0?
- Quantify the CI behaviour at the boundary.

Candidate remedies to evaluate (pick from the data):
1. Gate the reported MLE on the presence test: when host-presence is not detected,
   report "< LoD" / "not detected" instead of a misleading positive point estimate.
2. Report a one-sided upper bound ("host < X%") at the boundary rather than a point.
3. Subtract the per-sample contamination/error floor from the host fraction before
   reporting, so the zero-host point estimate collapses to ~0.
4. A boundary likelihood-ratio test for host > 0 (this is essentially the presence test;
   the cleanest fix may be to surface the presence call alongside the MLE in the report).
If none cleanly fixes it, WRITE IT UP: the MLE point is "how much given host present",
the presence test is "is any host present", and the raw MLE point must not be read as a
detection at/near zero. Add the S12 caption sentence explaining the endpoints.

## Status / outputs
- Experiment scripts + data: `output/figure_review/` (gitignored).
- Findings: `claude/srp434573_figure_review_findings.md` (written by the overnight run).
- No paper prose committed from these without Dave's review (the discussion CI-emphasis
  point from the earlier CI check is also pending his call).
