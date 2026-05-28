# Plan: Host-Presence Detection at Donor-Homozygous Markers

Status: proposed, not implemented.

## TL;DR

Add a dedicated detection test for "is the host (minor contributor) present at all?",
separate from the fraction MLE. The test uses only the markers where the donor is
homozygous and the host carries the allele the donor lacks. At those markers the
donor-absent allele should appear at essentially zero frequency in a pure-donor
sample (sequencing error only), so any reads of it are evidence of host. The
statistic is a one-sided count test (and a likelihood-ratio test giving a host-fraction
estimate) against a per-site error background, combining evidence across markers.

This is the same read counts the MLE already sees, but reframed as a detection
question and freed from two things that blunt the MLE at very low fractions: the
single shared overdispersion parameter and the global error rate. It is reported
alongside the MLE, not folded into it (route A); a unified two-component likelihood
is sketched as a follow-up (route B).

**Hard dependency:** Step 14 (empirical per-site error rates) gives the per-site,
per-direction background this test needs. Without it the test falls back to the
global `--error-rate`, which caps the achievable detection limit at the error floor
and removes most of the benefit. Step 14 should land first.

**Soft dependency:** Step 12 (per-marker context refactor) makes the plumbing cleaner
if it lands first, but this plan does not strictly require it because the detector is
a separate code path from the estimator.

## Background

### What the current pipeline does, and why low-fraction detection is under-served

allomix estimates a donor fraction `f` by maximum likelihood over a beta-binomial with
a single shared concentration `rho` (`src/allomix/chimerism.py:595`), using a 4-state
error model with one global `error_rate` (`src/allomix/chimerism.py:144`). Detection is
derived from estimation: `detection_limit()` computes EP17 LoB/LoD from the Fisher
information of the fit (`src/allomix/chimerism.py:550`), and the LoD-validation rule is
"estimate exceeds LoB" (`claude/issue_8_lod_markers_plan.md`). There is no standalone
"is the minor contributor present?" hypothesis test.

The markers this plan targets are already in the model. A site where the donor is
homozygous and the host carries the donor-absent allele is exactly a Vynck type 0 or 1
(full contrast) or type 10 or 11 (host heterozygous). So the reads are not new data.
But three things mean a dedicated detector extracts more signal at very low host
fractions, and the third is a reason to have it regardless:

1. **The shared `rho` taxes exactly these markers.** Overdispersion is real at
   het/intermediate markers (amplification bias). At a donor-homozygous marker the
   donor-absent allele is essentially absent, so there is nothing to be biased and the
   only noise is sequencing error (near-Poisson). The beta-binomial inflates every
   marker's variance by `1 + (n-1)/(rho+1)` (`src/allomix/chimerism.py:537`). At 1000x
   with a fitted `rho ~ 100` that is ~11x variance inflation (SE up ~3.3x). Because
   `detection_limit()` uses the pooled `rho_mle` (`src/allomix/chimerism.py:713`), the
   het markers' overdispersion is charged against the clean near-zero markers,
   degrading the reported LoD by roughly that factor. A detector restricted to the
   donor-homozygous markers and modelled as error-background-only does not pay that tax.

2. **The global error rate is the binding constraint here and is far too high.**
   Detection power at these markers is set by the donor-absent-allele error background.
   The global rate gives an effective per-direction floor of `e/3 ~ 3.3e-3`, one to two
   orders of magnitude above a clean SNP's real per-site rate (~1e-4 to 5e-4).
   Worked example at 1000x with ~20 donor-homozygous markers: with the global rate the
   background is `sum(n_i e_i) ~ 20*1000*0.0033 ~ 66` stray reads, which swamps a 0.1%
   host signal (~20 reads). With realistic per-site rates the background drops to ~10
   reads and 0.1% becomes clearly detectable. This is exactly why Step 14 is a
   prerequisite.

3. **Detection is a different question from estimation.** For relapse monitoring the
   clinical question is "is host coming back?", a one-sided test (H0: f_host = 0)
   against a near-zero background. The MLE's "estimate crossed LoB" is an indirect
   proxy. A direct p-value is more powerful at the low end and more interpretable, and
   is worth reporting even if the likelihood were perfectly specified.

## Marker set

Use markers where the **donor is homozygous** and the host carries the donor-absent
allele. In Vynck types (host genotype, donor genotype):

| Type | Host | Donor | Donor-absent allele | Host dose of it |
|:----:|:----:|:-----:|:--------------------|:---------------:|
| 1    | 1/1  | 0/0   | ALT                 | 2               |
| 0    | 0/0  | 1/1   | REF                 | 2               |
| 10   | 0/1  | 0/0   | ALT                 | 1               |
| 11   | 0/1  | 1/1   | REF                 | 1               |

Exclude types 20 and 21: there the donor is heterozygous and carries both alleles, so
there is no donor-absent allele to count against a clean background.

Full-contrast types (0, 1) carry twice the per-read signal of host-het types (10, 11)
and are the strongest. Both are usable.

**Multi-donor:** the donor-absent allele must be absent from *every* donor, not just one,
otherwise a low-level second donor carrying that allele would trip the test. Restrict to
markers where all donor genotypes are homozygous for the same allele and the host carries
the other.

## The statistic

For each selected marker `i`:

- `y_i` = observed read count of the donor-absent allele, `n_i` = depth.
- `e_i` = per-site, per-direction background rate of observing the donor-absent allele
  under a pure-donor (homozygous) background. This is `e_refalt` when the donor is
  hom-ref (absent allele is ALT) or `e_altref` when the donor is hom-alt (absent allele
  is REF). These are exactly the two directions produced by Step 14
  (`claude/14_empirical_error_rates_plan.md`). Apply the per-site floor from that plan.
- `h_i in {1, 2}` = host dose of the donor-absent allele.

Expected donor-absent-allele frequency given host fraction `f_h` (small-rate
approximation):

```
q_i(f_h) ~= e_i + (h_i / 2) * f_h
```

**Null (no host present), f_h = 0:** `y_i ~ Binomial(n_i, e_i)`.

**Pooled count test (primary, transparent):**

```
Y    = sum_i y_i
Lam  = sum_i n_i * e_i           # expected stray-read count under H0
p    = P(Poisson(Lam) >= Y)      # one-sided
```

Poisson is a good approximation because `e_i` is tiny and `n_i` large; use an exact
Binomial-sum or a saddlepoint correction if `Lam` is small and accuracy at the tail
matters.

**Likelihood-ratio test (gives an estimate):**

```
loglik(f_h) = sum_i [ y_i*log(q_i(f_h)) + (n_i - y_i)*log(1 - q_i(f_h)) ]
f_h_hat     = argmax_{f_h >= 0} loglik(f_h)
D           = 2 * (loglik(f_h_hat) - loglik(0))
```

The parameter is bounded at 0, so under H0 the statistic follows a chi-bar-square
mixture `0.5*chi2_0 + 0.5*chi2_1`; the one-sided p-value is `0.5 * P(chi2_1 >= D)` for
`D > 0` (and 1 when `f_h_hat = 0`). Report `f_h_hat` and a profile-likelihood CI.

The LRT is slightly more powerful when depths and host doses vary across markers and
yields a host-fraction estimate specialised to the low-fraction regime; the pooled
Poisson test is a robust cross-check. Report both.

**Cross-check against the global MLE:** `f_h_hat` should agree with `1 - donor_fraction`
from `estimate_single_donor_bb`. A material disagreement (presence test says host is
there, MLE says ~100% donor, or vice versa) is itself a QC signal worth surfacing.

### Expected sensitivity

LoD of this detector is governed by `Lam`: a host signal is detectable when
`sum_i n_i (h_i/2) f_h` is a few times `sqrt(Lam)` above background. With ~20 clean
donor-homozygous markers (`e_i ~ 5e-4`) at 1000x, `Lam ~ 10`, so f_h around 0.05-0.1%
is detectable, below what the global-error MLE achieves. The lever is the per-site error
rate, which is why Step 14 gates the headline number.

## Fix scope

New module plus reporting and CLI plumbing. The detector is a separate path from the
estimator, so it does not touch the MLE likelihood functions.

### 1. `src/allomix/detect.py` (NEW)

```python
"""Host-presence detection at donor-homozygous markers.

Tests whether a minor contributor (the host, re-occurring post-HSCT) is present
at all, using the markers where the donor is homozygous and the host carries the
donor-absent allele. At those markers the donor-absent allele is expected at the
sequencing-error background in a pure-donor sample, so its read counts give a
one-sided detection test against that background.

This is complementary to the fraction MLE in ``chimerism``: the MLE estimates the
magnitude, this test guards the low end and answers "is host present?" directly.
"""
```

- `HostPresenceResult` dataclass: `n_markers`, `n_donor_absent_reads` (Y), `expected_background` (Lam), `poisson_pval`, `lrt_pval`, `f_host_mle`, `f_host_ci`, `used_per_site_error` (bool), `error_rate_source` (str: "per-site" | "global-fallback").
- `select_donor_hom_markers(markers) -> list[...]`: filter to types 0/1/10/11 (multi-donor: absent from all donors), returning per-marker `(y_i, n_i, h_i, direction)` where `direction` selects `e_refalt` vs `e_altref`.
- `host_presence_test(markers, marker_errors=None, error_rate=0.01, error_floor=1e-5) -> HostPresenceResult`: assemble `e_i` from the Step 14 table (per direction) with a floor; fall back to `error_rate/3` per direction when a site is missing, setting `error_rate_source` accordingly. Compute pooled Poisson p-value and the LRT (MLE + profile CI + chi-bar-square p-value).

Reuse the Step 14 error-table loader (`allomix.error_rates.load_error_table`) rather than
re-reading the TSV.

### 2. `src/allomix/cli.py`

- `monitor` and `timeline` already accept `--error-table` once Step 14 lands. Run the
  presence test by default in `_run_single_sample` (it is cheap) and attach the
  `HostPresenceResult` to the per-sample output. Add `--no-host-presence` to disable.
- When no error table is supplied, run with the global fallback and mark
  `error_rate_source="global-fallback"` so the report shows the LoD is error-floor-limited.

### 3. `src/allomix/report.py`

- Add columns to the TSV summary: `host_present_p` (the LRT p-value), `host_f_est`,
  `host_f_ci_lo`, `host_f_ci_hi`, `host_detect_markers` (count), `host_err_source`.
  Append, do not rename existing columns (downstream `scripts/run_xls_batch.py` and paper
  figure scripts parse the TSV by header).
- Mirror the same fields in the JSON output under a `host_presence` object.

### 4. `src/allomix/qc.py`

- Optional: when `host_present_p < 0.01` but the MLE donor fraction is ~1.0 (or the
  presence estimate disagrees materially with `1 - donor_fraction`), add a REVIEW-level
  warning ("low-level host signal detected below the fraction estimate's resolution" or
  "presence test and fraction estimate disagree"). This is the clinically useful output
  for relapse monitoring. Keep it a warning, not a hard status change, in v1.

### 5. Tests

- `tests/test_detect.py` (NEW):
  - `select_donor_hom_markers` picks types 0/1/10/11 and rejects 20/21; multi-donor case
    requires absence from all donors.
  - Pure-donor synthetic input (donor-absent reads at the background rate only) gives a
    non-significant p-value at the expected false-positive rate.
  - Spiked input (a few donor-absent reads above background at several markers) gives a
    significant p-value and an `f_host_mle` close to the spiked fraction.
  - Per-site error table changes the background and shifts the p-value in the expected
    direction; missing-site fallback to global rate sets `error_rate_source` correctly.
  - LRT MLE agrees with `1 - donor_fraction` from `estimate_single_donor_bb` on a shared
    synthetic sample.
- `tests/test_cli.py`: smoke test that `monitor` emits the new columns and that
  `--no-host-presence` suppresses them.

## Validation plan

1. **Unit sanity** (cheapest): `pytest tests/test_detect.py -x -q`.
2. **Synthetic LoD comparison**: extend the Step 8 / issue #8 LoD sweep
   (`paper/scripts/run_lod_validation.py`, `claude/issue_8_lod_markers_plan.md`) to also
   record, per replicate, the presence-test p-value and whether it called detection at
   alpha=0.05. Compute a presence-test LoD (lowest true fraction with >=95% detection)
   alongside the existing MLE-LoB LoD. Run with a per-site error table provided (perfect
   rates, matching how the issue #8 sweep provides perfect biases) to show the achievable
   gain, and again with the global rate to show the floor effect. Expectation: at low
   fractions and high depth the presence test detects below the MLE-LoB LoD when per-site
   errors are available; the two converge when only the global rate is used.
3. **Real-data smoke test**: run on the post-HSCT samples in
   `output/validation_run_new_bias2/` with the per-site error table from the
   bias-training cohort. Check that full-donor samples give non-significant presence
   p-values (no false relapse calls) and that any sample with known low-level host shows
   a significant p-value. Record concordance with the MLE fraction.

## Edge cases and risks

- **Index hopping / barcode bleed (assay-level, the dominant false positive).** The host
  genotyping library carries the donor-absent allele at high VAF, so hopped host reads
  mimic host re-occurrence. This is a wet-lab requirement, not a code fix: sequence host
  and admixture on separate runs, use unique dual indexes (UDIs), keep any library
  carrying those alleles (siblings, the patient's earlier samples) off the admixture run,
  and include a no-template/unrelated control to measure the hopping rate. Document this
  in the user guide as a precondition for trusting low-level detection. Note that physical
  run separation is orthogonal to bioinformatic joint genotyping (`doc/joint_calling.md`):
  samples can still be joint-called across separate runs, so this does not break the
  2-element-AD requirement.
- **Donor heterozygote miscalled as homozygous.** If a donor het is called hom (shallow
  donor coverage), the allele is genuinely present and every read is a false positive.
  For the selected markers require stricter donor evidence than the default GQ>=20:
  high donor GQ and a donor-sample VAF at that site consistent with homozygous
  (the donor's own AD for the absent allele should be ~0, not borderline). Deep donor
  genotyping pays off directly here.
- **Per-site error overconfidence.** A site with zero observed errors in training would
  give `e_i = 0` and make a single stray read infinitely significant. Apply the per-site
  floor from Step 14 (`error_floor`, e.g. 1e-5).
- **Sites missing from the error table.** Fall back to the global rate per direction and
  flag `error_rate_source="global-fallback"`. Do not silently treat missing as zero.
- **Double-counting with the MLE.** Route A reports the two side by side; it does not sum
  the same reads into one combined likelihood, so there is no double counting. The unified
  model (route B) would replace, not add to, the current likelihood.
- **Strand bias and base quality.** Strand-skewed artifacts (deamination, OxoG) and
  low-BQ reads are the residual false-positive sources after per-site rates. A strand-bias
  filter needs F1R2/F2R1 or SB counts; per-marker BQ is Step 17 (`claude/17_bq_aware_plan.md`).
  Both are optional sharpeners, out of scope for v1.
- **Mapping artifacts / paralogy.** Handled at panel design (avoid low-mappability sites);
  the per-site error table also down-ranks sites with chronically elevated background.

## Out of scope (follow-ups)

- **Route B: unified two-component likelihood.** Give the per-marker likelihood a noise
  model that uses tight error-background variance at near-zero markers and the
  overdispersion `rho` at intermediate markers (a per-marker or two-class `rho`, or an
  explicit error-only component). This folds the detection signal into the single MLE
  properly and the LRT-against-f=0 becomes native. Larger change; depends on Step 12 and
  Step 14. Track separately once route A is validated.
- **Second-donor / generic minor-contributor detection.** The same construction detects a
  low-level second donor against a (host + donor1) homozygous background. Natural
  generalisation; keep v1 scoped to host presence.
- **Strand-bias filter** (needs SB / F1R2 FORMAT fields).
- **Per-marker BQ weighting of the background** (Step 17).

## File-by-file checklist

- [ ] `src/allomix/detect.py` (NEW): `HostPresenceResult`, `select_donor_hom_markers`,
      `host_presence_test` (pooled Poisson + LRT + profile CI), error-table assembly with
      per-direction fallback and floor.
- [ ] `src/allomix/cli.py`: run the presence test by default in `_run_single_sample`,
      attach result; add `--no-host-presence`. Reuse `--error-table` from Step 14.
- [ ] `src/allomix/report.py`: add `host_present_p`, `host_f_est`, `host_f_ci_lo`,
      `host_f_ci_hi`, `host_detect_markers`, `host_err_source` to TSV summary; mirror in
      JSON under `host_presence`.
- [ ] `src/allomix/qc.py`: optional REVIEW warning when presence is significant but the
      MLE fraction does not reflect it, or the two disagree.
- [ ] `tests/test_detect.py` (NEW): marker selection, false-positive rate on pure donor,
      power on spiked input, error-table effect and fallback, agreement with MLE.
- [ ] `tests/test_cli.py`: new columns present; `--no-host-presence` suppresses them.
- [ ] `paper/scripts/run_lod_validation.py`: record presence-test detection alongside the
      MLE-LoB rule; compute presence-test LoD with and without per-site errors.
- [ ] `claude/allomix_overall_plan.md`: add Step 20 entry linking here (done).

## Dependencies summary

- **Step 14 (empirical per-site error rates)** — `claude/14_empirical_error_rates_plan.md`.
  Hard prerequisite for the headline LoD. Land first.
- **Step 12 (per-marker context refactor)** — `claude/12_marker_context_refactor_plan.md`.
  Soft; only matters for route B.
- **Step 17 (per-base quality)** — `claude/17_bq_aware_plan.md`. Optional sharpener of the
  background, not required.
- **issue #8 LoD sweep** — `claude/issue_8_lod_markers_plan.md`. Reused to validate the
  detection-LoD gain.
