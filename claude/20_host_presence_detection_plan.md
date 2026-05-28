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

## Update (2026-05-28): overdispersion is now in the simulator — what it means here

A separate LoD investigation (Step 21) found that the in-silico LoD was optimistic
because the simulator drew reads from a pure binomial; it now confirmed empirically
that the per-sample LoD is dominated by overdispersion `rho` (the variance approaches
`p(1-p)/(rho+1)`, so depth saturates near an effective cap of `rho+1` reads). The
simulator gained a `rho` argument (`sample_allele_counts` / `blend_vcfs` /
`blend_from_genotype_dicts`, default `inf` = binomial). This affects this plan in three
ways:

1. **It is direct evidence for Background point 1** (the shared `rho` taxes exactly the
   markers this detector uses). The MLE's `detection_limit` pays the overdispersion tax;
   a detector restricted to donor-homozygous, error-background-only markers does not. The
   new figures `output/facts/fig_lod_saturation.png` and `fig_overdispersion_lod.png`
   (from `paper/scripts/plot_lod_saturation.py` and `run_overdispersion_lod.py`) quantify
   the tax: at `rho≈100`, LoD is ~3.5x the binomial value.

2. **The prototype below is unaffected and still correct as written.** It simulates with
   the binomial default (`rho=inf`) and the `e_i = e/3` self-consistency trick holds. Do
   **not** turn on a global `rho` in the prototype: the simulator applies one `rho` to
   *every* marker/allele uniformly, including the near-zero donor-absent allele, where
   overdispersion is not physical (it is a het/intermediate-marker amplification effect).
   A uniform `rho` would overdisperse the clean donor-absent background and miscalibrate
   the presence-test binomial null (inflated false positives), making the detector look
   bad for the wrong reason.

3. **Caveat for acceptance gate #3.** Because the prototype's binomial simulation gives
   the MLE no overdispersion to pay, the MLE-LoB LoD it measures is the optimistic
   binomial value and the presence detector's structural advantage from point 1 is
   invisible. The gate-#3 gap measured on binomial sim is therefore a **lower bound** on
   the real-data advantage — do not reject the detector if that gap is small; the gap is
   expected to widen once overdispersion is present. Fully demonstrating point 1 needs a
   **marker-type-aware** overdispersion model (apply `rho` only to het/intermediate
   markers, keep the donor-absent allele at the binomial error background), which is a
   Step 21 follow-up, not part of this prototype.

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

## Control data and calibration (build this first)

Before wiring the detector into the CLI and reports, generate dedicated positive and
negative controls and confirm the statistic is calibrated. This is the same EP17 LoB/LoD
construction already used for the MLE (`claude/issue_8_lod_markers_plan.md`), applied to
the presence statistic, and follows the protocol of generating extremely low-fraction
positives and running many replicates until detection stabilises. We do not currently have
this data: the existing LoD sweep
(`paper/scripts/run_lod_validation.py:75`) stops at a lowest positive fraction of 0.001
(0.1%) and uses a single global error rate, neither of which exercises the regime this
detector targets.

**Negative controls (the LoB / false-positive floor).** Many pure-donor replicates
(`f = 0`): donor-homozygous markers with simulated sequencing error only and no host. Each
yields a presence p-value and the donor-absent read count `Y`. Under H0 the p-values should
be approximately uniform; the 5% false-positive point (or the 95th percentile of `Y/Lam`)
sets the operating threshold. Run enough replicates to pin the tail (order 1,000), since the
whole claim is about rare false positives.

**Positive controls (the LoD).** Extremely low fractions below the current floor, e.g.
`f in {1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3}`, extended downward until the detection
rate hits zero, at deep coverage (1000x, 2000x, and a 5000x cell to mimic deep relapse
monitoring), with many replicates (>=60, more at the bottom because host reads are rare).
"Keep running replicates until a signal" operationalises as: measure detection rate at each
fraction and depth; LoD is the lowest fraction with >=95% detection above the LoB threshold.

**Generation.** Reuse `simulate.generate_related_genotypes`, `blend_vcfs`, and
`sample_allele_counts`. Nothing in the simulator blocks `f < 0.001`: at `f = 1e-4`,
2000x, a full-contrast marker has an expected 0.2 host reads, so detection relies on
pooling across many markers and replicates, which is exactly why the protocol runs many
replicates. No new simulator capability is needed for the first (global-rate) pass.

**Calibration pass (cheap, gate before CLI/report work).** A standalone script
(`paper/scripts/run_presence_lod.py`, or a prototype next to the issue #8 sweep) that
generates the controls, runs `host_presence_test`, and reports: (a) negative-control
false-positive rate vs nominal alpha, (b) the positive-control detection-rate curve and the
fitted 95% LoD, reusing the `fit_lod` / `compute_lob` helpers from
`paper/scripts/run_lod_validation.py`. Acceptance gate before building the reporting and CLI
integration: false-positive rate ~= alpha on negatives (calibrated), detection rate monotone
in `f` and depth, and presence-test LoD below the global-error MLE LoB at matched cells.

**Error-model consistency (the one real gap).** `sample_allele_counts` injects a single
global, symmetric error rate (`src/allomix/simulate.py:244`, scalar `error_rate`). For the
first calibration pass this is fine and self-consistent: simulate with rate `e`, give the
detector rate `e`, and check calibration and power against a known truth. No simulator change
needed. To validate the *per-site* advantage (the headline of this plan), the simulator must
inject per-site, per-direction error matching an error table, which is a small extension to
let `blend_vcfs` / `sample_allele_counts` take a per-marker (per-direction) rate instead of a
scalar. This is the same simulator update that Step 14 lists as an out-of-scope follow-up
(`claude/14_empirical_error_rates_plan.md`, "Out of scope"). Sequence: global-rate
calibration first (no sim change), per-site validation after that extension lands.

**Quality scores are not needed for any of this.** The detector uses AD counts plus per-site
error rates (Step 14), not per-read base qualities (Step 17). For synthetic controls the
error rate is declared, not derived from a quality score; for real-data negative controls the
empirical per-site error table captures the background directly. Step 17 (FORMAT/QS) remains
an optional sharpener only and is itself flagged "skeptical, may not ship" in the master plan,
so this validation and the detector can both proceed without the QS pipeline work.

### Prototype spec (the actual first task, self-contained)

A throwaway experiment script, no `src/allomix/` changes. The statistic is small enough to
inline here; promote it to `src/allomix/detect.py` only after the gate passes.

**File:** `paper/scripts/run_presence_lod.py` (paper-tree, sits next to the issue #8 sweep).
Outputs `output/facts/presence_lod_raw.csv` and `output/facts/presence_lod_summary.csv`,
plus a console summary. Gitignore-friendly (`output/` is already gitignored).

**Key trick that removes the Step 14 dependency for this pass:** the simulator injects a
single global symmetric error rate `e` via `sample_allele_counts`
(`src/allomix/simulate.py:244`), giving an effective per-direction floor of `e/3` at a
donor-homozygous marker (pure donor, `f_h=0`, ALT appears at `~e/3`). So if we simulate with
rate `e` and tell the detector `e_i = e/3` for every marker, the test is exactly calibrated
by construction. Sweeping `e` (e.g. 1e-2, 3e-3, 1e-3, 3e-4) is then a clean proxy for "what a
clean per-site error background would buy", and directly demonstrates the headline claim
(LoD is set by the error floor) without any per-site machinery.

**Grid:**
- relatedness: `unrelated`, `sibling` (the two clinically bracketing cases)
- n_markers: 76 (our panel), optionally 200
- depth: 1000, 2000, 5000
- error rate `e`: 1e-2, 3e-3, 1e-3, 3e-4
- host fraction `f_h`: `0.0` (negative controls) plus `{1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3}`
- replicates: >=200 at `f_h=0` (need the FP tail), >=60 at each positive (more at the bottom)

**Per replicate:**
1. `gts = simulate.generate_related_genotypes(n_markers, relatedness, rng, maf_range=(0.2,0.5))`
   (`src/allomix/simulate.py:508`) — gives host and donor genotypes.
2. Build the admixture sample by blending at **donor fraction = `1 - f_h`** (host re-occurrence
   means host is the minor contributor). Use `simulate.blend_from_genotype_dicts(...)`
   (`src/allomix/simulate.py:914`, avoids temp files) or `blend_vcfs` (`:689`) with
   `error_rate=e`, and the existing realistic-noise knobs (`realistic_biases`,
   `locus_dropout_rate=0.016`, `depth_cv=0.43`) to match the issue #8 sweep. Confirm the exact
   signature when implementing.
3. Select donor-homozygous markers where host carries the donor-absent allele (Vynck types 0,
   1, 10, 11; see "Marker set"). For each, take `y_i` = donor-absent-allele read count
   (`admix_ad_alt` when donor is hom-ref, `admix_ad_ref` when donor is hom-alt), `n_i` = depth,
   `h_i` = host dose of that allele (1 if host het, 2 if host homozygous-opposite).
4. Background `e_i = e/3` for every marker (self-consistent with the simulator).
5. Statistic (inline): pooled Poisson `Y=sum y_i`, `Lam=sum n_i*e_i`,
   `p_pois = scipy.stats.poisson.sf(Y-1, Lam)`; plus the LRT in "The statistic" giving
   `f_h_hat`, `D`, `p_lrt = 0.5*chi2.sf(D,1)`. Record both.
6. For the MLE comparison, also run `chimerism.estimate_single_donor_bb` on the same markers
   and store `1 - donor_fraction` and the donor_fraction itself.

**Analysis:**
- Negative controls: report the empirical false-positive rate at alpha=0.05 (fraction of
  `f_h=0` replicates with `p_lrt < 0.05`); it should be ~0.05. Also report the realised
  decision threshold (95th percentile of `Y/Lam` or of the statistic) as the LoB analogue.
- Positive controls: detection rate per `f_h` = fraction with `p_lrt < 0.05` (and a variant
  using the LoB threshold from negatives). Fit the 95% LoD with `fit_lod`/`compute_lob` from
  `paper/scripts/run_lod_validation.py` (import or copy). Report presence-LoD per
  (relatedness, depth, e).
- MLE comparison: at matched cells, compute the MLE-LoB-based LoD (estimate-exceeds-LoB rule,
  as issue #8 does) and tabulate presence-LoD vs MLE-LoD.

**Acceptance gate (decides whether to build Step 14 + `detect.py` + CLI):**
1. Calibrated: FP rate ~= 0.05 on negatives across cells.
2. Detection rate monotone in `f_h`, depth, and decreasing `e`.
3. Presence-LoD below the MLE-LoB LoD at low `e` / high depth; the gap shrinks as `e` rises
   toward 1e-2 (showing the error floor is the lever, i.e. Step 14 is worth doing).
If 1-2 hold and 3 shows a real gap at clean error rates, proceed with the full plan. If not,
write up why in this file and stop.

> Note (see "Update (2026-05-28)" above): gate #3 runs on binomial simulation, which
> charges the MLE no overdispersion, so the measured presence-vs-MLE gap is a **lower
> bound**. A small gap here does not by itself justify stopping; the gap widens under the
> real overdispersion the MLE pays and this detector does not. Do not "fix" this by turning
> on a global `rho` in the prototype — that miscalibrates the presence null (it overdisperses
> the donor-absent background too). The honest demonstration needs marker-type-aware
> overdispersion (Step 21 follow-up).

**Sanity checks before the full grid:** run a 10-replicate pilot first;
`f_h=0` should give few/no detections, `f_h=1e-3` at 2000x with `e=3e-4` should detect almost
always, and `f_h_hat` should track `1 - donor_fraction` from the MLE on the same sample.

## Fix scope

New module plus reporting and CLI plumbing. The detector is a separate path from the
estimator, so it does not touch the MLE likelihood functions. Build the control data and
clear the calibration gate above before the CLI/report integration below.

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

1. **Control data + calibration gate** (do first): build the positive and negative controls
   and clear the calibration gate described in "Control data and calibration (build this
   first)" above. This is the primary evidence for the LoD claim and must pass before the
   detector is wired into the CLI.
2. **Unit sanity** (cheapest): `pytest tests/test_detect.py -x -q`.
3. **Synthetic LoD comparison**: fold the presence-test detection into the issue #8 LoD sweep
   (`paper/scripts/run_lod_validation.py`, `claude/issue_8_lod_markers_plan.md`) so the
   presence-test LoD sits on the same axes as the existing MLE-LoB LoD. Run with a per-site
   error table provided (perfect rates, matching how the issue #8 sweep provides perfect
   biases) to show the achievable gain, and again with the global rate to show the floor
   effect. Expectation: at low fractions and high depth the presence test detects below the
   MLE-LoB LoD when per-site errors are available; the two converge under the global rate.
4. **Real-data smoke test**: run on the post-HSCT samples in
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

- [ ] **Control data + calibration first** — `paper/scripts/run_presence_lod.py` (NEW, or a
      prototype): generate extremely-low positive controls and error-only negative controls,
      run `host_presence_test`, report false-positive rate vs alpha and the fitted 95% LoD
      (reuse `compute_lob` / `fit_lod` from `run_lod_validation.py`). Clear the calibration
      gate before the integration items below.
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
  background, **not required**: the detector uses AD counts plus per-site error rates, and
  the controls declare their error rate rather than deriving it from quality scores.
- **issue #8 LoD sweep** — `claude/issue_8_lod_markers_plan.md`. Reused to build the
  positive/negative controls and validate the detection-LoD gain. The control-data and
  calibration step is sequenced first, before the detector is wired into the CLI.
- **Simulator per-site error extension** — co-requisite with the Step 14 simulator
  follow-up. Only needed to validate the per-site advantage; the first calibration pass runs
  on the existing scalar global error model.
