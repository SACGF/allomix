# Plan: Per-site dropout rate in the MLE likelihood

Status: proposed, not implemented.

## TL;DR

allomix currently treats every informative marker that happens to call in the
admixture sample as fully trustworthy. Some panel sites drop out (./.) far
more often than others, and a call at a flaky site is informative but less
so. The bias-training cohort already gives us per-site no-call rates (see
`scripts/measure_panel_bias.py:340-358`, which reports them). We want to
estimate a per-site dropout probability `d_s`, store it in a TSV alongside
the existing bias table, load it optionally into `monitor`, and downweight
each marker's likelihood contribution by `w_s = 1 - d_s`. Sites missing from
the table behave as `d_s = 0` (backwards-compatible).

This is a quality-of-fit improvement, not a point-estimate improvement on
its own. It reduces the influence of chronically unreliable sites on both
the MLE and the profile-likelihood CIs, and feeds naturally into the
planned beta-binomial GoF (Step 13).

## Background

### What "dropout" means here

`scripts/measure_panel_bias.py` distinguishes two signals:

1. **Locus dropout / no-call rate** — fraction of samples where the site is
   `./.` (`measure_panel_bias.py:164-168`). Empirically ~1.6% on the IDT
   rhAmpSeq panel, with a long tail (see `paper/scripts/generate_supp_synthetic.py`
   references).
2. **Allele dropout** — het calls under-represented vs HWE expectation
   (`measure_panel_bias.py:390-399`). This is a separate, already-measured
   but currently-unused signal.

Step 15 in `claude/allomix_overall_plan.md` is about (1): "per-site no-call
rates ... flaky sites are automatically downweighted rather than treated as
fully informative when they happen to call". Allele-dropout modelling is a
distinct follow-up and is out of scope here.

### Why it matters

At present the MLE loops over informative markers in
`total_log_likelihood_bb` (`src/allomix/chimerism.py:185-211`) and sums
per-marker log-likelihoods with equal weight. A site that calls on this
sample but has historically dropped out 25% of the time contributes as much
to the fit as a site that always calls cleanly. The 25%-dropout site is
more likely to produce noisy or aberrant calls even when it does call
(e.g. low-complexity region, secondary-structure amplification artefact,
off-target capture). We should downweight it.

The training inputs already exist:

- `scripts/measure_panel_bias.py:340-358` — reports per-marker no-call rate
  across the cohort, writes it into `<prefix>_per_marker.tsv` as the
  `nocall_rate` column.
- `scripts/qc_bias_samples.py` — filters samples by whole-sample no-call
  rate before they enter bias training.
- `src/allomix/bias.py` — the existing "estimate from a cohort, save a TSV,
  load at monitor time" pattern we will mirror.

## Modelling choice

Two candidate integrations; pick (A) for the first cut.

### (A) Weighted likelihood (preferred)

Scale each marker's log-likelihood by a reliability weight `w_s = 1 - d_s`:

```
total_ll = sum_s  (1 - d_s) * ll_marker_bb(k_s, n_s, p_s, rho)
```

This treats `d_s` as a confidence discount. Sites with zero historical
dropout get full weight; a site with `d_s = 0.3` contributes 70% of its
normal log-likelihood; `d_s >= 0.95` sites are effectively excluded.
Trivial to implement, composes with the existing bias and (future) per-site
error-rate plumbing, and has no identifiability issues.

The weighting is not a true probabilistic density (the weighted sum is not
a log-density), so strictly it is a pseudo-likelihood. That is fine for
point estimation; for CIs via profile likelihood we keep the same threshold
`chi2.ppf(0.95, df=1)/2`, accepting that coverage may drift slightly. We
already tolerate this kind of slack — bias correction is also pseudo —
and the paper notes it (`paper/discussion.md`).

### (B) Dropout mixture (alternative, heavier)

Treat each call at a flaky site as a mixture:

```
P(k | n, ...) = (1 - d_s) * BB(k | n, p(f), rho) + d_s * Noise(k | n)
```

where `Noise(k | n)` is a flat / uniform model for aberrant calls
(e.g. `Uniform[0, n]` or a wide beta with mean 0.5 and low concentration).
Proper likelihood, handles the "when it does call, sometimes the call is
junk" intuition directly.

Downsides: `d_s` now actively fights `rho` (overdispersion vs mixture
weight) during joint fitting; needs a noise model we don't really have
data to calibrate; numerically `log(a + b)` inside the sum means
`log_likelihood_marker_bb` can't just be scaled, it has to be rewritten
with `logsumexp`. Worth revisiting once the weighted-likelihood version is
in and we have real data showing it's insufficient. Document as follow-up,
do not implement now.

## Scope of the change

Four files touched, one new module, one new subcommand, one new optional
argument on `monitor`/`timeline`.

### 1. New module: `src/allomix/dropout.py`

Mirrors `src/allomix/bias.py`. Responsibilities:

- Dataclass for a per-marker dropout estimate.
- `estimate_dropouts(marker_lists, ...)` trainer.
- `save_dropout_table(...)` / `load_dropout_table(...)` I/O.

Proposed implementation:

```python
"""Per-marker locus-dropout rate estimation.

A locus-dropout event is a no-call ('./.')  for a sample at a marker that
is otherwise part of the panel. The per-marker dropout rate d_s is the
fraction of training samples that produce a no-call at site s. This file
is the training-side counterpart to src/allomix/bias.py and is consumed
by estimate_single_donor_bb / estimate_multi_donor via the
`marker_dropouts` argument.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

# Type alias kept local to avoid a circular import with bias.py
MarkerKey = tuple[str, int, str, str]


@dataclass
class MarkerDropout:
    """Per-marker locus-dropout estimate."""

    chrom: str
    pos: int
    ref: str
    alt: str
    dropout_rate: float  # fraction of samples with no-call at this site
    n_nocall: int        # count of no-calls in the training cohort
    n_total: int         # total samples that saw this marker
```

`estimate_dropouts` has a different input shape from `estimate_biases`:
bias estimation needs called-only hets, whereas dropout estimation needs
the no-call signal which `parse_vcf` currently *discards*
(`src/allomix/genotype.py:103-105`). Two workable approaches:

1. **Parse VCFs directly in `dropout.py`** using cyvcf2, mirroring the
   no-call detection loop in `scripts/measure_panel_bias.py:148-168`.
   Keeps `parse_vcf` simple but duplicates a small amount of VCF plumbing.

2. **Extend `parse_vcf` to optionally emit a no-call sentinel** (e.g. a
   `MarkerData` with `gt=(-1,-1)` and `dp=0`). Cleaner but risks
   non-trivially breaking every downstream caller that assumes GT is a
   real diploid call.

Option 1 is the lower-risk choice; keep `parse_vcf` as-is and write a
dedicated `_iter_nocalls(vcf_path, sample)` helper inside `dropout.py`:

```python
from cyvcf2 import VCF

def _iter_genotype_presence(
    vcf_path: Path | str,
    sample: str | int = 0,
) -> list[tuple[MarkerKey, bool]]:
    """Yield (marker_key, is_called) for every biallelic record, including no-calls."""
    vcf = VCF(str(vcf_path))
    if isinstance(sample, str):
        if sample not in vcf.samples:
            raise ValueError(f"Sample '{sample}' not in VCF. Available: {list(vcf.samples)}")
        sample_idx = list(vcf.samples).index(sample)
    else:
        sample_idx = sample

    out: list[tuple[MarkerKey, bool]] = []
    for variant in vcf:
        if len(variant.ALT) > 1:
            continue
        alt = variant.ALT[0] if variant.ALT else "."
        key: MarkerKey = (variant.CHROM, variant.POS, variant.REF, alt)
        gt = variant.genotypes[sample_idx]
        is_called = gt[0] >= 0 and gt[1] >= 0
        out.append((key, is_called))
    vcf.close()
    return out


def estimate_dropouts(
    vcf_paths_and_samples: list[tuple[Path | str, str | int]],
    min_samples: int = 5,
) -> dict[MarkerKey, MarkerDropout]:
    """Estimate locus-dropout rate per marker from a training cohort.

    Args:
        vcf_paths_and_samples: List of (vcf_path, sample) pairs. One pair
            per training sample. The sample may be a string name or an
            integer column index (matches parse_vcf semantics).
        min_samples: Minimum training samples a marker must appear in
            before a rate is reported. Markers below this threshold are
            omitted from the output (caller treats them as "unknown"
            downstream, i.e. d_s = 0).

    Returns:
        Dict mapping (chrom, pos, ref, alt) to MarkerDropout.
    """
    nocall_counts: dict[MarkerKey, int] = {}
    total_counts: dict[MarkerKey, int] = {}
    info: dict[MarkerKey, tuple[str, int, str, str]] = {}

    for vcf_path, sample in vcf_paths_and_samples:
        for key, is_called in _iter_genotype_presence(vcf_path, sample):
            total_counts[key] = total_counts.get(key, 0) + 1
            if not is_called:
                nocall_counts[key] = nocall_counts.get(key, 0) + 1
            info[key] = key

    out: dict[MarkerKey, MarkerDropout] = {}
    for key, n_total in total_counts.items():
        if n_total < min_samples:
            continue
        n_nocall = nocall_counts.get(key, 0)
        chrom, pos, ref, alt = info[key]
        out[key] = MarkerDropout(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            dropout_rate=n_nocall / n_total,
            n_nocall=n_nocall,
            n_total=n_total,
        )
    return out
```

TSV schema (stable, compatible with `save_bias_table` style):

```
chrom    pos    ref    alt    dropout_rate    n_nocall    n_total
chr1     12345  A      T      0.016           3           187
```

`save_dropout_table` / `load_dropout_table` are copy-paste from `bias.py`
with the field names changed. `load_dropout_table` returns
`dict[MarkerKey, float]` (the rate alone) so the MLE only sees the number
it needs.

Reference implementations to mirror exactly: `src/allomix/bias.py:95-127`.

### 2. `src/allomix/chimerism.py` — accept and apply dropout weights

Two choices about *where* in the likelihood to apply the weight. The
cheaper and clearer choice is inside the two `total_log_likelihood_*`
functions, symmetric with the existing `marker_biases` plumbing.

#### `log_likelihood_marker_bb` — no signature change needed

Leave this function alone. Weighting is a property of the *marker in the
panel*, not of the per-read likelihood. Applying it at the aggregator
keeps the per-read math clean and lets us unit-test it independently.

#### `total_log_likelihood_bb` (lines 185-211)

Current:

```python
def total_log_likelihood_bb(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    ll = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        ll += log_likelihood_marker_bb(m.admix_ad_ref, m.admix_ad_alt, w, error_rate, rho)
    return ll
```

Proposed:

```python
def total_log_likelihood_bb(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_dropouts: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """...
    marker_dropouts: Optional per-marker locus-dropout probability dict.
        When provided, each marker's log-likelihood contribution is
        multiplied by (1 - d_s). Missing keys are treated as d_s = 0
        (full weight).
    """
    ll = 0.0
    for m in markers:
        key = (m.chrom, m.pos, m.ref, m.alt)
        bias = marker_biases.get(key, 0.0) if marker_biases is not None else 0.0
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        ll_marker = log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w, error_rate, rho
        )
        if marker_dropouts is not None:
            d = marker_dropouts.get(key, 0.0)
            # Defensive clamp: training can spit out rates at 0 or 1.
            d = max(0.0, min(1.0, d))
            ll_marker *= (1.0 - d)
        ll += ll_marker
    return ll
```

Same change pattern in `total_log_likelihood_multi_bb` (lines 214-240),
adding `marker_dropouts` as a kwarg and applying the same `(1 - d)`
scaling per marker.

#### Estimators: `estimate_single_donor_bb` and `estimate_multi_donor`

Add `marker_dropouts` kwarg, thread through to every `total_log_likelihood_*`
call, and return a ChimerismResult/MultiDonorResult unchanged in shape.
No change to `expected_weight` / `expected_weight_multi`.

`estimate_single_donor_bb` (line 341) adds a kwarg and passes it to:

- The grid-search `total_log_likelihood_bb` call (inside `minimize_scalar` at
  line 386).
- `neg_ll_joint` (line 399).
- The profile-likelihood `minimize_scalar` at line 427.

`estimate_multi_donor` (line 476) adds a kwarg and passes it to:

- The grid-search `total_log_likelihood_multi_bb` call (line 541).
- `neg_ll` (line 556).
- `_profile_likelihood_cis_multi` (line 573), which also needs the kwarg
  plumbed through both inner optimisers (lines 635 and 654).

Keep the `marker_biases` and `marker_dropouts` kwargs positioned next to
each other; they're logically paired.

#### Per-marker residuals path

`_compute_per_marker_results` (line 278) and `_per_marker_results_multi`
(line 691) do NOT need the dropout weights. Residuals are a report-level
diagnostic and are already outlier-flagged via a 3-SD rule
(`chimerism.py:317-336`). Applying dropout weights there would double-count
the discount (the 3-SD test would shift because reported residuals are
unscaled). Leave them alone. Flag in a docstring that the per-marker
residuals do not reflect dropout weighting.

### 3. `src/allomix/cli.py`

Two additions.

#### (a) `--dropout-table` on `monitor` / `timeline`

Extend `_add_common_args` (line 17) near the existing `--bias-table`
block:

```python
parser.add_argument(
    "--dropout-table",
    default=None,
    help="Per-marker locus-dropout rate TSV (from allomix estimate-dropout). "
         "Sites missing from the table are treated as having zero dropout.",
)
parser.add_argument(
    "--no-dropout-correction",
    action="store_true",
    help="Disable dropout weighting even when a dropout table is provided",
)
```

Add a sibling helper:

```python
def _load_dropouts(args: argparse.Namespace) -> dict | None:
    if args.dropout_table and not args.no_dropout_correction:
        from allomix.dropout import load_dropout_table
        return load_dropout_table(args.dropout_table)
    return None
```

Thread into `_run_single_sample` as a new kwarg and pass to both
estimators:

```python
def _run_single_sample(
    host, donors, vcf_path, admix_sample,
    min_dp, min_gq, error_rate,
    marker_biases=None,
    marker_dropouts=None,
):
    ...
    if len(donors) == 1:
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
            marker_dropouts=marker_dropouts,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            marker_biases=marker_biases,
            marker_dropouts=marker_dropouts,
        )
```

Call sites `cmd_monitor` (line 128) and `cmd_timeline` (line 165) both
call `_load_dropouts(args)` once up front and pass the result into every
`_run_single_sample` call, symmetrical with the existing bias path.

#### (b) `estimate-dropout` subcommand

Mirror of the `estimate-bias` subcommand at lines 261-294. Accepts either
`--vcfs file1.vcf file2.vcf ...` (one sample per file) or
`--vcf joint.vcf --samples name1 name2 ...`:

```python
def cmd_estimate_dropout(args: argparse.Namespace) -> int:
    if args.vcfs and args.vcf:
        raise SystemExit("Use either --vcfs or --vcf/--samples, not both")
    if not args.vcfs and not args.vcf:
        raise SystemExit("One of --vcfs or --vcf is required")
    if args.vcf and not args.samples:
        raise SystemExit("--samples is required when using --vcf")

    pairs: list[tuple[str, str | int]] = []
    if args.vcfs:
        pairs = [(v, 0) for v in args.vcfs]
        n_source = f"{len(args.vcfs)} VCFs"
    else:
        _validate_sample_names(args.vcf, args.samples)
        pairs = [(args.vcf, s) for s in args.samples]
        n_source = f"{len(args.samples)} samples from {args.vcf}"

    from allomix.dropout import estimate_dropouts, save_dropout_table
    dropouts = estimate_dropouts(pairs, min_samples=args.min_samples)
    save_dropout_table(dropouts, args.output)
    print(
        f"Estimated dropout rate for {len(dropouts)} markers from {n_source} "
        f"-> {args.output}",
        file=sys.stderr,
    )
    return 0
```

Argparse block mirrors the `estimate-bias` block exactly, with
`--min-samples` (default 5) replacing `--min-het`, and default output
`dropout_table.tsv`. Register in `main()` alongside the existing
subcommands.

### 4. `src/allomix/qc.py` / reporting

No required changes. Two small niceties worth considering:

- In the QC report, include a line for `mean_dropout_used` when a dropout
  table was applied. Helps users sanity-check that the table actually
  matched most sites. Plumb the dict into `assess_quality` and aggregate.
  Borderline — consider for Step-15-follow-up rather than the core fix.
- In the TSV header, if the dropout table was supplied, emit a `# dropout_table=...`
  comment so the output is self-describing. Same pattern would work for
  `bias_table` but is not done today, so leave alone for consistency.

### 5. Documentation

- Add a paragraph to `paper/methods.md` near the bias-correction
  description noting the weighted-likelihood dropout model: "Each
  marker's log-likelihood contribution is scaled by `1 - d_s`, where `d_s`
  is the per-marker locus-dropout probability estimated from an
  independent training cohort."
- Update README CLI examples (after `estimate-bias`):

  ```
  # Estimate per-marker dropout rate from training VCFs
  allomix estimate-dropout --vcfs train/*.vcf.gz -o dropout_table.tsv

  # Monitor with dropout weighting (combines with bias correction)
  allomix monitor --host ... --donor ... --sample ... \
      --bias-table bias.tsv --dropout-table dropout_table.tsv
  ```
- `scripts/measure_panel_bias.py` already reports no-call rates but does
  not write a dropout table. Consider a small `--dropout-output` flag
  there so the ad-hoc cohort-scan script can emit the table in one pass,
  but this is optional — `estimate-dropout` is the primary path.

## Tests

New file `tests/test_dropout.py`, modelled on `tests/test_bias.py`.

### Unit tests on the trainer

1. **No no-calls in cohort → rate 0.0.** Construct two VCFs (use the
   `tests/test_data/synthetic_*.vcf` fixtures or write a small `_make_vcf`
   helper similar to the one in `tests/test_simulate.py`) in which every
   marker is called in every sample. `estimate_dropouts` returns one
   entry per marker, all with `dropout_rate == 0.0`, `n_nocall == 0`.

2. **One marker always missing → rate 1.0.** Single marker emitted as
   `./.:.:.:.` in every training sample. `dropout_rate == 1.0`,
   `n_nocall == n_total`.

3. **Mixed cohort.** 10 samples, marker A called in 10/10, marker B
   called in 8/10. `d_A == 0.0`, `d_B == 0.2`, `n_total == 10` in both.

4. **`min_samples` threshold.** A marker present in only 3/10 samples
   (missing from the VCFs entirely in the other 7) is omitted when
   `min_samples=5` and included when `min_samples=1`. (Note: "missing
   from the VCF" and "present but no-called" are different — only the
   latter increments `n_nocall`; the former simply isn't seen. This is
   the same semantics as the existing `bias.estimate_biases`.)

5. **TSV round-trip.** `save_dropout_table` → `load_dropout_table` ==
   input rates.

### Unit tests on the likelihood plumbing

1. **`marker_dropouts=None` unchanged.** `total_log_likelihood_bb(m, f)`
   with and without a `None` dict gives the same value to machine
   precision.

2. **All-zero dropout is a no-op.** Passing `{key: 0.0 for key in keys}`
   matches `marker_dropouts=None`.

3. **Dropout=1.0 fully excludes.** Passing `{key: 1.0}` for a single
   marker reduces `total_log_likelihood_bb` by exactly that marker's
   unweighted contribution.

4. **Linearity.** For a single marker at `d=0.5`, the total LL decreases
   by exactly 0.5 * `log_likelihood_marker_bb(...)` compared to `d=0`.

### Integration test

One end-to-end test using the synthetic multi-donor VCFs in
`tests/test_data/multidonor/`: estimate with and without a handcrafted
dropout table that marks 10% of the informative markers at
`dropout_rate=0.9`. Assert that (a) both estimators still run, (b) point
estimates are within 1% of the no-dropout version (small data; the
weighted version should shift only slightly), (c) CIs from the weighted
version are wider-or-equal. A CI-narrowing would suggest an implementation
bug.

### Regression

Run the full suite. The existing 261 tests must stay green.
`tests/test_chimerism.py` and `tests/test_multidonor.py` should be
unaffected because `marker_dropouts` defaults to `None`.

## Verification plan

Ordered cheapest to most expensive.

### 1. Unit tests

```bash
pytest tests/test_dropout.py tests/test_chimerism.py tests/test_multidonor.py -x -q
```

### 2. Synthetic CLI sanity

Use the existing multi-donor VCFs:

```bash
.venv/bin/allomix estimate-dropout \
    --vcfs tests/test_data/multidonor/*.vcf \
    --output /tmp/dropout_table.tsv \
    --min-samples 1

.venv/bin/allomix monitor \
    --vcf tests/test_data/multidonor/sample_h_d1_d2_f1_20_f2_10.vcf \
    --host-sample HOST --donor-sample DONOR1 --donor-sample DONOR2 \
    --sample ADMIX \
    --dropout-table /tmp/dropout_table.tsv \
    --format json | jq '.donor_fractions,.dropout_applied'
```

The synthetic cohort has no no-calls, so every `dropout_rate` in the
table is 0. The output must match a no-dropout run bit-for-bit.

### 3. Empirical cohort

The IDT rhAmpSeq cohort already reports per-marker no-call rates via
`scripts/measure_panel_bias.py`. Pipe those into the new trainer (or add
the `--dropout-output` side flag mentioned above):

```bash
find /tau/data/clinical_hg38/idt_rhampseq_sid/ \
    -path '*/2_variants/*.gatk.hg38.vcf.gz' \
    -not -path '*/gatk_per_sample/*' > /tmp/vcf_list.txt

# Example: call estimate-dropout on the full cohort
.venv/bin/allomix estimate-dropout \
    --vcfs $(cat /tmp/vcf_list.txt | head -n 200) \
    --output output/dropout_training/dropout_table.tsv \
    --min-samples 20
```

Compare the resulting table's `dropout_rate` column to the `nocall_rate`
column in the output of `scripts/measure_panel_bias.py` for the same
cohort. They should be equal up to the `min_samples` filter.

Then rerun the April-24 validation batch with and without the dropout
table:

```bash
python scripts/run_xls_batch.py output/Chimerism\ project\ patient\ list.xlsx \
    --vcf output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz \
    --host-column="NGS Sample ID" \
    --donor-column="NGS sample ID TP1" \
    --test-sample-column "NGS sample ID TP2" \
    --output-dir output/validation_run_dropout \
    --copy-columns="Donor,Chimerism result TP2" \
    --bias-table-tsv output/bias_training/bias_table.tsv \
    --dropout-table-tsv output/dropout_training/dropout_table.tsv
```

Compare `donor_fraction`, CI width, and `gof_pval` vs
`output/validation_run_bb_gof/batch.tsv` (the Step-13 baseline). Expected:

- Point estimates within ~0.5% of the no-dropout baseline.
- Slightly wider CIs on samples where the fit leaned on flaky sites.
- Slightly *higher* `gof_pval` on problematic samples once Step 13 is in.

(`scripts/run_xls_batch.py` would need matching `--dropout-table-tsv`
plumbing. Flag as a downstream follow-up if that script isn't touched as
part of this change.)

### 4. Paper benchmarks

Re-run `paper/Snakefile` once the dropout table exists. None of the
validation scripts (`paper/scripts/run_depth_validation.py`,
`run_relatedness_validation.py`, etc.) currently consume a dropout table;
they operate on fully simulated data where true dropout rate is known.
No paper regeneration needed unless we want a new figure showing the
effect of the dropout-weighted likelihood on real data — this belongs in
the real-sample validation (Step 11).

## Edge cases and risks

- **Missing table entries.** A marker present in the admixture sample but
  not in the training cohort (e.g. panel update) gets `d = 0.0` via the
  `.get(key, 0.0)` fallback. Behaviour identical to today. Safe default.
- **Rate = 1.0 but site called.** Can happen if the training cohort
  somehow never called the site but this admixture did (e.g. new capture
  chemistry). The weight becomes zero and the marker contributes nothing
  — equivalent to a hard exclusion. Acceptable.
- **Training cohort too small.** `min_samples` filters these out so
  they inherit the `d=0` default. Document clearly in the CLI help:
  "Sites missing from the table are treated as having zero dropout."
- **Composition with bias correction.** No interaction: bias shifts
  `expected_weight`, dropout scales the per-marker log-likelihood. They
  multiply cleanly.
- **Composition with Step-13 beta-binomial GoF.** The GoF test
  (`src/allomix/qc.py:74-138`) does *not* currently see per-marker
  weights. If we introduce dropout weighting, the GoF chi-sq is still
  computed over *all* included markers with equal weight. That is
  inconsistent with the MLE but is fine as a first cut: the GoF is a
  diagnostic, not a fit objective. If it becomes a concern, the natural
  change is `chi_sq += (1 - d_s) * residual^2 / var`, with
  `df = sum(1 - d_s) - n_fitted_params`, but this is follow-up. The
  beta-binomial plan explicitly leaves "per-marker weighted GoF" out of
  scope.
- **Composition with Step 16 (GQ-weighted contributions).** Both dropout
  weighting and GQ weighting scale the per-marker log-likelihood by a
  reliability factor. The obvious combination is `w_s = (1 - d_s) * gq_scale(m.gq)`.
  No conflict, and the single-scalar-per-marker plumbing here is what
  Step 16 will extend.
- **Profile-likelihood CIs.** The weighted log-likelihood is a
  pseudo-likelihood, so the chi2(df=1) threshold is approximate. For
  real panels with typical dropout <5% the distortion is negligible.
  Document in `paper/methods.md` when the change lands.
- **Simulate path.** `src/allomix/simulate.py:698-813` already models
  locus dropout; tests that use `locus_dropout_rate > 0` will produce
  VCFs where the affected markers are simply absent from the output (the
  simulator skips writing them, `simulate.py:813-815`). That's the
  correct signal for `estimate_dropouts` if we feed such simulated VCFs
  to the trainer — but the dropout information disappears if the
  synthetic sample has only one replicate per marker. To exercise the
  trainer in a simulated setting, generate N replicates per marker and
  aggregate. Doable inside the integration test.

## File-by-file checklist

- [ ] `src/allomix/dropout.py` — new module with `MarkerDropout`,
      `estimate_dropouts`, `save_dropout_table`, `load_dropout_table`.
- [ ] `src/allomix/chimerism.py` — add `marker_dropouts` kwarg to
      `total_log_likelihood_bb` (line 185), `total_log_likelihood_multi_bb`
      (line 214), `estimate_single_donor_bb` (line 341), and
      `estimate_multi_donor` (line 476); apply `(1 - d)` scaling at the
      aggregator; do not modify `log_likelihood_marker_bb` or the
      per-marker-residuals functions.
- [ ] `src/allomix/cli.py` — add `--dropout-table` and
      `--no-dropout-correction` to `_add_common_args` (line 17), add
      `_load_dropouts` helper, thread through `_run_single_sample`
      (line 73) and `cmd_monitor` / `cmd_timeline` (lines 128, 165), add
      `cmd_estimate_dropout` + argparse block (alongside lines 261-294),
      register in `main()` (around line 296).
- [ ] `tests/test_dropout.py` — unit tests on the trainer + likelihood
      plumbing; integration test against `tests/test_data/multidonor/`.
- [ ] README — add the two CLI examples from the Documentation section.
- [ ] `paper/methods.md` — add one sentence on weighted-likelihood
      dropout.
- [ ] Optional follow-up: expose `--dropout-table-tsv` in
      `scripts/run_xls_batch.py` for the real-cohort validation rerun.

## Out of scope / follow-ups

- Allele-dropout modelling (het-looks-like-hom via minor-allele loss).
  This is the het/HWE-ratio signal in `measure_panel_bias.py:390-399`
  and is a separate modelling change (mixture over genotype classes, not
  just a weight). Track as its own step.
- Dropout mixture (option B in Modelling choice). Revisit if the
  weighted-likelihood version fails to improve real-data GoF once Step
  13 is live.
- Surfacing `d_s` in the `--verbose` per-marker TSV output as an extra
  column. Useful once real-data results are in; low value pre-launch.
- Weighted GoF in `qc._compute_gof_pval` to stay consistent with the
  weighted MLE. See interaction note above; fold into Step 13 follow-up.
- `scripts/measure_panel_bias.py --dropout-output` side flag so the
  single-pass cohort scan can emit both the bias and dropout tables.
  Pure ergonomics.
