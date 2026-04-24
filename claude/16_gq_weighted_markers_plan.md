# Plan: GQ-weighted marker contributions

Status: proposed, not implemented.

## TL;DR

`--min-gq 20` is currently a hard pass/fail in two places (`parse_vcf`,
`classify_markers`). A GQ=19 host call is thrown away entirely; a GQ=21
call gets full weight. The Phred interpretation of GQ is continuous
(`P(call wrong) = 10^(-GQ/10)`), so the natural generalisation is to
keep borderline-confidence calls but scale each marker's log-likelihood
contribution by a reliability weight derived from the GQ.

The concrete change is small and composes cleanly with Step 15's dropout
weighting:

- `classify_markers` keeps markers above the active hard floor and
  attaches `host_gq` + `donor_gqs` to every `InformativeMarker`. The
  existing `--min-gq 20` default is unchanged; when `--gq-weighted` is
  on, the hard floor is set by a separate `--gq-floor` (default 10)
  instead, so GQ 10–20 calls are retained at reduced weight.
- A new `gq_weight(host_gq, donor_gqs, gq_floor, scheme, gq_full)`
  helper returns `w_gq in [0, 1]` per marker.
- The estimators (`estimate_single_donor_bb`, `estimate_multi_donor`)
  own user-facing kwargs `use_gq_weight`, `gq_floor`,
  `gq_weight_scheme`, precompute a `gq_weights: list[float]` once
  per call, and pass it to the aggregators. The aggregators
  (`total_log_likelihood_bb`, `total_log_likelihood_multi_bb`) take
  that precomputed list and multiply each marker's LL by the
  corresponding weight. Composes multiplicatively with Step 15's
  dropout weights (the combined per-marker factor is
  `(1 - d_s) * w_gq(m)`).
- CLI gains `--gq-weighted`, `--gq-weight-scheme`, and `--gq-floor`.
  `--gq-weighted` is a boolean opt-in; default off preserves current
  behaviour exactly. When on, the effective hard GQ floor is
  `args.gq_floor` (i.e. `--gq-weighted` replaces `--min-gq` for the
  purposes of the hard cut; `--min-gq` is ignored for that run).

Point-estimate gains are expected to be small (most real-data markers
are GQ >> 30). The motivation is: (1) keep more markers on low-quality
samples where we currently lose informative sites to the 20 cutoff, and
(2) treat GQ monotonically rather than as a cliff, which is the right
thing in principle and shows up on real-data edge cases (patients with
low host-genotyping coverage).

## Background

### What GQ means

GQ in GATK is defined as the Phred-scaled ratio of the two highest
genotype likelihoods: `GQ = PL[second-best] - PL[best]`. Under a flat
prior over candidate genotypes, the posterior probability that the
reported genotype is correct is

```
P(GT correct) ≈ 1 / (1 + 10^(-GQ/10))
```

This is the form used below. A common shortcut in the literature is
`1 - 10^(-GQ/10)`, which agrees to three decimals for GQ ≥ 20 but
diverges at low GQ (GQ=5: proper 0.760, shortcut 0.684; GQ=0: proper
0.5, shortcut 0). Use the proper posterior since we specifically want
GQ=10–20 markers to behave sensibly.

Values we see in practice on the idt_rhampseq_sid panel at >1000x:
- Typical: GQ 99 (capped in GATK), P(wrong) ≈ 10⁻¹⁰
- Borderline het call: GQ 20–40
- Low-depth / marginal call: GQ 5–20

For chimerism estimation the *genotype* enters the likelihood via
`expected_weight(host_gt, donor_gt, f_donor)` in
`src/allomix/chimerism.py:97`. If the host or donor genotype is wrong,
the expected weight is wrong, and the marker contributes misleading
evidence. So the confidence we should assign to the marker's LL
contribution tracks the joint probability that *all* contributing
genotypes (host + every donor) are correct:

```
P(all GTs correct) = P(host GT correct) * prod_i P(donor_i GT correct)
                   = (1/(1 + 10^(-host_gq/10))) * prod_i (1/(1 + 10^(-donor_i_gq/10)))
```

Under independence across contributors this is a product of per-sample
posteriors. Admixture-sample GQ is *not* used: for a mixture the
"called genotype" isn't a meaningful quantity, and the likelihood
already models mixture uncertainty through `f_donor` and `rho`.
Current code is consistent with this — `classify_markers` filters only
host/donor GQ, not admix GQ (`src/allomix/genotype.py:264-270`).

### Current hard-pass/fail behaviour

Two GQ cutoffs exist:

1. **`parse_vcf(..., min_gq=...)`**
   (`src/allomix/genotype.py:66-67,132`) — drops the call entirely
   before it enters any downstream code. Called from `cmd_monitor` /
   `cmd_timeline` with `args.min_gq` (`cli.py:136,137,173,174`).
2. **`classify_markers(..., min_gq=20)`**
   (`src/allomix/genotype.py:211,264-270`) — after joining host +
   donors + admix by position, re-checks host/donor GQ and drops the
   whole marker if any fail.

(1) is redundant when (2) is also applied with the same threshold
(both points only drop, never add information). (2) is the meaningful
check. With GQ-weighting enabled we want (2) to keep markers above a
low floor (say GQ 10) rather than the current 20.

### Why expect small gains

The clinical panel (76 SNPs, >1000x) mostly yields GQ=99 calls. The
difference between GQ 99 and GQ 30 is a posterior of ~1.0 vs 0.999 —
numerically indistinguishable for the MLE. The lift comes from two
places:

- **Recovered markers** — samples with weaker host genotyping (lower
  coverage or contamination) currently lose a handful of informative
  sites to the GQ 20 cutoff. With the floor dropped to 10 and weights
  applied, those sites come back with 0.9–0.99 weight, adding a small
  amount of information.
- **Dissuasion from marginal calls** — on edge-case samples where a
  host GT at a key informative site is borderline, the weighted
  likelihood self-protects by contributing less evidence from that
  marker rather than all-or-nothing.

This is Step 16 and the overall plan notes "likely small gains
relative to Steps 14 and 15, but cheap to add once the per-marker
likelihood is already being modified". The architecture here is
designed so Step 15 (dropout) and Step 16 (GQ) can land in either
order and compose multiplicatively.

## Modelling choice

Pick (A) with the Phred weight function. (B) and (C) are recorded for
completeness.

### (A) Weighted pseudo-likelihood with Phred-derived weight (preferred)

Scale each marker's log-likelihood by the joint probability that all
contributing genotypes are correct:

```
w_gq(m) = (1 / (1 + 10^(-host_gq/10))) * prod_i (1 / (1 + 10^(-donor_i_gq/10)))
total_ll = sum_m  w_gq(m) * ll_marker_bb(m, f, rho)
```

With a lower floor to guard against numerically-bad GQs (0, NA, or
unrealistically low):

```
def gq_weight(host_gq, donor_gqs, gq_floor=10, gq_full=30):
    def phred_posterior(g):
        # GQ -> P(GT correct) under flat-prior two-genotype model
        return 1.0 / (1.0 + 10.0 ** (-g / 10.0))

    # None GQ => treat as fully-confident (legacy behaviour)
    if host_gq is None:
        h_w = 1.0
    elif host_gq < gq_floor:
        return 0.0
    else:
        h_w = phred_posterior(host_gq)

    w = h_w
    for dgq in donor_gqs:
        if dgq is None:
            continue
        if dgq < gq_floor:
            return 0.0
        w *= phred_posterior(dgq)
    return w
```

The `gq_full` argument is unused in the Phred form above (kept in the
signature for consistency with the linear alternative). Retain it so
flipping to the linear form is a one-line change.

**Why the proper posterior and not the `1 - 10^(-GQ/10)` shortcut?**
The two agree to three decimals at GQ ≥ 20 and diverge at low GQ
(GQ=5: 0.760 vs 0.684; GQ=0: 0.5 vs 0). Since the whole point of
`--gq-weighted` is to retain GQ 10–20 markers at reduced weight, and
the floor defaults to 10, the low-GQ behaviour matters. Use the
proper posterior.

**Why Phred and not a linear ramp?** GQ is an ordinal but its
Phred semantics give a principled mapping to a probability. A linear
ramp (`w = clamp((gq - 10) / 20, 0, 1)`) is equally valid as a
first-cut heuristic but less defensible in print. The probability
interpretation also lets us describe the composite weight with dropout
as "the joint probability the marker is informative and all genotypes
are correct", which reads well in the paper.

### (B) Linear ramp between `gq_floor` and `gq_full` (alternative)

```
w_gq(m) = clamp((m.min_gq - gq_floor) / (gq_full - gq_floor), 0, 1)
```

Uses the minimum GQ across host+donors for a single-sample rule.
Simpler, doesn't assume independence across contributors. Useful as a
fallback if the Phred product is empirically too aggressive (i.e.
`w_gq` concentrates near 1.0 for almost every marker and the method
degenerates to the status quo). Accept as a hidden `--gq-weight-scheme
linear` option; default is Phred.

### (C) Dropout-style mixture at the genotype level (out of scope)

Treat each marker as a mixture over the possible miscalled genotypes:

```
P(k|n,...) = P(GT correct) * BB(k|n, p(f, GT_correct), rho)
           + P(GT wrong)   * sum_{GT'} P(GT'|miscall) * BB(k|n, p(f, GT'), rho)
```

Proper likelihood. Needs a miscall-prior over neighbouring genotypes
(e.g. a het miscalled as hom-ref) which we don't have calibrated data
for. Revisit only if (A) underperforms on real data. Do not implement
now.

## Scope of the change

Four files touched, one new helper module, one new CLI flag group. No
new subcommand (unlike Steps 14/15, GQ-weighting doesn't need a
training step — it's computed at runtime from the VCFs already being
parsed).

### 1. `src/allomix/gq_weight.py` (NEW)

Small module so the weight function has one authoritative home and is
trivially unit-testable. Could live inside `chimerism.py`, but keeping
it separate matches the `bias.py` / (planned) `dropout.py` pattern and
makes future schemes (B/C) easy to drop in.

```python
"""Per-marker genotype-quality weights for the chimerism likelihood.

GQ is the Phred-scaled probability that a sample's reported genotype
is correct. At every informative marker, the MLE uses the host and
donor genotypes to compute an expected reference-allele weight. When
a genotype is borderline, that expected weight is itself uncertain,
and the marker's log-likelihood contribution should be discounted
accordingly.

The canonical scheme is the joint probability that every contributing
genotype is correct, using the flat-prior posterior
``P(GT correct) = 1 / (1 + 10^(-GQ/10))``:

    w_gq(m) = (1 / (1 + 10^(-host_gq / 10)))
            * prod_i (1 / (1 + 10^(-donor_i_gq / 10)))

with a hard-zero floor below ``gq_floor`` to guard against garbage
GQs. The admixture sample's GQ is intentionally ignored: for a mixture
the called genotype isn't meaningful.

Consumed by ``chimerism.total_log_likelihood_bb`` /
``chimerism.total_log_likelihood_multi_bb`` via a ``use_gq_weight``
kwarg. Missing GQ (None) is treated as fully confident, matching the
pre-Step-16 legacy behaviour.
"""

from __future__ import annotations

import math


def gq_weight(
    host_gq: int | None,
    donor_gqs: list[int | None],
    gq_floor: int = 10,
    scheme: str = "phred",
    gq_full: int = 30,
) -> float:
    """Per-marker reliability weight derived from host and donor GQs.

    Args:
        host_gq: Host-sample GQ at this marker. ``None`` means the VCF
            did not emit a GQ FORMAT field — treat as fully confident
            (weight contribution = 1.0) to preserve legacy behaviour.
        donor_gqs: Per-donor GQ at this marker, in the same order as
            ``InformativeMarker.donor_gts``. ``None`` entries are
            skipped (treated as 1.0).
        gq_floor: GQs strictly below this return a zero weight —
            effectively a hard exclusion for the worst calls. This is
            the replacement for the old ``min_gq`` hard cutoff.
        scheme: ``"phred"`` (default) uses Phred-interpreted posteriors;
            ``"linear"`` uses a (gq - gq_floor) / (gq_full - gq_floor)
            ramp, taking min across host+donors. Retained for ablation
            experiments and as a fallback.
        gq_full: Upper GQ for the ``"linear"`` scheme (clamped to 1.0
            at and above this value). Ignored by ``"phred"``.

    Returns:
        Weight in [0.0, 1.0].
    """
    if scheme == "phred":
        def phred_posterior(g: int) -> float:
            return 1.0 / (1.0 + 10.0 ** (-g / 10.0))

        if host_gq is not None and host_gq < gq_floor:
            return 0.0
        h_w = 1.0 if host_gq is None else phred_posterior(host_gq)
        w = h_w
        for dgq in donor_gqs:
            if dgq is None:
                continue
            if dgq < gq_floor:
                return 0.0
            w *= phred_posterior(dgq)
        return w

    if scheme == "linear":
        gqs = [host_gq] + list(donor_gqs)
        valid = [g for g in gqs if g is not None]
        if not valid:
            return 1.0
        min_gq = min(valid)
        if min_gq < gq_floor:
            return 0.0
        if min_gq >= gq_full:
            return 1.0
        return (min_gq - gq_floor) / (gq_full - gq_floor)

    raise ValueError(f"Unknown GQ weight scheme: {scheme!r}")
```

### 2. `src/allomix/genotype.py`

Two changes.

#### 2a. Carry host+donor GQ on `InformativeMarker`

Current (lines 32-47):

```python
@dataclass
class InformativeMarker:
    ...
    admix_dp: int
    marker_types: list[int | None] | None = None
    informative_for: list[bool] | None = None
```

Append two optional fields (defaults None for backward compat with
every existing test fixture):

```python
@dataclass
class InformativeMarker:
    """A marker where host and at least one donor have different genotypes."""
    chrom: str
    pos: int
    ref: str
    alt: str
    host_gt: tuple[int, int]
    donor_gts: list[tuple[int, int]]
    marker_type: int
    admix_ad_ref: int
    admix_ad_alt: int
    admix_dp: int
    marker_types: list[int | None] | None = None
    informative_for: list[bool] | None = None
    host_gq: int | None = None
    donor_gqs: list[int | None] | None = None
```

#### 2b. `classify_markers` — keep borderline markers, attach GQs

Current (lines 263-270):

```python
# Filter: host/donor GQ
if min_gq > 0:
    if h.gq is not None and h.gq < min_gq:
        n_filtered += 1
        continue
    if any(d.gq is not None and d.gq < min_gq for d in ds):
        n_filtered += 1
        continue
```

Proposed — split into a **hard floor** (always applied, exclude
garbage calls) and propagation of raw GQ through to the marker. The
soft weighting is applied later in the likelihood.

```python
# Filter: host/donor GQ (hard floor only — soft weighting happens
# in total_log_likelihood_bb when --gq-weighted is on).
if min_gq > 0:
    if h.gq is not None and h.gq < min_gq:
        n_filtered += 1
        continue
    if any(d.gq is not None and d.gq < min_gq for d in ds):
        n_filtered += 1
        continue
```

(No code change here — semantics are unchanged. The behavioural shift
comes from lowering the CLI default, documented in §3.)

Then, where the `InformativeMarker` is constructed (lines 287-302),
add the two new fields:

```python
informative.append(
    InformativeMarker(
        chrom=key[0],
        pos=key[1],
        ref=key[2],
        alt=key[3],
        host_gt=h.gt,
        donor_gts=donor_gts,
        marker_type=mtype_first,
        admix_ad_ref=a.ad_ref,
        admix_ad_alt=a.ad_alt,
        admix_dp=a.dp,
        marker_types=mtypes,
        informative_for=[mt is not None for mt in mtypes],
        host_gq=h.gq,
        donor_gqs=[d.gq for d in ds],
    )
)
```

### 3. `src/allomix/cli.py`

#### 3a. Add `--gq-weighted`, `--gq-weight-scheme`, `--gq-floor` to common args

Extend `_add_common_args` (line 17), placed next to the existing
bias-correction flags for symmetry. A single boolean `--gq-weighted`
is sufficient (no `--no-gq-weighting`) because weighting is off by
default and requires no auxiliary state to disable:

```python
parser.add_argument(
    "--gq-weighted",
    action="store_true",
    help="Weight per-marker log-likelihood by the joint probability "
         "that host and donor genotypes are correct, derived from GQ "
         "(Phred). When enabled, --min-gq is ignored and --gq-floor "
         "takes its place as the hard GQ floor; markers above the "
         "floor contribute at the per-contributor posterior "
         "1 / (1 + 10^(-gq/10)).",
)
parser.add_argument(
    "--gq-weight-scheme",
    choices=["phred", "linear"],
    default="phred",
    help="GQ weight function when --gq-weighted is on (default: phred).",
)
parser.add_argument(
    "--gq-floor",
    type=int,
    default=10,
    help="Hard GQ floor when --gq-weighted is on (default: 10). "
         "Replaces --min-gq for that run. Ignored when --gq-weighted "
         "is off.",
)
```

Design decision: when `--gq-weighted` is on, `--gq-floor` fully
replaces `--min-gq` rather than taking the min/max of the two. This
keeps the opt-in semantics clean — users who want soft weighting down
to GQ 10 get exactly that by passing `--gq-weighted`, regardless of
what `--min-gq` defaults to or is set to. If a user explicitly passes
both `--min-gq 30 --gq-weighted`, the `--min-gq 30` is silently
ignored; emit a warning in that case (see §3b).

#### 3b. Apply effective floor before parsing

Both `cmd_monitor` (line 128) and `cmd_timeline` (line 165) resolve
the effective hard GQ floor once up front:

```python
if args.gq_weighted:
    if args.min_gq != 20:  # 20 is the default; non-default implies explicit
        print(
            f"Note: --min-gq={args.min_gq} is ignored because --gq-weighted "
            f"is on; using --gq-floor={args.gq_floor} as the hard floor.",
            file=sys.stderr,
        )
    effective_min_gq = args.gq_floor
else:
    effective_min_gq = args.min_gq
```

and use `effective_min_gq` in place of `args.min_gq` for every
`parse_vcf(...)` and `classify_markers(...)` call below.

#### 3c. Thread `use_gq_weight` + kwargs into `_run_single_sample`

`_run_single_sample` (lines 73-111) grows new kwargs and passes them
to whichever estimator it picks. Keep the argument names short:

```python
def _run_single_sample(
    host: list,
    donors: list[list],
    vcf_path: str,
    admix_sample: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    use_gq_weight: bool = False,
    gq_floor: int = 10,
    gq_weight_scheme: str = "phred",
) -> tuple:
    """..."""
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0)
    genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
    genotypes.sample_name = admix_sample

    if len(donors) == 1:
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
            use_gq_weight=use_gq_weight,
            gq_floor=gq_floor,
            gq_weight_scheme=gq_weight_scheme,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            marker_biases=marker_biases,
            use_gq_weight=use_gq_weight,
            gq_floor=gq_floor,
            gq_weight_scheme=gq_weight_scheme,
        )
    qc = assess_quality(result, genotypes)
    return result, qc, genotypes
```

`cmd_monitor` and `cmd_timeline` pass the new kwargs through:

```python
result, qc, genotypes = _run_single_sample(
    host, donors, args.vcf, sample_name,
    args.min_dp, effective_min_gq, args.error_rate,
    marker_biases=marker_biases,
    use_gq_weight=args.gq_weighted,
    gq_floor=args.gq_floor,
    gq_weight_scheme=args.gq_weight_scheme,
)
```

Nothing else in CLI needs changing. No new subcommand. `estimate-bias`
doesn't need GQ-weighting (it operates on per-sample hets at
fixed-confidence thresholds) and is unaffected.

### 4. `src/allomix/chimerism.py`

All changes below parallel the Step 15 `marker_dropouts` plan — same
shape, same plumbing points, same fallback-to-legacy-when-None rule.
If Step 15 lands first, apply these as edits on top of its diff; if
Step 16 lands first, Step 15 will extend this same kwarg list.

#### 4a. `total_log_likelihood_bb` — apply per-marker GQ weight

Current (lines 185-211):

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
    gq_weights: list[float] | None = None,
) -> float:
    """Sum of per-marker beta-binomial log-likelihoods.

    When ``gq_weights`` is provided, each marker's log-likelihood
    contribution is multiplied by the corresponding precomputed weight
    in ``[0, 1]``. The list must be aligned 1-to-1 with ``markers``.
    When None (default), no weighting is applied (preserves pre-Step-16
    behaviour).

    ``gq_weights`` is computed once per estimator call from
    ``gq_weight(...)``; it does not depend on ``f_donor``, ``rho``, or
    ``error_rate``, so precomputing avoids tens of thousands of
    redundant calls inside the grid/Nelder-Mead/profile loops.
    """
    if gq_weights is not None and len(gq_weights) != len(markers):
        raise ValueError(
            f"gq_weights length {len(gq_weights)} != markers length {len(markers)}"
        )
    ll = 0.0
    for i, m in enumerate(markers):
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        ll_marker = log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w, error_rate, rho
        )
        if gq_weights is not None:
            ll_marker *= gq_weights[i]
        ll += ll_marker
    return ll
```

The estimator (`estimate_single_donor_bb`) owns the three user-facing
kwargs (`use_gq_weight`, `gq_floor`, `gq_weight_scheme`), computes
`gq_weights` once from `allomix.gq_weight.gq_weight`, and threads the
list into every `total_log_likelihood_bb` call site.

Add the import at the top of `chimerism.py`:

```python
from allomix.gq_weight import gq_weight
```

#### 4b. `total_log_likelihood_multi_bb` — identical shape

Current (lines 214-240). Same edit: add `gq_weights: list[float] |
None = None` kwarg and apply it as in 4a. `m.donor_gqs` has multiple
entries for multi-donor, handled at the estimator level when building
the precomputed `gq_weights` list.

**Performance note.** `total_log_likelihood_bb` runs inside a 1001-point
grid × an inner `minimize_scalar` over rho × a Nelder-Mead refinement
× a profile-likelihood scan — tens of thousands of calls per sample.
The per-marker `gq_weight` is a function of the marker alone and does
not depend on `f`, `rho`, or `error_rate`, so precompute it once per
estimator call:

```python
# Inside estimate_single_donor_bb and estimate_multi_donor, before the
# grid search starts:
if use_gq_weight:
    gq_weights = [
        gq_weight(
            m.host_gq,
            m.donor_gqs if m.donor_gqs is not None else [],
            gq_floor=gq_floor,
            scheme=gq_weight_scheme,
        )
        for m in markers
    ]
else:
    gq_weights = None
```

Pass `gq_weights` (a plain list aligned with `markers`) into
`total_log_likelihood_bb` / `total_log_likelihood_multi_bb` instead of
`use_gq_weight + gq_floor + gq_weight_scheme`. Inside the aggregator:

```python
for i, m in enumerate(markers):
    ...
    ll_marker = log_likelihood_marker_bb(...)
    if gq_weights is not None:
        ll_marker *= gq_weights[i]
    ll += ll_marker
```

Caller (the estimator) keeps the three `use_gq_weight / gq_floor /
gq_weight_scheme` kwargs, computes the list once, and threads the list
through. Keeps the estimator signature ergonomic, keeps the aggregator
hot-loop fast, and keeps tests on the aggregator unchanged (they can
pass `gq_weights=[...]` directly if they want to exercise the path).

#### 4c. `estimate_single_donor_bb` — add kwargs and thread through

Current signature (line 341):

```python
def estimate_single_donor_bb(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> ChimerismResult:
```

Proposed:

```python
def estimate_single_donor_bb(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    use_gq_weight: bool = False,
    gq_floor: int = 10,
    gq_weight_scheme: str = "phred",
) -> ChimerismResult:
```

Precompute `gq_weights` once at the top of the function (see the
Performance note above), then thread it into every
`total_log_likelihood_bb` call site inside this estimator:

- The grid-search profiling `minimize_scalar` body (line 384-390).
- The Nelder-Mead `neg_ll_joint` inner function (line 399-406).
- The profile-likelihood `profile_ll_f` inner function (line 424-433).

Each call goes from:

```python
total_log_likelihood_bb(markers, f, error_rate, math.exp(log_r), marker_biases)
```

to:

```python
total_log_likelihood_bb(
    markers, f, error_rate, math.exp(log_r), marker_biases,
    gq_weights=gq_weights,
)
```

Keep `_compute_per_marker_results` (line 278, called at line 455)
unchanged *for the residual values themselves* — residuals are an
unweighted diagnostic. But surface `gq_weight` as a new field on
`MarkerResult` for the verbose TSV (see §5). Mark in a one-liner
comment that the 3-SD outlier flag is computed against unweighted
residuals; if this becomes a concern later (e.g. low-GQ markers
flagged as outliers that the MLE already discounted), add a
gq-weight-aware variant as a follow-up.

#### 4d. `estimate_multi_donor` — same treatment

Add the three kwargs (`use_gq_weight`, `gq_floor`, `gq_weight_scheme`)
to `estimate_multi_donor` (line 476), precompute `gq_weights` once at
the top, and thread the list into:

- `total_log_likelihood_multi_bb` in the grid search (line 541-543).
- `total_log_likelihood_multi_bb` inside `neg_ll` (line 549-556).
- `_profile_likelihood_cis_multi` (line 573), which gains a
  `gq_weights: list[float] | None = None` parameter and forwards it
  into the two inner optimisers (lines 634-642 and 646-656).

Leave `_per_marker_results_multi` (line 576) with the same "residuals
are unweighted diagnostic" note as the single-donor helper.

### 5. `src/allomix/qc.py` and reporting

Three changes in v1. The GoF inconsistency would otherwise be a
regression on real data where GQs are genuinely heterogeneous.

#### 5a. `MarkerResult.gq_weight` field and verbose TSV

Add an optional field on `MarkerResult` (default `1.0` for backward
compat):

```python
@dataclass
class MarkerResult:
    ...
    included: bool
    gq_weight: float = 1.0  # per-marker GQ-derived reliability (1.0 = no weighting)
```

Populated by `_compute_per_marker_results` / `_per_marker_results_multi`
from the same `gq_weights` list passed into the aggregator. In
`report.py`'s verbose TSV, emit it as a new column next to the
existing `residual` / `included` columns. This gives users a direct
view of which sites were downweighted and by how much — follow-up
item #2 in the original v1 plan, promoted to v1 because the GoF
change below produces asymmetric per-marker behaviour that the user
should be able to audit.

#### 5b. `n_markers_effective` in `QCReport`

Add `n_markers_effective: float` alongside the existing `n_informative`
and `n_markers_used` in the QC summary, computed as
`sum(mr.gq_weight for mr in result.per_marker if mr.included)`. Surface
in the TSV summary header. No warnings gated on it for v1 — it's
informational.

#### 5c. Weighted GoF chi-squared

`_compute_gof_pval` (`qc.py:~74-138`) currently sums squared
standardised residuals with equal weight per marker. With the MLE
downweighting some markers, the unweighted chi-squared becomes
inconsistent — a low-GQ marker whose genotype is probably wrong will
look like a large residual and inflate chi-squared, rejecting the fit
when the MLE itself has already discounted that marker.

Pass `gq_weights` into `_compute_gof_pval` and compute

```python
chi_sq = sum(w_s * (resid_s**2 / var_s) for s, w_s in zip(...))
df = max(sum(gq_weights) - n_fitted_params, 1.0)
```

Use `sum(gq_weights)` as the *effective degrees of freedom* rather
than the integer count of included markers. `chi2.sf(chi_sq, df)`
with a float `df` is supported by `scipy.stats.chi2`. Fall back to
the current integer-df path when `gq_weights` is None (pre-Step-16
behaviour, identical numerically). Guard `df >= 1.0` to avoid the
degenerate case where every weight is near zero.

This is the same logic the Step 15 plan defers; we promote it here
because GQ weighting *will* produce non-trivial weight variance on
real data (low-coverage host samples yield meaningful GQ 10–30
markers), while locus dropout on the current panel is ~1.6% and the
weighted/unweighted GoF will differ less. The Step 15 implementation
should pick up this same extension when it lands — add a note in
`qc.py` that the weighted GoF signature accepts both `gq_weights` and
(future) `dropout_weights` composed multiplicatively.

### 6. Documentation

- README: add an example under the bias/dropout examples:

  ```bash
  allomix monitor --vcf joint.vcf \
      --host-sample HOST --donor-sample DONOR --sample ADMIX_D30 \
      --bias-table bias.tsv --gq-weighted
  ```

- `paper/methods.md`: one paragraph near the bias-correction /
  dropout-weighting description:

  > Each marker's log-likelihood contribution is additionally scaled by
  > a reliability weight derived from the Phred-scaled genotype quality
  > (GQ) of the host and donor genotypes, `w_gq(m) = prod_c
  > (1 − 10^(−GQ_c/10))` across the host and all donor contributors.
  > Missing GQ values default to 1.0. A hard floor (GQ 10) excludes
  > calls whose genotype is effectively random. Under independence
  > this is the joint probability that every contributing genotype is
  > correct, so the weighted sum is the natural generalisation of the
  > current hard `--min-gq 20` cutoff.

## Tests

New file `tests/test_gq_weight.py` and surgical additions to
`tests/test_chimerism.py`, `tests/test_multidonor.py`, and
`tests/test_cli.py`.

### `tests/test_gq_weight.py` — unit tests on the helper

```python
import math
import pytest
from allomix.gq_weight import gq_weight


def _pp(g):
    """Reference posterior used by the phred scheme."""
    return 1.0 / (1.0 + 10 ** (-g / 10))


class TestGqWeightPhred:
    def test_all_none_is_one(self):
        assert gq_weight(None, [None, None]) == 1.0

    def test_gq_99_is_essentially_one(self):
        w = gq_weight(99, [99])
        assert 1.0 - w < 1e-9

    def test_gq_20_matches_posterior(self):
        # Single contributor at GQ 20 => 1/(1 + 0.01) = 0.990099...
        w = gq_weight(20, [])
        assert math.isclose(w, _pp(20), rel_tol=1e-12)

    def test_multiple_contributors_multiply(self):
        w = gq_weight(20, [30])
        assert math.isclose(w, _pp(20) * _pp(30), rel_tol=1e-12)

    def test_below_floor_is_zero(self):
        assert gq_weight(5, [30], gq_floor=10) == 0.0
        assert gq_weight(30, [5], gq_floor=10) == 0.0

    def test_none_contributor_skipped(self):
        # Donor GQ None should contribute factor 1.0
        assert math.isclose(gq_weight(30, [None]), _pp(30), rel_tol=1e-12)


class TestGqWeightLinear:
    def test_linear_scheme_ramp(self):
        # floor=10, full=30: gq=20 => 0.5
        w = gq_weight(20, [], scheme="linear", gq_floor=10, gq_full=30)
        assert math.isclose(w, 0.5, rel_tol=1e-12)

    def test_linear_clamps_above_full(self):
        assert gq_weight(99, [], scheme="linear", gq_floor=10, gq_full=30) == 1.0

    def test_linear_uses_min_across_contributors(self):
        # min GQ = 20 (host); donor 40 ignored
        w = gq_weight(20, [40], scheme="linear", gq_floor=10, gq_full=30)
        assert math.isclose(w, 0.5, rel_tol=1e-12)

    def test_linear_below_floor_zero(self):
        assert gq_weight(5, [40], scheme="linear", gq_floor=10, gq_full=30) == 0.0


class TestInvalidScheme:
    def test_unknown_scheme_raises(self):
        with pytest.raises(ValueError):
            gq_weight(30, [], scheme="bogus")
```

### `tests/test_chimerism.py` — likelihood-plumbing regression + new

The aggregator takes a precomputed `gq_weights` list, so these tests
target it directly rather than round-tripping through `use_gq_weight`.
The `use_gq_weight` path is exercised end-to-end via
`estimate_single_donor_bb` in the integration test below.

1. **`gq_weights=None` unchanged (regression).** Existing tests pass
   with the new kwarg defaulted off.
2. **Weights of all 1.0 match None, exactly.** For any marker list `M`,
   `total_log_likelihood_bb(M, ..., gq_weights=None)` equals
   `total_log_likelihood_bb(M, ..., gq_weights=[1.0]*len(M))`
   bit-for-bit — the code path is identical except for the multiplication,
   and `x * 1.0 == x` in IEEE-754.
3. **Zero-weight excludes.** A two-marker call with
   `gq_weights=[1.0, 0.0]` equals the single-marker call on just the
   first marker, bit-for-bit.
4. **Linearity.** With `gq_weights=[1.0, 0.5]`, the total equals
   `ll(m1) + 0.5 * ll(m2)` to `rel_tol=1e-12`, where each `ll(m_i)`
   is computed separately.
5. **Length mismatch raises.** `len(gq_weights) != len(markers)`
   raises `ValueError`.
6. **End-to-end estimator regression.** `estimate_single_donor_bb`
   with all GQ=99 host + donor and `use_gq_weight=True` gives a point
   estimate within `rel_tol=1e-6` of `use_gq_weight=False`. (Tolerance
   reflects `~2 * 10^-10` per-marker deviation accumulated through a
   nonlinear optimizer over many markers, not a machine-precision
   claim. This is the test that my earlier 1e-9 draft got wrong.)

Example skeleton (extending `_make_marker` in `test_chimerism.py:30`):

```python
def test_gq_weights_all_ones_matches_none():
    ms = _make_markers_for_fraction(f_donor=0.3, n_markers=10, dp=500)
    ll_none = total_log_likelihood_bb(ms, f_donor=0.3)
    ll_ones = total_log_likelihood_bb(ms, f_donor=0.3, gq_weights=[1.0] * len(ms))
    assert ll_none == ll_ones  # exact


def test_gq_weights_zero_excludes_marker():
    m_keep = _make_marker(host_gt=(0, 0), donor_gt=(1, 1), ad_ref=500, ad_alt=500, pos=100)
    m_drop = _make_marker(host_gt=(0, 0), donor_gt=(1, 1), ad_ref=500, ad_alt=500, pos=200)
    ll_weighted = total_log_likelihood_bb(
        [m_keep, m_drop], f_donor=0.5, gq_weights=[1.0, 0.0],
    )
    ll_keep_only = total_log_likelihood_bb([m_keep], f_donor=0.5)
    assert ll_weighted == ll_keep_only  # exact


def test_gq_weights_linearity():
    m1 = _make_marker(host_gt=(0, 0), donor_gt=(1, 1), ad_ref=500, ad_alt=500, pos=100)
    m2 = _make_marker(host_gt=(0, 0), donor_gt=(1, 1), ad_ref=400, ad_alt=600, pos=200)
    ll1 = total_log_likelihood_bb([m1], f_donor=0.5)
    ll2 = total_log_likelihood_bb([m2], f_donor=0.5)
    ll_combined = total_log_likelihood_bb([m1, m2], f_donor=0.5, gq_weights=[1.0, 0.5])
    assert math.isclose(ll_combined, ll1 + 0.5 * ll2, rel_tol=1e-12)


def test_gq_weights_length_mismatch_raises():
    ms = _make_markers_for_fraction(f_donor=0.3, n_markers=3, dp=500)
    with pytest.raises(ValueError):
        total_log_likelihood_bb(ms, f_donor=0.3, gq_weights=[1.0, 0.5])
```

### `tests/test_multidonor.py`

Mirror tests (1)–(5) above for `total_log_likelihood_multi_bb`.

The multi-donor integration test needs a scenario strong enough that
the CI widening is reliably detectable (not the GQ=15 draft, which
gives only ~9% per-contributor downweight and doesn't reliably move
the profile CI). Use an extreme scenario:

- Take a sibling-donor fixture (`tests/test_data/multidonor/*.vcf`).
- Construct two copies of the `InformativeMarker` list: `M_full` with
  all GQs at 99, and `M_weak` with the first donor's GQ lowered to
  **11** (just above the floor) on **every marker informative for
  donor 1**. At GQ=11 the posterior is `1/(1+10^-1.1) ≈ 0.926`, and
  concentrated on one donor means roughly halving the effective
  information for that donor's fraction.
- Run `estimate_multi_donor(M_full, use_gq_weight=True)` vs
  `estimate_multi_donor(M_weak, use_gq_weight=True)`.
- Assert: donor-1 point estimate within 1.5% between the two runs;
  donor-1 CI *width* is at least 1.5× in the weak run; donor-2 point
  estimate within 0.5% and CI widens by < 1.2×.

If the test is still too flaky at GQ=11, drop to GQ=10.5 (via
`gq_floor=10` and a float-ish GQ cast) or simply pass precomputed
`gq_weights` of 0.5 for donor-1-informative markers. The goal is a
deterministic assertion on widening, not a real-GQ test.

### `tests/test_cli.py`

- Smoke test `allomix monitor --gq-weighted` end-to-end on one of the
  synthetic joint VCFs. Assert the command returns 0 and produces a
  valid TSV.
- Default (no `--gq-weighted`) produces bit-identical output to the
  pre-Step-16 baseline, confirming no silent behaviour change when
  off. Pin with a stored-reference TSV checked into `tests/test_data/`
  if convenient, or just compare two runs in the same test.
- Warning emission: when both `--min-gq 30` (non-default) and
  `--gq-weighted` are passed, stderr contains the "ignored" note.

### Regression

Full suite:

```bash
pytest -x -q
```

261+ tests must remain green. New GQ fields on `InformativeMarker`
default to `None`, so the fixtures in `test_qc.py:305`,
`test_multidonor.py:66`, and `test_chimerism.py:40` continue to work
unchanged.

## Verification plan

Ordered cheapest to most expensive.

### 1. Unit tests

```bash
.venv/bin/pytest tests/test_gq_weight.py tests/test_chimerism.py \
    tests/test_multidonor.py tests/test_cli.py -x -q
```

### 2. Synthetic CLI sanity (GQ=99 no-op check)

The existing `tests/test_data/multidonor/*.vcf` fixtures hardcode
GQ=99 (confirmed at `src/allomix/simulate.py:681` and
`src/allomix/simulate.py:841`). Step 1 of verification: confirm
`--gq-weighted` is an effective no-op on those fixtures, i.e. the
flag is wired up and doesn't break the run, but doesn't change the
numbers either because every `w_gq ≈ 1`:

```bash
.venv/bin/allomix monitor \
    --vcf tests/test_data/multidonor/sample_f1_20_f2_10.vcf \
    --host-sample HOST --donor-sample DONOR1 --donor-sample DONOR2 \
    --sample ADMIX --format json > /tmp/no_gq.json

.venv/bin/allomix monitor \
    --vcf tests/test_data/multidonor/sample_f1_20_f2_10.vcf \
    --host-sample HOST --donor-sample DONOR1 --donor-sample DONOR2 \
    --sample ADMIX --gq-weighted --format json > /tmp/with_gq.json

diff <(jq '.donor_fractions' /tmp/no_gq.json) \
     <(jq '.donor_fractions' /tmp/with_gq.json)
# Expect: identical to ~4 decimal places
```

This is a *wiring* check only. It does not exercise the weighted path.

### 3. Synthetic GQ-perturbation check (actual weighted path)

Because fixtures are GQ=99 everywhere, the weighted code path isn't
exercised by (2). Write a short verification script that patches the
GQ values in a fixture VCF and re-runs monitor:

```python
# scripts/verify_gq_weighting.py (single-use utility, do not commit
# unless we later want a persistent validation target).
import gzip, random, subprocess, sys
from pathlib import Path

IN = Path("tests/test_data/multidonor/sample_f1_20_f2_10.vcf")
OUT = Path("/tmp/sample_low_gq.vcf")

rng = random.Random(42)
with IN.open() as f_in, OUT.open("w") as f_out:
    for line in f_in:
        if line.startswith("#") or "\t" not in line:
            f_out.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        fmt = fields[8].split(":")
        gq_idx = fmt.index("GQ") if "GQ" in fmt else None
        if gq_idx is None:
            f_out.write(line)
            continue
        # Replace host (sample 9) and donor1 (sample 10) GQs with a
        # truncated-normal sample in [10, 99]; leave donor2 and admix
        # at 99 to isolate the effect on donor-1 inference.
        for s_idx in (9, 10):
            parts = fields[s_idx].split(":")
            parts[gq_idx] = str(max(10, min(99, int(rng.gauss(25, 10)))))
            fields[s_idx] = ":".join(parts)
        f_out.write("\t".join(fields) + "\n")

for label, args in [("no_gq", []), ("with_gq", ["--gq-weighted"])]:
    subprocess.run(
        ["allomix", "monitor", "--vcf", str(OUT),
         "--host-sample", "HOST", "--donor-sample", "DONOR1",
         "--donor-sample", "DONOR2", "--sample", "ADMIX",
         "--format", "json", "-o", f"/tmp/{label}.json", *args],
        check=True,
    )
```

(Sample column indices 9/10/11/12 above match the 4-sample fixture
layout; confirm with `grep '^#CHROM' tests/test_data/multidonor/sample_f1_20_f2_10.vcf`
before running.)

Then compare `/tmp/no_gq.json` vs `/tmp/with_gq.json` — the donor-1
fraction should shift slightly and its CI should widen noticeably
(donor-1 is the contributor whose GTs are now uncertain). Donor-2
should barely move. This is the smallest-possible real exercise of
the weighted code path without waiting for step (4) or changing the
simulator.

### 4. Optional: simulator extension for GQ variation

For paper-ready validation rather than smoke testing, add a
`--gq-mean` / `--gq-sd` knob to `scripts/generate_test_data.py` that
writes per-sample GQs from a truncated distribution (e.g. lognormal
mean 30 SD 15, clipped at 10). Feed that into a depth-validation-style
sweep:

```bash
.venv/bin/python scripts/generate_test_data.py \
    --output-dir output/gq_sim --depth 500 --gq-mean 30 --gq-sd 15
.venv/bin/allomix monitor \
    --vcf output/gq_sim/admix_f0.10.vcf \
    --host-sample HOST --donor-sample DONOR --sample ADMIX \
    --min-gq 20 --format json > /tmp/hardcut.json

.venv/bin/allomix monitor \
    --vcf output/gq_sim/admix_f0.10.vcf \
    --host-sample HOST --donor-sample DONOR --sample ADMIX \
    --gq-weighted --min-gq 20 --format json > /tmp/soft.json
```

Expectation: `soft` retains a handful more markers (GQ 10–19 recovered
at weight <1) and shows a modestly wider or slightly tighter CI
depending on which direction those markers point the MLE. Point
estimate shift <1 %.

### 5. Real-data smoke test

Use the April-24 validation batch (`output/validation_run_new_bias2/batch.tsv`)
as baseline. Re-run with `--gq-weighted`:

```bash
python scripts/run_xls_batch.py 'output/Chimerism project patient list.xlsx' \
    --vcf output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz \
    --host-column="NGS Sample ID" \
    --donor-column="NGS sample ID TP1" \
    --test-sample-column="NGS sample ID TP2" \
    --output-dir output/validation_run_gq_weighted \
    --copy-columns="Donor,Chimerism result TP2" \
    --bias-table-tsv output/bias_training/bias_table.tsv \
    --gq-weighted
```

`scripts/run_xls_batch.py` will need a thin `--gq-weighted` passthrough
(one kwarg, one subprocess-arg splice). Flag as a downstream follow-up
if this script is not being touched as part of this plan — it's small
and not strictly required for the core implementation.

Compare `donor_pct` and `ci_width` columns between `validation_run_new_bias2/batch.tsv`
and `validation_run_gq_weighted/batch.tsv`. Expectation: mean |delta
donor_pct| < 0.5%, mostly dominated by samples with low host-genotyping
GQs. No samples should flip PASS ↔ FAIL.

### 6. Paper benchmarks

Optional. `paper/scripts/run_depth_validation.py` at very low depth
(e.g. 50x) produces lower GQs in GATK-called data; a `--gq-weighted`
variant there might show a measurable MAE improvement at 50x where the
hard `min_gq` cut currently loses markers. Out of scope for the v1
implementation — add as a follow-up if early real-data results
(step 5) justify a paper figure.

## Edge cases and risks

- **VCFs without a GQ FORMAT field.** `parse_vcf` already returns
  `gq=None` in that case (`src/allomix/genotype.py:127`). `gq_weight`
  treats `None` as 1.0 per contributor, so `use_gq_weight=True` on a
  GQ-less VCF is a no-op. Safe default, no silent error.
- **Classification/runtime floor are the same under the new design.**
  With `effective_min_gq = args.gq_floor` when `--gq-weighted` is on,
  `classify_markers` filters below `gq_floor` and the runtime
  `gq_weight` also treats `gq < gq_floor` as weight zero — so the two
  never disagree. The legacy "double floor" risk from my earlier draft
  is removed by the `--min-gq` replacement semantics in §3b.
- **Pseudo-likelihood CI coverage.** Multiplying per-marker
  log-likelihoods by weights in [0, 1] changes the profile likelihood
  from a true likelihood to a pseudo-likelihood. The `chi2.ppf(0.95,
  df=1)/2` threshold is exact only asymptotically under the true
  model; under weighting, expected coverage drifts below nominal in
  proportion to the weight variance. This plan's §5c weighted GoF
  partially mitigates the diagnostic side (so a bad fit still looks
  bad), but nominal-95% coverage is no longer guaranteed. We already
  tolerate this from bias correction and will tolerate it again from
  dropout weighting, but the *compounded* effect of three approximate
  likelihoods (bias + dropout + GQ) is worth spot-checking against
  the empirical CI-coverage experiments (see `paper/scripts/run_depth_validation.py`)
  once steps 14-16 land. If coverage craters, revisit: bootstrap CIs
  or drop one of the weight mechanisms. Document explicitly in
  `paper/discussion.md`.
- **Interaction with Step 15 dropout weighting.** Anticipated in the
  Step 15 plan: the combined per-marker weight is the product
  `(1 - d_s) * w_gq(m)`. If both Step 15 and Step 16 are on, both
  kwargs are passed to `total_log_likelihood_bb`, and the function
  multiplies each LL by both factors independently. No ordering
  concern, no interaction to unit-test beyond "both on" smoke.
- **Interaction with Step 14 empirical error rates.** Orthogonal —
  Step 14 changes the `p_alt` computation inside
  `log_likelihood_marker_bb`; Step 16 scales the whole per-marker LL
  *after* it's computed. They compose without interference.
- **Per-donor CI widening for one donor.** When one donor's GT is
  low-confidence across many markers, the marker LLs that *inform*
  that donor shrink toward zero — effectively less data on that
  donor. The profile-likelihood CI correctly widens. This is the
  right behaviour but may surprise users used to tight CIs; call it
  out in the release note.
- **Numerical underflow at extreme GQ.** The Phred formula uses
  `10 ** (-gq / 10)`. At GQ 99, `10**-9.9 ≈ 1.26e-10`; at GQ 200,
  `1e-20`. No underflow risk in single precision; GQ rarely exceeds
  99 anyway (GATK caps there).
- **Fixture compatibility.** Every existing test that constructs
  `InformativeMarker` directly relies on the default values of
  unmentioned fields. Adding two more fields with defaults of `None`
  keeps those constructions working unchanged. Verified by inspection
  of `tests/test_chimerism.py:40`, `tests/test_multidonor.py:66`,
  `tests/test_multidonor.py:693`, `tests/test_bias.py:58`,
  `tests/test_qc.py:305`.
- **Simulator consistency.** `src/allomix/simulate.py` currently
  does not vary GQ. If we later add a `--gq-sd` knob (verification
  step 3), tests that rely on "GQ=99 everywhere" remain valid because
  the default stays GQ=99. Document in the simulator's docstring.

## Out of scope (follow-ups)

- Paper ablation: add a "no GQ weighting" baseline to Figure S4
  alongside "no bias" / "no overdispersion". Revisit after step-5
  real-data results.
- `scripts/run_xls_batch.py` — `--gq-weighted` passthrough flag.
- `scripts/generate_test_data.py` — `--gq-mean` / `--gq-sd` knobs to
  drive variation-in-GQ experiments.
- Genotype-miscall mixture (Modelling option C) if the
  pseudo-likelihood weighting underperforms on real data.

## File-by-file checklist

- [ ] `src/allomix/gq_weight.py` — NEW. `gq_weight(host_gq,
      donor_gqs, gq_floor, scheme, gq_full)` helper covering the
      Phred (proper posterior) and linear schemes.
- [ ] `src/allomix/genotype.py` — add `host_gq` and `donor_gqs`
      optional fields to `InformativeMarker` (lines 32-47); populate
      them in `classify_markers` at the `InformativeMarker(...)`
      construction site (lines 287-302). No change to `parse_vcf` or
      to the hard-floor GQ filter (lines 263-270).
- [ ] `src/allomix/chimerism.py` —
      (a) add `gq_weights: list[float] | None = None` kwarg to
      `total_log_likelihood_bb` (line 185) and
      `total_log_likelihood_multi_bb` (line 214); multiply per-marker
      LL by `gq_weights[i]` when provided; raise on length mismatch.
      (b) add user-facing kwargs `use_gq_weight`, `gq_floor`,
      `gq_weight_scheme` to `estimate_single_donor_bb` (line 341) and
      `estimate_multi_donor` (line 476); precompute the `gq_weights`
      list once at the top of each estimator and thread it through to
      every `total_log_likelihood_*` call site.
      (c) add the same `gq_weights` kwarg to
      `_profile_likelihood_cis_multi` (line 605) and forward it into
      the two inner optimisers.
      (d) extend `_compute_per_marker_results` and
      `_per_marker_results_multi` to accept `gq_weights` and set
      `MarkerResult.gq_weight` per marker (residual values themselves
      stay unweighted; 3-SD flag stays on unweighted residuals).
      Import `from allomix.gq_weight import gq_weight`. Do not modify
      `log_likelihood_marker_bb`.
- [ ] `src/allomix/qc.py` —
      (a) add `n_markers_effective: float` to `QCReport`; populate
      from `sum(mr.gq_weight for mr in result.per_marker if mr.included)`.
      (b) extend `_compute_gof_pval` to accept `gq_weights` (derived
      from `result.per_marker`) and compute weighted chi-squared with
      float `df = max(sum(gq_weights) - n_fitted_params, 1.0)`.
      Legacy path (`gq_weights` None) unchanged.
- [ ] `src/allomix/report.py` — add `gq_weight` column to verbose
      per-marker TSV; add `n_markers_effective` to summary TSV /
      JSON when present.
- [ ] `src/allomix/cli.py` — add `--gq-weighted`,
      `--gq-weight-scheme`, `--gq-floor` to `_add_common_args`
      (line 17); compute `effective_min_gq = args.gq_floor` when
      `--gq-weighted`, else `args.min_gq`, at the top of both
      `cmd_monitor` and `cmd_timeline` with a stderr note if both
      `--min-gq` (non-default) and `--gq-weighted` are set; extend
      `_run_single_sample` (line 73) with the three kwargs; pass them
      through estimator calls.
- [ ] `tests/test_gq_weight.py` — NEW. Phred (using proper posterior)
      and linear schemes, floor behaviour, None handling, invalid
      scheme raises.
- [ ] `tests/test_chimerism.py` — 5 aggregator tests on the
      `gq_weights` kwarg (None == all-ones exact, zero-weight exclusion
      exact, subset-equivalence exact, linearity to `rel_tol=1e-12`,
      length mismatch raises) plus one end-to-end estimator regression
      (`use_gq_weight=True` at all-GQ=99 matches `False` to
      `rel_tol=1e-6`).
- [ ] `tests/test_multidonor.py` — mirror the 5 aggregator cases on
      `total_log_likelihood_multi_bb`; add the extreme-GQ end-to-end
      test (GQ=11 on donor-1-informative markers, assert donor-1 CI
      width ≥ 1.5× baseline, donor-2 CI widens < 1.2×).
- [ ] `tests/test_cli.py` — smoke test `allomix monitor
      --gq-weighted`; assert bit-identical output when the flag is
      off vs unset; smoke test the `--min-gq` ignored warning when
      both are set.
- [ ] `tests/test_qc.py` — weighted-GoF regression (`gq_weights=None`
      matches current behaviour exactly); one test with heterogeneous
      weights showing lower chi-squared than the unweighted version
      when low-weight markers have large residuals.
- [ ] README — add `--gq-weighted` example block after the bias /
      dropout examples.
- [ ] `paper/methods.md` — one paragraph on the Phred posterior
      per-marker weight; cite as the natural generalisation of the
      `--min-gq` cutoff. One sentence in `paper/discussion.md` on
      compounded pseudo-likelihood coverage across bias + dropout +
      GQ weighting.
