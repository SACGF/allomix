# Plan: Empirical Per-Site Error Rates

Status: proposed, not implemented.

## TL;DR

Replace the global `--error-rate 0.01` constant with per-site error rates measured from the bias-training cohort. At hom-ref sites the observed ALT-read rate is the direct empirical estimate of `P(observe ALT | true REF)` (call this `e_refalt`); at hom-alt sites the observed REF-read rate gives `e_altref`. The two are not generally equal, so store both and switch the per-marker likelihood to an asymmetric form. Falls back to the existing global `--error-rate` when a site is missing from the table or not enough observations exist for a confident per-site estimate.

The architecture mirrors `bias.py` / `--bias-table` exactly: a new module `error_rates.py`, a new `allomix estimate-errors` subcommand, a new `--error-table` flag on `monitor` / `timeline`, and a new `marker_errors` parameter threaded through the chimerism estimators.

## Background

The current 4-state symmetric error model (`src/allomix/chimerism.py:159-165`) computes:

```python
p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
p_ref = w * (1.0 - e) + (1.0 - w) * e / 3.0
p_alt = p_alt / (p_ref + p_alt)
```

`e` is supplied by the caller (default `0.01`). Conceptually: an error happens with rate `e`, and when it does, it is uniform over the 3 other bases, so per-direction substitution is `e/3`. After conditional renormalisation onto `{REF, ALT}` (since N-base or third-allele observations are not represented in AD), the effective per-direction rate is approximately `(e/3) / (1 - 2e/3) ≈ e/3` for small `e`.

This model has two limitations:

1. **It is a single global constant.** Real per-site error rates vary by orders of magnitude — flanking sequence, repetitive context, strand bias and chemistry-specific substitution rates all matter. A site with mean ALT-rate at hom-ref of 5x10⁻⁴ is treated identically to a site with rate 5x10⁻³, both clamped to the global `e/3 = 3.3x10⁻³`. Likelihoods at low donor fractions (where the signal is in the per-marker error tail) lose information.
2. **It assumes symmetry.** Real sequencing-error rates differ between REF→ALT and ALT→REF, especially in oxidative-damage-driven contexts (`G→T`, `C→A`). The 4-state model cannot represent this.

The bias-training cohort (`/tau/data/clinical_hg38/idt_rhampseq_sid/` — 18,047 samples across 210 joint-called VCFs) gives both quantities directly: at every hom-ref call, observed `ad_alt/dp` is `e_refalt` for that site; at every hom-alt call, observed `ad_ref/dp` is `e_altref`. Pool across the cohort, store per site, plug into the likelihood. Should reduce dependence on a hand-tuned constant and improve fits at low-fraction samples.

References already in the codebase (informative, not to change):

- `src/allomix/bias.py` — analogous estimator/IO pattern this plan mirrors.
- `src/allomix/chimerism.py:133-182` — current likelihood with symmetric `e`.
- `src/allomix/simulate.py:280-307` — simulator's symmetric error model (must be kept consistent or updated alongside; see Out of Scope).

## Asymmetric error model

Drop the 4-state hack and model REF/ALT directly. With weight `w = P(read sampled from a REF allele)`:

```
p_alt = w * e_refalt + (1 - w) * (1 - e_altref)
p_ref = w * (1 - e_refalt) + (1 - w) * e_altref
```

These sum to 1 exactly, no renormalisation needed. Endpoints:

- `w = 1` (hom-ref): `p_alt = e_refalt`, `p_ref = 1 - e_refalt`.
- `w = 0` (hom-alt): `p_alt = 1 - e_altref`, `p_ref = e_altref`.
- `w = 0.5` (het): `p_alt = 0.5 + 0.5 · (e_refalt - e_altref)`.

Backward compatibility with the symmetric path: the legacy 4-state model with rate `e` is equivalent (for small `e`) to the asymmetric model with `e_refalt = e_altref = e/3`. So when no error table is provided, fall through to the existing 4-state code path unchanged. When an error table is provided, use the asymmetric form.

## Estimator design

For each marker key `(chrom, pos, ref, alt)` in the training cohort, accumulate four counters:

- `n_alt_at_homref` — sum of `ad_alt` across samples where `gt == (0, 0)` and quality filters pass
- `n_total_at_homref` — sum of `ad_ref + ad_alt` across the same samples
- `n_ref_at_homalt` — sum of `ad_ref` across samples where `gt == (1, 1)`
- `n_total_at_homalt` — sum of `ad_ref + ad_alt` across the same samples

Pooled MLE (one read = one trial):

```
e_refalt = n_alt_at_homref / n_total_at_homref
e_altref = n_ref_at_homalt / n_total_at_homalt
```

This is read-pooled, not sample-averaged. Sample-averaging would give equal weight to a 50x sample and a 5000x sample; read-pooling correctly weights by depth.

### Quality filters at training time

Borderline het calls miscalled as hom-ref/hom-alt would inflate the rate by orders of magnitude. Two safeguards:

1. **Use `min_gq`** when parsing training VCFs (default 20, matching the rest of the tool).
2. **Per-sample sanity check**: if `ad_alt / dp > 0.10` at a hom-ref call, drop that observation; if `ad_ref / dp > 0.10` at a hom-alt call, drop. The 0.10 threshold is well above any realistic error rate (>30σ at typical depth) so anything above it is almost certainly a miscalled het, contamination, or somatic event. Configurable via `--max-vaf-homref` / `--min-vaf-homalt`.

### Minimum-observations filter

Sites with too few homozygous training observations have noisy rate estimates. Default `--min-reads 1000` (configurable). Sites below threshold are omitted from the output; the runtime `monitor` will fall through to the global `--error-rate` for those sites.

Per direction: a site with abundant hom-ref but no hom-alt (e.g. rare ALT allele) gets `e_refalt` populated and `e_altref = NA`. The runtime should handle each direction independently — at a hom-ref call site we only need `e_refalt`; at a hom-alt call site only `e_altref`. Het sites need both. Storing NA per direction (rather than dropping the row) preserves partial information.

## Fix scope

Five files, plus tests. New module + new CLI subcommand + plumbing through the estimators.

### 1. `src/allomix/error_rates.py` (NEW)

```python
"""Per-site empirical sequencing error rate estimation.

Estimates per-marker, per-direction substitution rates from a training
cohort of joint-called VCFs. At hom-ref calls the observed ALT-read
rate estimates ``P(observe ALT | true REF)``; at hom-alt calls the
observed REF-read rate estimates ``P(observe REF | true ALT)``.

The two are not generally equal: oxidation damage, strand bias, and
flanking sequence context all produce direction-specific error rates.
The output table is consumed by ``chimerism.estimate_single_donor_bb``
and ``chimerism.estimate_multi_donor`` via the ``marker_errors``
parameter.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from allomix.genotype import MarkerData

MarkerKey = tuple[str, int, str, str]


@dataclass
class MarkerError:
    """Per-marker, per-direction empirical error rates."""

    chrom: str
    pos: int
    ref: str
    alt: str
    e_refalt: float | None  # ALT-read rate at hom-ref calls; None if no data
    e_altref: float | None  # REF-read rate at hom-alt calls; None if no data
    n_reads_homref: int
    n_reads_homalt: int


def _marker_key(m: MarkerData) -> MarkerKey:
    return (m.chrom, m.pos, m.ref, m.alt)


def estimate_error_rates(
    marker_lists: list[list[MarkerData]],
    min_reads: int = 1000,
    max_vaf_homref: float = 0.10,
    min_vaf_homalt: float = 0.90,
) -> dict[MarkerKey, MarkerError]:
    """Estimate per-marker, per-direction error rates from training samples.

    Pooled across reads (not averaged across samples), so high-depth
    samples carry more weight, which is the correct MLE under the
    assumption that all reads at a site share the same per-direction
    error rate.

    Args:
        marker_lists: List of MarkerData lists, one per training sample.
            Apply ``min_gq`` at parse time (e.g. ``parse_vcf(..., min_gq=20)``)
            to exclude low-confidence calls.
        min_reads: Minimum total reads required *per direction* to retain
            a site's estimate. Sites with fewer reads in a direction get
            ``None`` for that direction's rate (the runtime falls through
            to the global ``--error-rate``). Default 1000.
        max_vaf_homref: Drop hom-ref observations where ``ad_alt/dp``
            exceeds this threshold. Protects against miscalled hets and
            contamination inflating the rate. Default 0.10.
        min_vaf_homalt: Drop hom-alt observations where ``ad_ref/dp``
            falls below this threshold (i.e. ``vaf < min_vaf_homalt``).
            Default 0.90.

    Returns:
        Dict mapping (chrom, pos, ref, alt) to MarkerError. Sites with no
        usable observations in either direction are omitted.
    """
    # Per-site accumulators
    n_alt_homref: dict[MarkerKey, int] = {}
    n_tot_homref: dict[MarkerKey, int] = {}
    n_ref_homalt: dict[MarkerKey, int] = {}
    n_tot_homalt: dict[MarkerKey, int] = {}
    info: dict[MarkerKey, tuple[str, int, str, str]] = {}

    for markers in marker_lists:
        for m in markers:
            dp = m.ad_ref + m.ad_alt
            if dp <= 0:
                continue
            key = _marker_key(m)
            info[key] = (m.chrom, m.pos, m.ref, m.alt)
            if m.gt == (0, 0):
                vaf = m.ad_alt / dp
                if vaf > max_vaf_homref:
                    continue
                n_alt_homref[key] = n_alt_homref.get(key, 0) + m.ad_alt
                n_tot_homref[key] = n_tot_homref.get(key, 0) + dp
            elif m.gt == (1, 1):
                vaf = m.ad_alt / dp
                if vaf < min_vaf_homalt:
                    continue
                n_ref_homalt[key] = n_ref_homalt.get(key, 0) + m.ad_ref
                n_tot_homalt[key] = n_tot_homalt.get(key, 0) + dp
            # Hets are ignored — they're used by bias estimation, not error.

    out: dict[MarkerKey, MarkerError] = {}
    for key, (chrom, pos, ref, alt) in info.items():
        nh_tot = n_tot_homref.get(key, 0)
        na_tot = n_tot_homalt.get(key, 0)
        e_ra: float | None = None
        e_ar: float | None = None
        if nh_tot >= min_reads:
            e_ra = n_alt_homref[key] / nh_tot
        if na_tot >= min_reads:
            e_ar = n_ref_homalt[key] / na_tot
        if e_ra is None and e_ar is None:
            continue
        out[key] = MarkerError(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            e_refalt=e_ra,
            e_altref=e_ar,
            n_reads_homref=nh_tot,
            n_reads_homalt=na_tot,
        )
    return out


def save_error_table(errors: dict[MarkerKey, MarkerError], path: Path | str) -> None:
    """Write error-rate estimates to a TSV file.

    Format:

        chrom  pos  ref  alt  e_refalt  e_altref  n_reads_homref  n_reads_homalt

    Missing per-direction rates are written as ``NA``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["chrom", "pos", "ref", "alt", "e_refalt", "e_altref",
             "n_reads_homref", "n_reads_homalt"]
        )
        for key in sorted(errors.keys()):
            me = errors[key]
            w.writerow(
                [
                    me.chrom,
                    me.pos,
                    me.ref,
                    me.alt,
                    "NA" if me.e_refalt is None else f"{me.e_refalt:.6e}",
                    "NA" if me.e_altref is None else f"{me.e_altref:.6e}",
                    me.n_reads_homref,
                    me.n_reads_homalt,
                ]
            )


def load_error_table(path: Path | str) -> dict[MarkerKey, tuple[float | None, float | None]]:
    """Load an error-rate table. Returns key -> (e_refalt, e_altref).

    NA entries become ``None``.
    """
    out: dict[MarkerKey, tuple[float | None, float | None]] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key: MarkerKey = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            e_ra = None if row["e_refalt"] == "NA" else float(row["e_refalt"])
            e_ar = None if row["e_altref"] == "NA" else float(row["e_altref"])
            out[key] = (e_ra, e_ar)
    return out
```

### 2. `src/allomix/chimerism.py`

#### `log_likelihood_marker_bb` — extend to accept asymmetric rates

Current (lines 133-182):

```python
def log_likelihood_marker_bb(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
) -> float:
    e = error_rate
    p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
    p_ref = w * (1.0 - e) + (1.0 - w) * e / 3.0
    p_alt = p_alt / (p_ref + p_alt)
    p_alt = max(1e-6, min(1.0 - 1e-6, p_alt))
    ...
```

Proposed:

```python
def log_likelihood_marker_bb(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    e_refalt: float | None = None,
    e_altref: float | None = None,
) -> float:
    """Per-marker beta-binomial log-likelihood.

    When ``e_refalt`` and ``e_altref`` are both supplied, uses the
    asymmetric REF/ALT-only error model:

        p_alt = w * e_refalt + (1 - w) * (1 - e_altref)

    Otherwise falls back to the legacy 4-state symmetric model with rate
    ``error_rate``. Per-direction asymmetric rates may come from
    ``error_rates.estimate_error_rates``. Either rate may be ``None``
    individually; the legacy fallback is used in that case as well.

    Args:
        ad_ref: REF read count.
        ad_alt: ALT read count.
        w: Expected reference allele weight (after bias correction).
        error_rate: Symmetric 4-state rate, used only when asymmetric
            rates are not both provided.
        rho: Beta-binomial concentration.
        e_refalt: ``P(observe ALT | true REF base)``. Empirical per-site.
        e_altref: ``P(observe REF | true ALT base)``. Empirical per-site.
    """
    if e_refalt is not None and e_altref is not None:
        p_alt = w * e_refalt + (1.0 - w) * (1.0 - e_altref)
    else:
        e = error_rate
        p_alt_raw = (1.0 - w) * (1.0 - e) + w * e / 3.0
        p_ref_raw = w * (1.0 - e) + (1.0 - w) * e / 3.0
        p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)

    p_alt = max(1e-6, min(1.0 - 1e-6, p_alt))

    n = ad_ref + ad_alt
    k = ad_alt
    if n == 0:
        return 0.0

    a = max(p_alt * rho, 1e-10)
    b = max((1.0 - p_alt) * rho, 1e-10)
    return (
        lgamma(k + a) + lgamma(n - k + b) - lgamma(n + rho)
        - lgamma(a) - lgamma(b) + lgamma(rho)
    )
```

#### `total_log_likelihood_bb` — accept and dispatch a per-marker error table

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
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> float:
    ll = 0.0
    for m in markers:
        key = (m.chrom, m.pos, m.ref, m.alt)
        bias = marker_biases.get(key, 0.0) if marker_biases is not None else 0.0
        e_ra: float | None = None
        e_ar: float | None = None
        if marker_errors is not None:
            erra = marker_errors.get(key)
            if erra is not None:
                e_ra, e_ar = erra
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        ll += log_likelihood_marker_bb(
            m.admix_ad_ref,
            m.admix_ad_alt,
            w,
            error_rate=error_rate,
            rho=rho,
            e_refalt=e_ra,
            e_altref=e_ar,
        )
    return ll
```

`total_log_likelihood_multi_bb` (lines 214-240) gets the same `marker_errors` parameter and same dispatch — same code shape, just the multi-donor `expected_weight_multi` call instead.

#### `_compute_per_marker_results` — propagate `marker_errors` for residual / GoF

Current (lines 278-338) does not need the error table for the *residual* itself (residual is observed_vaf − expected_vaf, no probability scaling). But to keep `expected_vaf` consistent with what the likelihood scored, no change is needed unless future GoF wants per-marker error-aware variance. Leave for now; flagged as follow-up.

#### `estimate_single_donor_bb` — accept and thread `marker_errors`

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
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> ChimerismResult:
```

Inside the function, every call to `total_log_likelihood_bb(...)` (lines 386, 406, 428) needs the new `marker_errors=marker_errors` argument:

```python
# line 384-390 (grid search)
opt_rho = minimize_scalar(
    lambda log_r: (
        -total_log_likelihood_bb(
            markers, f, error_rate, math.exp(log_r), marker_biases, marker_errors
        )
    ),
    bounds=(math.log(1.0), math.log(10000.0)),
    method="bounded",
)

# line 399-406 (Nelder-Mead joint refinement)
def neg_ll_joint(x):
    f_val, log_rho_val = x
    if f_val < 0.0 or f_val > 1.0:
        return 1e30
    rho_val = math.exp(log_rho_val)
    if rho_val < 0.5 or rho_val > 50000:
        return 1e30
    return -total_log_likelihood_bb(
        markers, f_val, error_rate, rho_val, marker_biases, marker_errors
    )

# line 424-433 (profile_ll_f for CIs)
def profile_ll_f(f_val: float) -> float:
    opt_rho = minimize_scalar(
        lambda log_r: (
            -total_log_likelihood_bb(
                markers, f_val, error_rate, math.exp(log_r), marker_biases, marker_errors
            )
        ),
        bounds=(math.log(1.0), math.log(50000.0)),
        method="bounded",
    )
    return -float(opt_rho.fun)
```

The `_compute_per_marker_results` call at line 455 does not need updating for v1 (residuals are model-agnostic).

#### `estimate_multi_donor` — analogous

Add `marker_errors` parameter; thread into:

- `total_log_likelihood_multi_bb` calls at lines 542 and 556.
- `_profile_likelihood_cis_multi` (line 573) — needs the same parameter added to that helper, and forwarded in lines 636 and 654.
- `_per_marker_results_multi` (line 576) — no change needed for v1.

### 3. `src/allomix/cli.py`

#### Add `--error-table` and `--no-error-correction` to common args

`_add_common_args` (lines 17-54) — append:

```python
parser.add_argument(
    "--error-table",
    default=None,
    help="Per-site empirical error-rate table TSV "
         "(from `allomix estimate-errors`). When provided, sites with "
         "per-direction rates override --error-rate; missing sites or "
         "missing directions fall back to --error-rate.",
)
parser.add_argument(
    "--no-error-correction",
    action="store_true",
    help="Disable empirical error-rate correction even when an error "
         "table is provided",
)
```

#### Add a `_load_errors` helper next to `_load_biases` (line 121)

```python
def _load_errors(args: argparse.Namespace):
    """Load per-site error table if specified and not disabled."""
    if args.error_table and not args.no_error_correction:
        return load_error_table(args.error_table)
    return None
```

Add the import at the top:

```python
from allomix.error_rates import (
    estimate_error_rates,
    load_error_table,
    save_error_table,
)
```

#### Thread `marker_errors` into `_run_single_sample`

Current (lines 73-111): add `marker_errors` parameter, pass to both `estimate_single_donor_bb` and `estimate_multi_donor`. Both `cmd_monitor` and `cmd_timeline` load the table once before the per-sample loop and pass it in (mirroring `marker_biases`).

```python
def _run_single_sample(
    host: list,
    donors: list[list],
    vcf_path: str,
    admix_sample: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    marker_biases=None,
    marker_errors=None,
) -> tuple:
    ...
    if len(donors) == 1:
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
            marker_errors=marker_errors,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            marker_biases=marker_biases,
            marker_errors=marker_errors,
        )
    ...
```

`cmd_monitor` and `cmd_timeline` each gain:

```python
marker_errors = _load_errors(args)
...
result, qc, genotypes = _run_single_sample(
    host, donors, args.vcf, sample_name,
    args.min_dp, args.min_gq, args.error_rate,
    marker_biases=marker_biases,
    marker_errors=marker_errors,
)
```

#### New `estimate-errors` subcommand

Mirror `cmd_estimate_bias` (lines 202-230) and `bias_parser` (lines 261-294). New function:

```python
def cmd_estimate_errors(args: argparse.Namespace) -> int:
    """Run the estimate-errors subcommand."""
    if args.vcfs and args.vcf:
        raise SystemExit("Use either --vcfs or --vcf/--samples, not both")
    if not args.vcfs and not args.vcf:
        raise SystemExit("One of --vcfs or --vcf is required")
    if args.vcf and not args.samples:
        raise SystemExit("--samples is required when using --vcf")

    marker_lists = []
    if args.vcfs:
        for vcf_path in args.vcfs:
            markers = parse_vcf(vcf_path, min_dp=0, min_gq=args.min_gq)
            marker_lists.append(markers)
        n_source = f"{len(args.vcfs)} VCFs"
    else:
        _validate_sample_names(args.vcf, args.samples)
        for sample in args.samples:
            markers = parse_vcf(args.vcf, sample=sample, min_dp=0, min_gq=args.min_gq)
            marker_lists.append(markers)
        n_source = f"{len(args.samples)} samples from {args.vcf}"

    errors = estimate_error_rates(
        marker_lists,
        min_reads=args.min_reads,
        max_vaf_homref=args.max_vaf_homref,
        min_vaf_homalt=args.min_vaf_homalt,
    )
    save_error_table(errors, args.output)
    print(
        f"Estimated error rates for {len(errors)} sites from {n_source} "
        f"-> {args.output}",
        file=sys.stderr,
    )
    return 0
```

And the parser block in `main()`:

```python
err_parser = subparsers.add_parser(
    "estimate-errors",
    help="Estimate per-site empirical error rates from VCFs",
)
err_input = err_parser.add_mutually_exclusive_group()
err_input.add_argument("--vcfs", nargs="+", metavar="VCF",
    help="Per-sample VCFs, one per file (reads first sample from each)")
err_input.add_argument("--vcf", metavar="VCF",
    help="Joint-called multi-sample VCF (use with --samples)")
err_parser.add_argument("--samples", nargs="+", metavar="SAMPLE_NAME",
    help="Sample names to extract from --vcf")
err_parser.add_argument("--output", "-o", default="error_table.tsv",
    help="Output error table TSV (default: error_table.tsv)")
err_parser.add_argument("--min-reads", type=int, default=1000,
    help="Minimum total reads per direction to retain a site's estimate "
         "(default: 1000)")
err_parser.add_argument("--max-vaf-homref", type=float, default=0.10,
    help="Drop hom-ref training observations with vaf > this (default: 0.10)")
err_parser.add_argument("--min-vaf-homalt", type=float, default=0.90,
    help="Drop hom-alt training observations with vaf < this (default: 0.90)")
err_parser.add_argument("--min-gq", type=int, default=20,
    help="Minimum GQ for training calls (default: 20)")
```

Plus the dispatch line in the `if args.command == ...` chain:

```python
if args.command == "estimate-errors":
    return cmd_estimate_errors(args)
```

### 4. `src/allomix/qc.py`

The existing `_compute_gof_pval` uses `result.error_rate` (single global value) for the variance floor (lines 121-125 in current code). With per-site errors it would be more accurate to use the per-site rates here too, but the GoF test is a coarse diagnostic and the beta-binomial `ρ` already absorbs site-to-site error variance. **No change for v1.** Mark as follow-up.

If reviewers push back: add an optional `marker_errors` parameter to `_compute_gof_pval` and use the per-site rate when computing `ev_raw` for the variance floor. Same fallback to the global rate when missing.

### 5. Tests

#### New file `tests/test_error_rates.py`

Mirror `tests/test_bias.py` structure. Cover:

- `TestEstimateErrorRates`
  - `test_no_homozygous_calls_returns_empty` — only het calls in input
  - `test_clean_homref_recovers_input_rate` — synthetic data with known per-direction error rate, check recovered rate within tolerance
  - `test_max_vaf_homref_drops_miscalled_het` — observation with `vaf > 0.10` at a hom-ref call is excluded
  - `test_min_reads_filter` — site with too few reads returns no estimate (or `None` per direction)
  - `test_per_direction_independence` — site with hom-ref obs but no hom-alt obs gets `e_refalt` populated and `e_altref=None`
  - `test_pooling_weights_by_depth` — two samples (depths 100 and 10000) at the same site should produce a pooled rate biased toward the deeper sample
- `TestSaveLoadRoundtrip`
  - Write a table, reload it, check identity; check `NA` round-trips to `None`
- `TestLikelihoodIntegration` (in `test_chimerism.py` extension or here)
  - Synthetic markers where `marker_errors` shifts the MLE in a predictable direction (e.g. asymmetric error suppresses ALT signal at one site, MLE at next-higher fraction)
  - Symmetric `e_refalt = e_altref = e/3` should approximately match the legacy 4-state model with rate `e` (loose tolerance, since the symmetric model also allocates probability to two non-REF/ALT bases)

#### Update `tests/test_chimerism.py` (no current breakage)

The `marker_errors` parameter has a default of `None`, so existing tests pass unchanged. Add at least one test exercising the asymmetric path:

```python
def test_asymmetric_error_overrides_symmetric():
    markers = [...]  # constructed with f=0.05 known truth
    me = {(m.chrom, m.pos, m.ref, m.alt): (1e-4, 5e-3) for m in markers}
    res_with = estimate_single_donor_bb(markers, marker_errors=me)
    res_without = estimate_single_donor_bb(markers)  # default e=0.01 -> e/3
    # Two should differ noticeably; assert finite + sane
    assert 0.0 <= res_with.donor_fraction <= 1.0
    assert res_with.log_likelihood > res_without.log_likelihood  # expected if synth uses asymmetric truth
```

Pick the synthetic data to make the inequality robust; back-of-envelope check before committing the assertion.

#### Update `tests/test_cli.py` if it covers `monitor` end-to-end

Add a smoke test that invokes `monitor` with `--error-table` pointing to a small fixture TSV and asserts a non-trivial result. Ensure `--no-error-correction` short-circuits loading.

## Verification plan

Cheapest first.

### 1. Unit-test sanity

```bash
.venv/bin/pytest tests/test_error_rates.py tests/test_chimerism.py tests/test_cli.py -x -q
```

Existing tests pass unchanged (default `marker_errors=None` keeps the legacy 4-state path). New tests pass.

### 2. Synthetic round-trip

Generate a training cohort with the simulator at known `error_rate=0.01` and confirm the estimator recovers the expected per-direction rate of `≈ 0.01 / 3 ≈ 3.3e-3` (the symmetric 4-state model assigns `e/3` to each per-direction substitution).

`scripts/generate_test_data.py` currently hardcodes `error_rate=0.01` (lines 144 and 162) and doesn't expose it as a CLI flag. Either temporarily edit it for the round-trip, add a `--error-rate` argument as part of this work, or write a short ad-hoc generation script that calls `simulate.sample_allele_counts(...)` directly. Recommended: add the CLI flag (one-line change) since it's useful beyond this verification.

```bash
# After adding --error-rate to generate_test_data.py:
.venv/bin/python scripts/generate_test_data.py --output-dir output/error_train \
    --error-rate 0.01 --depth 1000 --n-replicates 50

.venv/bin/allomix estimate-errors --vcfs output/error_train/*.vcf \
    -o output/error_train/error_table.tsv --min-reads 1000

.venv/bin/python -c "
from allomix.error_rates import load_error_table
import statistics
t = load_error_table('output/error_train/error_table.tsv')
ras = [v[0] for v in t.values() if v[0] is not None]
ars = [v[1] for v in t.values() if v[1] is not None]
print('e_refalt mean=', statistics.mean(ras), 'median=', statistics.median(ras))
print('e_altref mean=', statistics.mean(ars), 'median=', statistics.median(ars))
print('expected ≈ 0.0033 (= 0.01/3)')
"
```

### 3. End-to-end monitor with synthetic chimerism

Use `tests/test_data/joint_single_donor.vcf` (already exists). Estimate errors over the unrelated training cohort, then run monitor with and without `--error-table`. Compare donor fractions and CIs:

```bash
.venv/bin/allomix monitor --vcf tests/test_data/joint_single_donor.vcf \
    --host-sample HOST --donor-sample DONOR --sample ADMIX_F0.10 \
    --format tsv > /tmp/monitor_no_err.tsv

.venv/bin/allomix monitor --vcf tests/test_data/joint_single_donor.vcf \
    --host-sample HOST --donor-sample DONOR --sample ADMIX_F0.10 \
    --error-table output/error_train/error_table.tsv \
    --format tsv > /tmp/monitor_with_err.tsv

diff /tmp/monitor_no_err.tsv /tmp/monitor_with_err.tsv
```

Expectation: small differences in donor_pct (< 0.5%); CIs may shrink slightly at low fractions.

### 4. Real-data smoke test

Generate the empirical error table from the bias-training cohort (~50–200 unrelated samples from `/tau`). Then re-run the April 24 validation batch (`output/validation_run_new_bias/batch.tsv`) with `--error-table` and compare donor_pct deltas vs the existing run.

```bash
.venv/bin/allomix estimate-errors \
    --vcf output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz \
    --samples $(cat output/training_samples.txt) \
    -o output/error_training/error_table.tsv

python scripts/run_xls_batch.py output/Chimerism\ project\ patient\ list.xlsx \
    --vcf output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz \
    --host-column="NGS Sample ID" \
    --donor-column="NGS sample ID TP1" \
    --test-sample-column "NGS sample ID TP2" \
    --output-dir output/validation_run_per_site_err \
    --copy-columns="Donor,Chimerism result TP2" \
    --bias-table-tsv output/bias_training/bias_table.tsv \
    --error-table-tsv output/error_training/error_table.tsv
```

(Note: `scripts/run_xls_batch.py` will need a thin wrapper for the new flag — out of scope for this plan if not already supported, but small.)

Compare `donor_pct` columns between `validation_run_new_bias/batch.tsv` and `validation_run_per_site_err/batch.tsv`. Largest deltas should be at low-fraction samples.

### 5. Paper benchmarks

`paper/scripts/run_depth_validation.py` and `paper/scripts/compare_bias_correction.py` should be extended to optionally take an error table. Out of scope for the v1 implementation; add a paragraph to `paper/discussion.md` once empirical error rates are in.

## Edge cases and risks

- **Sites missing from training cohort**: Brand-new panel markers, or rare ALT alleles never observed homozygous, will be absent from the error table. The dispatch in `total_log_likelihood_bb` handles this — `marker_errors.get(key)` returns `None`, both `e_refalt` and `e_altref` stay `None`, and the legacy 4-state path is used. No crash, no silent miscalculation.
- **Only one direction populated**: A site may have abundant hom-ref observations but no hom-alt (rare ALT allele). The plan stores `(e_refalt, None)` and the dispatch falls through to legacy when *either* is `None`. **Decision needed**: alternative is to use `e_refalt` for hom-ref-call markers and global for hom-alt-call markers, mixing per-marker. Cleaner v1: require *both* directions; user can lower `min_reads` or fall back globally. Document in the table header comment.
- **Overconfidence at very low estimated rates**: A site with 10⁻⁶ observed rate (e.g. 0 ALT reads in 10⁵ depth at hom-ref) will give `e_refalt ≈ 0`, and if a real admixture sample has any ALT reads at that site under hom-ref expectation, it will have a likelihood near zero. Add a minimum floor (e.g. `e_refalt = max(observed, 1e-5)`) to prevent infinite log-likelihood penalties on a single noisy read. Implement in the loader (cleaner than baking into the estimator).
- **Het contamination at the training threshold**: `max_vaf_homref=0.10` is a hard cutoff. A real site where the true error rate is 8% (extreme but possible at long-read homopolymer ends) would have legitimate hom-ref observations rejected. The current panel doesn't include such regions, but worth noting as a follow-up to validate against real per-site rate distribution before fixing the threshold.
- **GT call quality dependence**: Per-site error rates depend on the variant caller's hom-ref/hom-alt classification. If GATK marginal-likelihood thresholds change between training and admixture runs, rates could drift. Mitigation: regenerate the table whenever the variant calling pipeline changes. Document.
- **Pooling assumption**: Pooling reads across samples assumes the per-direction error rate is shared across samples at a site. If a per-sample run effect dominates (e.g. one bad batch), pooled rates will be biased toward the bad batch in proportion to its read share. Mitigation: per-sample weighting variant in a follow-up; for v1 trust the cohort.
- **Test fixtures**: `tests/test_qc.py:54`, `tests/test_report.py:54`, `tests/test_multidonor.py:98` construct `ChimerismResult`/`MultiDonorResult` directly. The plan does not change those dataclasses — `marker_errors` is an *input* to the estimator, not stored on the result. So fixtures don't need updating.

## Out of scope (follow-ups to track separately)

- Update `simulate.py` to use a per-site error table for synthetic data generation. Necessary before the paper claims empirical error rates improve simulation realism.
- Update `_compute_gof_pval` (`qc.py`) to use per-site error rates for the variance floor. Small accuracy gain; defer until needed.
- Optional shrinkage of per-site rates toward a global prior when `n_reads` is moderate (e.g. empirical-Bayes Beta prior). Improves stability for sites with 1000–5000 reads. v1 uses a hard `min_reads` cutoff instead.
- Per-sample run-level error rate as a fallback (between per-site and global). Captures batch effects when a site has too few cohort observations but the current sample has enough total depth.
- Expose the per-site rates in `--verbose` output of `monitor` so users can see which sites used empirical vs fallback rates.
- Cross-validate per-site rates against MAF-stratified subsets — common ALT alleles will have abundant hom-alt obs for `e_altref`, rare ones won't. May want to report per-site coverage statistics alongside the rate.
- `scripts/run_xls_batch.py` extension to take `--error-table-tsv` (in line with the existing `--bias-table-tsv` option).
- Paper methods section: add a subsection on empirical per-site error rates after the bias-correction subsection in `paper/methods.md`.

## File-by-file checklist

- [ ] `src/allomix/error_rates.py`: NEW. `MarkerError` dataclass, `estimate_error_rates`, `save_error_table`, `load_error_table`.
- [ ] `src/allomix/chimerism.py`: extend `log_likelihood_marker_bb` with `e_refalt`, `e_altref`. Add `marker_errors` to `total_log_likelihood_bb`, `total_log_likelihood_multi_bb`, `estimate_single_donor_bb`, `estimate_multi_donor`, `_profile_likelihood_cis_multi`. Thread through every call site inside the estimators.
- [ ] `src/allomix/cli.py`: import from `allomix.error_rates`. Add `--error-table` and `--no-error-correction` to `_add_common_args`. Add `_load_errors` helper. Thread `marker_errors` into `_run_single_sample`, `cmd_monitor`, `cmd_timeline`. Add `cmd_estimate_errors` and the `estimate-errors` subparser.
- [ ] `tests/test_error_rates.py`: NEW. Cover estimator, save/load round-trip, likelihood integration.
- [ ] `tests/test_chimerism.py`: add at least one test exercising the asymmetric path.
- [ ] `tests/test_cli.py`: smoke test for `monitor --error-table` and `estimate-errors`.
- [ ] Run synthetic round-trip and real-data smoke test (Verification §2 and §4); record before/after donor_pct deltas.
