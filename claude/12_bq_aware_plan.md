# Plan: Per-Base Quality-Aware Likelihood (Step 12)

Status: proposed, not implemented.

## TL;DR

Today every read at every marker contributes to the likelihood with the same flat
sequencing error rate (`--error-rate`, default 0.01). Real base qualities vary a lot:
a typical Q30 base has e=1e-3, a Q20 base e=1e-2, and soft-clipped / low-complexity
bases run into Q10 or lower (e=0.1). Conpair's contamination estimator (Bergmann et al.
2016) sets a per-read error `e_i = 10^(-Q_i/10)` and takes the likelihood as a product
over reads. Replicating that approach in allomix should tighten CIs at borderline
samples (low depth, near detection limit, low donor fraction) and make the point
estimate less sensitive to the `--error-rate` constant.

**Complication:** the current production VCF has no per-base quality data. GATK
HaplotypeCaller output exposes neither per-read BQs nor per-sample aggregated BQ
sums — `ls` of `output/joint_called/` confirms only GT/AD/DP/GQ/PL/SB FORMAT
fields. **Solution: fix this at variant-calling time.** Add a one-line
`bcftools mpileup -a FORMAT/QS` + `bcftools annotate` step to the upstream
pipeline that folds `FORMAT/QS` (phred-summed BQ per allele per sample) into
the GATK joint-called VCF. The VCF remains the sole primary input to allomix;
no BAM fallback or sidecar table is needed. This keeps the VCF-first design
from Step 2 intact and makes the CLI plumbing trivial (one new flag).

### Why call the BAM twice (brief rationale)

It is a kludge, but it is the standard kludge. GATK HaplotypeCaller cannot emit
per-allele base quality as a FORMAT field. The `BaseQuality` annotation (FORMAT
tag `MBQ`) in GATK source explicitly implements `StandardMutectAnnotation`, so
it is wired to Mutect2 only; passing `-A BaseQuality` to HaplotypeCaller is
silently ignored (GATK docs: *"several annotations that are currently not
hooked up to HaplotypeCaller ... no error or warning message will be
provided"*). The alternatives are all worse:

- Switch to Mutect2 for calling — overkill, wrong tool (somatic caller with
  tumor-only caveats), invasive pipeline change.
- `GATK VariantAnnotator -A BaseQuality` against the HC VCF + BAM — undocumented
  for HC output since `BaseQuality` is a Mutect2 annotation class; worth a
  quick test but not a path to bet the plan on.
- Read BAMs directly from inside allomix — what Conpair/Demixtify/somalier do;
  brings pysam into the core package and contradicts the VCF-first design
  (Step 2) and CLAUDE.md's no-lazy-imports rule.

A second pass with `bcftools mpileup -a FORMAT/QS` restricted to the panel
sites (`-R panel.bed`) runs in seconds per BAM (76 sites, not whole-genome),
adds two Snakemake rules, and produces a well-specified, portable FORMAT tag
(`QS`, Number=R, bcftools docs). That is the right trade: a minor pipeline
addition in exchange for clean VCF-first downstream semantics.

## Background

Crysup & Woerner (2022, Formula 5) is the per-marker contribution to the total log
likelihood under a mixture genotype. Allomix currently evaluates it as a
beta-binomial on the aggregated (ad_ref, ad_alt) counts, with a 4-state error model
setting the mean:

```
p_alt = w(1-e) + (1-w) e/3        # chimerism.py:160
p_ref = (1-w)(1-e) + w e/3
```

where `w = (1-f) * host_ref_dose/2 + f * donor_ref_dose/2` and `e` is the global
flat error rate.

A per-read version replaces the aggregate with a product over reads, each with its
own e_i:

```
L_marker = ∏_{reads i}  [ I(read i is REF) * p_ref(e_i)
                         + I(read i is ALT) * p_alt(e_i) ]
```

where `p_ref(e_i) = w(1-e_i) + (1-w) e_i/3` and symmetric for p_alt. In practice
we don't need to know which *specific* reads were ALT and which REF; we only need
the two phred-summed quality totals:

```
Q_alt = sum of BQs at reads that observed the ALT base
Q_ref = sum of BQs at reads that observed the REF base
```

From these we can compute `e_alt_mean = 10^(-Q_alt/(10*ad_alt))`, i.e. the geometric
mean of per-read error rates (more precisely, the mean of phred scores). That
aggregated form is exactly what bcftools mpileup's `QS`/`QA`/`QR` tags carry.

Three levels of fidelity are possible, each cheaper than the last:

| Tier | What's modelled | Data needed | Estimate change |
|------|-----------------|-------------|----------------|
| T1 | Per-marker aggregated error rate `e_m = mean BQ across reads at marker m` | FORMAT QS or BAM pileup | Biggest improvement vs. flat e; negligible compute cost |
| T2 | Per-allele aggregated: separate e_ref_m and e_alt_m | FORMAT QA+QR (freebayes) or two BAM-pileup sums | Small incremental improvement; captures asymmetric miscalls |
| T3 | Per-read, exact product over Bernoullis | Pileup of individual reads | Most expensive; marginal gain when depth is uniformly sampled |

The plan below implements T1 as default, exposes T2 via a flag when the sidecar
data is available, and treats T3 as a future extension that is unlikely to be
worth the added BAM-parsing complexity.

## Detailed design

### 1. Data source: `FORMAT/QS` in the primary VCF

The joint-called VCF already fed to `allomix monitor --vcf ...` gains a new
`FORMAT/QS` field of `Number=R` (one entry per allele: REF, ALT). Produced
upstream by adding two lines to the variant-calling pipeline (Snakemake rule
or equivalent):

```bash
# After GATK HaplotypeCaller + GenotypeGVCFs produce joint_called.vcf.gz:
bcftools mpileup -a FORMAT/QS -f ref.fa -l sites.bed -Oz \
    -o pileup.vcf.gz --bam-list bams.txt
bcftools annotate -a pileup.vcf.gz -c FORMAT/QS \
    -Oz -o joint_called.qs.vcf.gz joint_called.vcf.gz
```

The annotated VCF then has, per sample per site:

```
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Sum of quality scores per allele">
chr1  1234  .  A  G  .  PASS  .  GT:AD:DP:QS  0/1:62,38:100:1960,1140
```

`QS[0]` is the summed BQ over reads reporting REF, `QS[1]` over ALT. Converting
to a per-marker effective error rate:

```
n = ad_ref + ad_alt
mean_phred = (qs_ref + qs_alt) / n
e = 10 ** (-mean_phred / 10)
```

This is the only data source. No BAM dependency, no sidecar TSV, no separate
`--bq-vcf` flag. The VCF either has QS (→ `--bq-aware` works) or it doesn't
(→ `--bq-aware` errors with a clear message pointing at the mpileup step).

### 2. Module layout

Simplest path: extend `src/allomix/genotype.py` to read `FORMAT/QS` inline
when parsing the admixture VCF, then add a small helper
`src/allomix/bq.py` with just the error-rate conversion. No separate module
for VCF parsing, no loader for a TSV sidecar, no `MarkerBQ` type.

Extend `parse_vcf` (genotype.py:62) to also extract QS:

```python
# src/allomix/genotype.py  (patch inside parse_vcf's record loop)
qs_arr = variant.format("QS")
qs_ref = qs_alt = None
if qs_arr is not None:
    qs_vals = qs_arr[sample_idx]
    if len(qs_vals) >= 2 and qs_vals[0] >= 0 and qs_vals[1] >= 0:
        qs_ref = int(qs_vals[0])
        qs_alt = int(qs_vals[1])
```

`MarkerData` gains `qs_ref: int | None = None` and `qs_alt: int | None = None`.
`classify_markers` then copies these onto `InformativeMarker.admix_qs_ref/alt`
from the admixture sample's `MarkerData` (no separate `admix_bq` dict
parameter needed).

The only new module is `src/allomix/bq.py`:

```python
# src/allomix/bq.py
"""Per-marker base quality handling.

Converts phred-summed per-allele base qualities (VCF FORMAT/QS produced by
bcftools mpileup -a FORMAT/QS) into per-marker effective error rates for
the chimerism likelihood.
"""

from __future__ import annotations


def effective_error_rate(
    ad_ref: int,
    ad_alt: int,
    qs_ref: int,
    qs_alt: int,
    floor_q: float = 2.0,
    cap_q: float = 45.0,
) -> float:
    """Return the per-marker effective error rate from phred-summed BQs.

    Uses the depth-weighted mean phred score over all reads:

        mean_phred = (qs_ref + qs_alt) / (ad_ref + ad_alt)
        e = 10^(-mean_phred/10)

    ``floor_q`` and ``cap_q`` bracket the mean phred before conversion;
    real BAMs sometimes carry Q0 soft-clipped bases and Q50+ pseudo-qualities
    that are not physically meaningful.

    Args:
        ad_ref: REF allele depth at this marker.
        ad_alt: ALT allele depth at this marker.
        qs_ref: Phred-summed BQ over REF-calling reads.
        qs_alt: Phred-summed BQ over ALT-calling reads.
        floor_q: Lower bound on mean phred (default 2).
        cap_q: Upper bound on mean phred (default 45).

    Returns:
        Effective per-read error rate in [10^(-cap_q/10), 10^(-floor_q/10)].
        Falls back to the caller-provided default if ad_ref+ad_alt == 0;
        callers should handle that case before calling.
    """
    n = ad_ref + ad_alt
    if n == 0:
        raise ValueError("n=0 at effective_error_rate; caller should skip this marker")
    mean_phred = (qs_ref + qs_alt) / n
    mean_phred = max(floor_q, min(cap_q, mean_phred))
    return 10.0 ** (-mean_phred / 10.0)
```

The floor/cap defaults are conservative (Q2–Q45). Q2 ≈ 0.63 error rate, anything
below that is already a soft-clip artefact. Q45 ≈ 3e-5, which is the limit where
BQs meaningfully discriminate between reads.

### 3. Attach BQ data to `InformativeMarker`

`MarkerData` and `InformativeMarker` both gain optional `qs_ref` / `qs_alt`
(and `admix_qs_ref` / `admix_qs_alt` on the informative variant). Defaults
`None` so fixtures and old paths keep working (same backward-compat pattern
used for `rho: float = float("inf")`).

```python
# src/allomix/genotype.py
@dataclass
class MarkerData:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int, int]
    ad_ref: int
    ad_alt: int
    dp: int
    gq: int | None = None
    filter: str = "PASS"
    qs_ref: int | None = None   # NEW: phred-summed BQ at REF-calling reads
    qs_alt: int | None = None   # NEW

@dataclass
class InformativeMarker:
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
    admix_qs_ref: int | None = None   # NEW
    admix_qs_alt: int | None = None   # NEW
```

`parse_vcf` reads `FORMAT/QS` when present (see section 2). `classify_markers`
copies `a.qs_ref`/`a.qs_alt` from the admixture `MarkerData` onto each new
`InformativeMarker`:

```python
# src/allomix/genotype.py — inside classify_markers informative branch
if any_informative:
    informative.append(
        InformativeMarker(
            ...,
            admix_qs_ref=a.qs_ref,
            admix_qs_alt=a.qs_alt,
        )
    )
```

No new parameters on `classify_markers` and no new dict to pass around. BQ data
rides through the existing parse → classify → estimate pipeline as ordinary
fields on the existing types. Host/donor QS is not used (we're estimating a
sequencing-error rate on the admixture sample only).

### 4. Likelihood change (T1: per-marker error rate)

In `src/allomix/chimerism.py:133-182`, `log_likelihood_marker_bb` takes a scalar
`error_rate`. Extend it to accept a per-call override:

```python
def log_likelihood_marker_bb(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
) -> float:
    ...
```

Leave this function untouched. Instead, have the per-marker callers resolve the
effective error rate **before** calling it. Concretely, change
`total_log_likelihood_bb` and `total_log_likelihood_multi_bb`:

```python
# src/allomix/chimerism.py
from allomix.bq import effective_error_rate  # NEW import

def total_log_likelihood_bb(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    use_bq: bool = False,                               # NEW
) -> float:
    ll = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)

        # Per-marker effective error rate when BQ data is attached.
        e = error_rate
        if use_bq and m.admix_qs_ref is not None and m.admix_qs_alt is not None:
            if m.admix_ad_ref + m.admix_ad_alt > 0:
                e = effective_error_rate(
                    m.admix_ad_ref, m.admix_ad_alt,
                    m.admix_qs_ref, m.admix_qs_alt,
                )

        ll += log_likelihood_marker_bb(m.admix_ad_ref, m.admix_ad_alt, w, e, rho)
    return ll
```

The same pattern applies to `total_log_likelihood_multi_bb`.

Then `estimate_single_donor_bb` and `estimate_multi_donor` gain a `use_bq: bool`
parameter that they plumb through to every call of the total-log-likelihood
function — including inside the profile-likelihood scans, which use the same
total functions via closure.

```python
# src/allomix/chimerism.py — example patch inside estimate_single_donor_bb
def estimate_single_donor_bb(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    use_bq: bool = False,                               # NEW
) -> ChimerismResult:
    ...
    # Step 1: Grid search — pass use_bq into every total_log_likelihood_bb call
    for f in grid:
        opt_rho = minimize_scalar(
            lambda log_r: (
                -total_log_likelihood_bb(
                    markers, f, error_rate, math.exp(log_r), marker_biases,
                    use_bq=use_bq,                      # NEW
                )
            ),
            bounds=(math.log(1.0), math.log(10000.0)),
            method="bounded",
        )
        ...
    # Nelder-Mead, profile LL, CI code all receive use_bq through closure:
    def neg_ll_joint(x):
        f_val, log_rho_val = x
        ...
        return -total_log_likelihood_bb(
            markers, f_val, error_rate, rho_val, marker_biases, use_bq=use_bq,
        )
    ...
    def profile_ll_f(f_val: float) -> float:
        opt_rho = minimize_scalar(
            lambda log_r: (
                -total_log_likelihood_bb(
                    markers, f_val, error_rate, math.exp(log_r), marker_biases,
                    use_bq=use_bq,
                )
            ),
            bounds=(math.log(1.0), math.log(50000.0)),
            method="bounded",
        )
        return -float(opt_rho.fun)
    ...
```

The `ChimerismResult.error_rate` field now becomes ambiguous: is it the global
default, or is it per-marker? Recommend keeping `error_rate` as "the default /
fallback rate", and adding a small summary field:

```python
@dataclass
class ChimerismResult:
    ...
    error_rate: float
    rho: float = float("inf")
    bq_aware: bool = False                  # NEW: True if use_bq was on
    mean_effective_error: float | None = None  # NEW: diagnostic mean e across markers
```

`mean_effective_error` is for the TSV/JSON diagnostic only; it summarises the
spread of per-marker effective error rates. Set it at the end of
`estimate_single_donor_bb` from the per-marker values that were actually used at
the MLE:

```python
if use_bq:
    effective_errors = []
    for m in markers:
        if m.admix_qs_ref is not None and m.admix_qs_alt is not None:
            if m.admix_ad_ref + m.admix_ad_alt > 0:
                effective_errors.append(effective_error_rate(
                    m.admix_ad_ref, m.admix_ad_alt,
                    m.admix_qs_ref, m.admix_qs_alt,
                ))
    mean_eff = (
        float(np.mean(effective_errors)) if effective_errors else None
    )
else:
    mean_eff = None

return ChimerismResult(
    ...,
    error_rate=error_rate,
    rho=rho_mle,
    bq_aware=use_bq,
    mean_effective_error=mean_eff,
)
```

### 5. GoF variance in `qc.py`

`_compute_gof_pval` already takes `error_rate` as a scalar. Under T1 the per-marker
variance floor should use each marker's effective error rate, not the global one.
Two options:

#### Option A: Pass per-marker effective errors as an optional list

```python
# src/allomix/qc.py
def _compute_gof_pval(
    per_marker: list[MarkerResult],
    rho: float = float("inf"),
    n_fitted_params: int = 2,
    error_rate: float = 0.0,
    per_marker_error_rate: list[float] | None = None,   # NEW
) -> float | None:
    ...
    for idx, m in enumerate(included):
        n = m.dp
        if n <= 0:
            continue
        e_m = (
            per_marker_error_rate[idx]
            if per_marker_error_rate is not None
            else error_rate
        )
        if e_m > 0:
            ev_raw = _error_adjusted_p_alt(m.expected_vaf, e_m)
        else:
            ev_raw = m.expected_vaf
        ev = max(1e-6, min(1.0 - 1e-6, ev_raw))
        ...
```

And `assess_quality` (`qc.py:141`) computes the list from the result's stored
BQ diagnostics. **But**: `MarkerResult` does not today carry BQ data, and we do
not want to add a new `effective_error_rate` field to every `MarkerResult` just
for T1. Cheaper to attach it there once we need it:

```python
# src/allomix/chimerism.py
@dataclass
class MarkerResult:
    ...
    included: bool
    effective_error_rate: float | None = None   # NEW
```

In `_compute_per_marker_results` (chimerism.py:278-339), set this when BQ data is
available. Then `_compute_gof_pval` can ignore the new qc-param entirely and read
`m.effective_error_rate` directly:

```python
for m in included:
    ...
    e_for_floor = m.effective_error_rate if m.effective_error_rate is not None else error_rate
    if e_for_floor > 0:
        ev_raw = _error_adjusted_p_alt(m.expected_vaf, e_for_floor)
    else:
        ev_raw = m.expected_vaf
    ...
```

This is cleaner and keeps the signature of `_compute_gof_pval` unchanged. Prefer
this path.

### 6. CLI

One flag on `monitor` and `timeline`:

```python
# src/allomix/cli.py — _add_common_args
parser.add_argument(
    "--bq-aware",
    action="store_true",
    help="Use per-marker base qualities (FORMAT/QS) in the likelihood. "
         "Requires the input VCF to carry FORMAT/QS (produced by "
         "bcftools mpileup -a FORMAT/QS + bcftools annotate).",
)
```

`_run_single_sample` gains one new kwarg (`use_bq`) and passes it through. No
QS loading plumbing in `cmd_monitor` — `parse_vcf` already reads QS inline
when the field is present:

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
    use_bq: bool = False,           # NEW
) -> tuple:
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0)
    genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
    genotypes.sample_name = admix_sample
    ...
    if len(donors) == 1:
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
            use_bq=use_bq,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            marker_biases=marker_biases,
            use_bq=use_bq,
        )
    qc = assess_quality(result, genotypes)
    return result, qc, genotypes
```

When `--bq-aware` is set but no `InformativeMarker` carries QS, raise a
SystemExit with a message pointing at the upstream mpileup step, e.g.:

```python
# in estimate_single_donor_bb / estimate_multi_donor, or in cmd_monitor
if use_bq and not any(
    m.admix_qs_ref is not None and m.admix_qs_alt is not None
    for m in markers
):
    raise SystemExit(
        "--bq-aware requested but no FORMAT/QS found in VCF. "
        "Re-run the variant-calling pipeline with "
        "`bcftools mpileup -a FORMAT/QS` + `bcftools annotate -c FORMAT/QS`, "
        "or drop --bq-aware to use the flat --error-rate."
    )
```

Per-marker fallback behaviour (when most markers have QS but a few don't) still
uses the global `error_rate` at those markers; the guard above is only for the
"QS is entirely absent" case.

### 7. Report columns

Add two columns next to the existing TSV summary and one field in the JSON:

```python
# src/allomix/report.py
summary_header = (
    "sample\tdonor_pct\tci_lo\tci_hi\tn_informative\tn_used\t"
    "mean_depth\tgof_pval\tbq_aware\tmean_eff_error\tqc_pass"
)
...
bq_aware_str = "Y" if getattr(result, "bq_aware", False) else "N"
mean_eff = getattr(result, "mean_effective_error", None)
mean_eff_str = f"{mean_eff:.4f}" if mean_eff is not None else "NA"
summary_line = (
    f"{sample_name}\t..."
    f"{bq_aware_str}\t"
    f"{mean_eff_str}\t"
    f"{qc_pass_str}"
)
```

JSON mirrors the same two fields at the top level of the result object. Avoid
renaming existing columns — downstream `scripts/run_xls_batch.py` and paper
figure scripts parse the TSV.

### 8. Simulator support (for in-silico validation)

To validate the new likelihood we must also simulate per-read BQ. Extend
`simulate.py`:

```python
# src/allomix/simulate.py
def sample_allele_counts_with_bq(
    vaf: float,
    depth: int,
    rng: random.Random,
    bq_distribution: list[int] | None = None,  # empirical phred histogram
    mean_phred: float = 30.0,
    sd_phred: float = 5.0,
) -> tuple[int, int, int, int]:
    """Like sample_allele_counts, but also emits (qs_ref, qs_alt).

    For each simulated read:
        1. Draw a base quality Q_i (from bq_distribution if given, else
           from a clipped Normal(mean_phred, sd_phred)).
        2. Compute e_i = 10^(-Q_i/10).
        3. Apply the 4-state error model with per-read e_i.
        4. Emit REF or ALT + add Q_i to qs_ref or qs_alt.

    Returns (ref_count, alt_count, qs_ref, qs_alt).
    """
    if depth <= 0:
        return (0, 0, 0, 0)
    ref_count = alt_count = qs_ref = qs_alt = 0
    for _ in range(depth):
        if bq_distribution:
            q = rng.choice(bq_distribution)
        else:
            q = max(2, min(45, int(round(rng.gauss(mean_phred, sd_phred)))))
        e = 10.0 ** (-q / 10.0)
        p_alt = vaf * (1.0 - e) + (1.0 - vaf) * e / 3.0
        p_ref_or_alt = p_alt + ((1.0 - vaf) * (1.0 - e) + vaf * e / 3.0)
        p_alt_cond = p_alt / p_ref_or_alt
        if rng.random() < p_alt_cond:
            alt_count += 1
            qs_alt += q
        else:
            ref_count += 1
            qs_ref += q
    return (ref_count, alt_count, qs_ref, qs_alt)
```

Note: this is O(depth) per marker (vs. the current O(1) `binomialvariate`), so at
1000x × 100 markers × 100 replicates it will take a few minutes instead of
seconds. Acceptable for validation runs; avoid regressing the fast-path.

`blend_vcfs` gains a `bq_mean`/`bq_sd`/`bq_distribution` parameter set. When
any is non-None, the output VCF should emit `FORMAT/QS`. Otherwise fall through
the existing path unchanged, producing no QS field.

### 9. BAM extraction script (mode C producer)

Sketch:

```python
# scripts/extract_bq_from_bam.py
"""Emit a per-marker BQ table (chrom, pos, ref, alt, sample, ad_ref, ad_alt,
qs_ref, qs_alt) by pileup over a BAM at the sites in a panel BED/VCF.

Usage:
    python scripts/extract_bq_from_bam.py \\
        --bam sample.bam --bed panel.bed --sample SAMPLE_NAME \\
        -o sample_bq.tsv
"""
from __future__ import annotations
import argparse
import pysam

# ... read sites from BED or VCF, pileup with min_base_quality=0,
#     ignore secondary/supplementary, sum BQ per allele, write TSV ...
```

Not shipped in the `allomix` package; lives in `scripts/`. Add a short
CLAUDE.md-compliant docstring and leave pysam as a scripts-only dependency.

## Validation plan

### 10. Unit tests

`tests/test_bq.py` (new):

- `effective_error_rate` edge cases: n=0 raises; all Q20 → e≈0.01; caps at
  floor_q/cap_q when either side is extreme.
- `parse_qs_from_vcf` on a small synthetic VCF with known QS values round-trips
  cleanly.
- `load_bq_table` round-trip against a written-out TSV.

Extend `tests/test_chimerism.py`:

- With flat mean BQ equal to `-10*log10(error_rate)`, the BQ-aware likelihood
  should equal the flat likelihood to within numerical tolerance.
- With heterogeneous BQs (half Q40, half Q10), the MLE should still recover the
  true donor fraction but the CI should differ from the flat case.

Extend `tests/test_cli.py`:

- `--bq-aware` on a VCF without QS should error with a clear message, not a
  silent fallback (catches the common mistake of forgetting `--bq-vcf`).
- `--bq-aware --bq-table path/to.tsv` round-trips and produces `bq_aware=Y` in
  the TSV output.

### 11. Synthetic-data validation

`paper/scripts/run_bq_validation.py` (new, modelled on `run_depth_validation.py`):

- Generate chimeric VCFs at f ∈ {0.01, 0.05, 0.1, 0.25, 0.5} across two BQ
  regimes: (a) uniform Q30 and (b) mixed Q10/Q30 (60/40 split).
- Run monitor with and without `--bq-aware`; record MAE, CI coverage, CI width.
- Expected: under regime (a) BQ-aware and flat are indistinguishable
  (mean_effective_error ≈ `--error-rate`). Under regime (b) BQ-aware should
  either (i) narrow CIs with equal coverage or (ii) widen them and improve
  coverage, depending on which direction the unmodelled noise was biasing the
  flat fit. Either outcome is informative.

### 12. Real-data smoke test

Rerun the April 24 validation batch once the pipeline has been updated to emit
FORMAT/QS. Steps:

1. Update the upstream Snakefile to run `bcftools mpileup -a FORMAT/QS` +
   `bcftools annotate` against the joint-called VCF, producing a QS-annotated
   replacement at `output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz`.
2. `allomix monitor --vcf output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz
   --host-sample H --donor-sample D --sample A --bq-aware -o out.tsv`
3. Compare the resulting donor fractions and gof_pvals to the
   `validation_run_new_bias2/` baseline.

Record which samples move most and how much. Where the flat-error fit was
already clean, we expect small movement. Where it was borderline (wide CIs,
low gof_pval), BQ-aware should push in one direction consistently.

## Scope boundary

Inside scope:

- T1 per-marker effective error from FORMAT/QS or a BQ table.
- CLI flags, report columns, `MarkerResult.effective_error_rate`.
- Simulator extension for validation.
- Unit + synthetic validation.

Out of scope (future work, flagged but not implemented):

- T2 per-allele asymmetric error (`e_ref_m ≠ e_alt_m`). Cheap to add once T1
  works; punt until we see evidence it matters.
- T3 per-read product likelihood. Requires BAM pileup in the core package,
  substantially more compute, and Poisson-binomial machinery. Only worth
  returning to if T1 validation shows that the per-marker aggregation loses
  real signal.
- Integration with `--bias-table`: when both are on, bias correction and
  per-marker error rates interact through the same `expected_weight` →
  `p_alt` path. The current patch composes them correctly (bias shifts w,
  BQ shifts e), but worth a dedicated validation.
- Paper update: methods section gains a sub-subsection on BQ weighting;
  a figure comparing flat vs BQ-aware MAE at matched CI coverage. Defer until
  in-silico validation lands.

## Edge cases and risks

- **Missing QS on some markers**: If QS is present on most markers but absent on
  a few, the per-marker code path should fall back to the global `error_rate`
  for those markers without silent failure. `effective_error_rate` raises on
  n=0; callers already skip those. For "n>0 but QS is None" treat as
  "use global rate" — this is the point of making `admix_qs_ref/alt`
  `int | None`.

- **QS at zero depth**: If the pileup reports ad_ref=0, ad_alt=0, don't divide
  by zero — skip BQ weighting at that marker. Already handled: we only call
  `effective_error_rate` when `ad_ref + ad_alt > 0`.

- **Fitted rho interaction**: rho absorbs extra-binomial variance. With a
  sharper error model the residuals shrink at high-quality markers, which can
  drive rho up (less overdispersion "left over") or leave it roughly unchanged
  (depending on whether BQ variation is the dominant residual noise source).
  Not a bug — just worth flagging when comparing ChimerismResult.rho values
  across runs with/without `--bq-aware`.

- **Cache warming and per-call overhead**: `effective_error_rate` is called
  inside every `total_log_likelihood_bb` evaluation, which runs thousands of
  times during grid search + profile CIs. Precompute per-marker effective
  errors once before the optimiser starts:

  ```python
  def estimate_single_donor_bb(...):
      if use_bq:
          precomputed_errors = []
          for m in markers:
              if (m.admix_qs_ref is not None and m.admix_qs_alt is not None
                      and m.admix_ad_ref + m.admix_ad_alt > 0):
                  precomputed_errors.append(effective_error_rate(
                      m.admix_ad_ref, m.admix_ad_alt,
                      m.admix_qs_ref, m.admix_qs_alt,
                  ))
              else:
                  precomputed_errors.append(error_rate)
      else:
          precomputed_errors = [error_rate] * len(markers)
  ```

  Then pass `precomputed_errors` through the total-log-likelihood functions as
  a pre-indexed list (same length as `markers`). This avoids re-doing the
  effective-error math ~1e5 times per sample.

  This suggests a cleaner function signature: instead of `use_bq: bool` on
  every call, pass `per_marker_error_rate: list[float] | None = None`:

  ```python
  def total_log_likelihood_bb(
      markers: list[InformativeMarker],
      f_donor: float,
      error_rate: float = 0.01,
      rho: float = 100.0,
      marker_biases: dict[tuple[str, int, str, str], float] | None = None,
      per_marker_error_rate: list[float] | None = None,   # NEW
  ) -> float:
      ll = 0.0
      for i, m in enumerate(markers):
          ...
          e = per_marker_error_rate[i] if per_marker_error_rate is not None else error_rate
          ll += log_likelihood_marker_bb(m.admix_ad_ref, m.admix_ad_alt, w, e, rho)
      return ll
  ```

  This is the recommended signature. Precompute the list once in
  `estimate_single_donor_bb` / `estimate_multi_donor`; pass it to every
  inner call.

- **VCF that claims QS but has all zeros**: bcftools emits `0,0` when no reads
  pass the filters. Treat as "no BQ data" at that marker. The
  `ad_ref + ad_alt > 0` guard handles this as long as DP is consistent with
  AD; add a safety check `if qs_ref + qs_alt == 0: e = error_rate`.

- **Mean-phred vs mean-error approximation**: `effective_error_rate` converts
  `QS/n` (mean phred) to a single error rate via `10^(-mean_phred/10)`. The
  quantity we actually want in the likelihood aggregation is `mean(10^(-Q_i/10))`
  (arithmetic mean of per-read error rates). By Jensen's inequality for the
  convex function `10^(-x/10)`:

      10^(-mean(Q)/10)  <=  mean(10^(-Q_i/10))

  so the mean-phred conversion systematically underestimates the effective
  error rate when BQs are heterogeneous. For a mix of half Q30, half Q10
  reads, the true mean e is ~0.05 but the mean-phred estimate gives 10^(-2) =
  0.01. That's a 5x underestimate at the extreme end. Mitigations:

  1. Accept the approximation. Most reads within a marker have similar BQs
     (recalibrated Illumina output is fairly tight around Q30±3), so the gap
     is typically under 20%. Document the approximation and move on.
  2. Ask bcftools mpileup to emit a sum-of-error-rates alongside QS (not a
     standard tag; would need a custom annotation). Not worth the engineering
     cost for T1.
  3. Fall back to the per-read product (T3) when accuracy matters — that's
     exact by construction.

  Recommend (1) for T1. The approximation has a defined bias direction
  (underestimates e, so overestimates confidence), which matters for CIs but
  less so for point estimates. Flag in methods. If the in-silico validation
  shows this eating into CI coverage, revisit.

## File-by-file checklist

- [ ] Upstream pipeline (Snakefile or equivalent): add a `bcftools mpileup
  -a FORMAT/QS` + `bcftools annotate -c FORMAT/QS` step that augments the
  joint-called VCF with FORMAT/QS before it lands in `output/joint_called/`.
- [ ] `src/allomix/bq.py` (NEW): single helper, `effective_error_rate(ad_ref,
  ad_alt, qs_ref, qs_alt, floor_q=2.0, cap_q=45.0) -> float`.
- [ ] `src/allomix/genotype.py`:
    - Add `qs_ref: int | None = None`, `qs_alt: int | None = None` to
      `MarkerData`.
    - In `parse_vcf`, read `FORMAT/QS` when present and populate those fields.
    - Add `admix_qs_ref: int | None = None`, `admix_qs_alt: int | None = None`
      to `InformativeMarker`.
    - In `classify_markers` informative branch, copy `a.qs_ref`/`a.qs_alt`
      onto the new `InformativeMarker`.
- [ ] `src/allomix/chimerism.py`:
    - Add `effective_error_rate: float | None = None` to `MarkerResult`.
    - Add `bq_aware: bool = False` and `mean_effective_error: float | None = None`
      to `ChimerismResult` and `MultiDonorResult`.
    - Change `total_log_likelihood_bb` / `total_log_likelihood_multi_bb` to take
      an optional `per_marker_error_rate` list; use it in place of scalar `e`.
    - In `estimate_single_donor_bb` / `estimate_multi_donor`: accept `use_bq`,
      precompute the per-marker error-rate list once (using `effective_error_rate`
      where QS is present, falling back to scalar `error_rate` where it isn't),
      and pass to all total-LL calls (including those inside `profile_ll_f`,
      `neg_ll_joint`, `_profile_likelihood_cis_multi`).
    - Raise a clear SystemExit when `use_bq=True` but no marker carries QS.
    - Populate `MarkerResult.effective_error_rate` in
      `_compute_per_marker_results` / `_per_marker_results_multi` from the
      precomputed list.
    - Set `bq_aware` and `mean_effective_error` on the returned result.
- [ ] `src/allomix/qc.py`: in `_compute_gof_pval`, use
  `m.effective_error_rate` if set, else `error_rate`. No new signature change.
- [ ] `src/allomix/report.py`: add `bq_aware` and `mean_eff_error` columns
  to TSV summary; add equivalent fields to JSON output.
- [ ] `src/allomix/cli.py`: add `--bq-aware` to `_add_common_args`; thread
  `use_bq` into `_run_single_sample`. No new VCF-loading plumbing.
- [ ] `src/allomix/simulate.py`: add `sample_allele_counts_with_bq` and an
  opt-in `bq_mean`/`bq_sd` arm on `blend_vcfs` that emits FORMAT/QS.
- [ ] `tests/test_genotype.py`: add tests for QS parsing and round-trip.
- [ ] `tests/test_bq.py` (NEW): unit tests for `effective_error_rate`.
- [ ] `tests/test_chimerism.py`: add tests exercising `use_bq` paths.
- [ ] `tests/test_cli.py` / `tests/test_integration.py`: CLI flag plumbing
  tests (including the "QS missing" error message).
- [ ] `paper/scripts/run_bq_validation.py` (NEW): in-silico validation script.
- [ ] `claude/allomix_overall_plan.md`: link Step 12 to this plan, mark
  status.

## Verification plan (run order)

1. Unit tests for `bq.py` (fastest).
2. `pytest tests/test_chimerism.py tests/test_qc.py -x -q` — ensure no
   regressions under flat-e path (default `use_bq=False`).
3. Run `tests/test_integration.py` with a hand-built synthetic QS VCF to
   verify CLI plumbing end-to-end.
4. Run `paper/scripts/run_bq_validation.py` at modest scale (5 replicates per
   cell) for a first-pass sanity check.
5. Scale up validation and regenerate any paper figures that compare error
   models (likely only if we decide to add a figure).
6. Real-data: rerun the April 24 batch with `--bq-aware` once a sidecar
   bcftools-mpileup VCF is available.
