# Step 7: Multi-Donor Chimerism Estimation — Implementation Plan

## Status

**Code: ✅ COMPLETE** (2026-04-08) — 261 tests pass (48 new multi-donor + 213 existing, zero regressions).

**Paper: 📝 TODO** — See `claude/allomix_overall_plan.md` Step 13 for detailed list of paper sections to update (methods formula, results validation, figure, discussion/intro/abstract/README updates).

## Overview

Extend allomix from single-donor (host + 1 donor) to multi-donor (host + 2 donors) chimerism estimation. This is needed clinically for dual cord blood transplants and cases where a patient receives cells from two donors. The current codebase is already architecturally prepared — `InformativeMarker.donor_gts` stores a list, the CLI accepts `--donor` multiple times, and `classify_markers()` joins across multiple donor VCFs.

The mathematical extension is clean: the 1D likelihood over donor fraction *f* becomes a 2D likelihood over (*f₁*, *f₂*) with the constraint *f₁ + f₂ ≤ 1*.

## Clinical Scenario for Validation

Three brothers (siblings with shared parents): Brother A is the host/recipient, Brothers B and C are donors. This is the hardest realistic case because:

1. **Sibling IBD reduces informative markers** — full siblings share genotype at ~25% of loci (IBD=2), reducing the markers where we can distinguish contributions
2. **Correlated donor genotypes** — both donors are siblings of each other AND of the host, so there are loci where all three are identical (uninformative) and loci where both donors differ from host in the same way (informative for total donor but not for distinguishing donors)
3. **Three-way relatedness** — unlike generating donor-vs-host pairs independently, the three siblings must be generated from shared parents to preserve realistic genotype correlations

This scenario will produce a paper figure demonstrating multi-donor estimation accuracy even under the challenging case of related donors.

---

## Part 1: Genotype Generation — Three Siblings from Shared Parents

### Problem with current approach

`generate_related_genotypes()` in `simulate.py` generates one donor conditional on a host. For 3 siblings, we can't generate donor1|host and donor2|host independently — that would ignore the correlation between donors (they share the same parents).

### Solution: Mendelian segregation from parents

Generate two parents from population allele frequencies, then derive each sibling's genotype by independent Mendelian segregation from those parents.

```python
# New function in simulate.py

def generate_sibling_trio_genotypes(
    n_markers: int,
    rng: random.Random,
    maf_range: tuple[float, float] = (0.2, 0.5),
) -> list[dict]:
    """Generate genotypes for 3 siblings (host + 2 donors) from shared parents.

    For each marker:
    1. Draw population ALT allele frequency
    2. Draw two parent genotypes from Hardy-Weinberg
    3. Derive each sibling independently by Mendelian segregation

    This preserves the correct 3-way sibling correlation structure:
    - Each pair has IBD distribution (0.25, 0.5, 0.25)
    - All three may be IBD=2 (identical) at some loci
    - Donor1 and donor2 are correlated even conditional on host

    Args:
        n_markers: Number of biallelic markers to generate.
        rng: Random instance for reproducibility.
        maf_range: (min, max) minor allele frequency range.

    Returns:
        List of dicts with keys: chrom, pos, ref, alt, host_gt, donor1_gt,
        donor2_gt, p_alt, parent1_gt, parent2_gt.
    """
    markers = []
    for i in range(n_markers):
        p_alt = rng.uniform(*maf_range)

        # Draw two parents from HWE
        parent1 = _draw_genotype(p_alt, rng)
        parent2 = _draw_genotype(p_alt, rng)

        # Each sibling gets one allele from each parent (independent segregation)
        host_gt = _mendelian_child(parent1, parent2, rng)
        donor1_gt = _mendelian_child(parent1, parent2, rng)
        donor2_gt = _mendelian_child(parent1, parent2, rng)

        markers.append({
            "chrom": f"chr{(i % 22) + 1}",
            "pos": 1_000_000 + i * 100_000,
            "ref": "A",
            "alt": "G",
            "host_gt": host_gt,
            "donor1_gt": donor1_gt,
            "donor2_gt": donor2_gt,
            "p_alt": p_alt,
            "parent1_gt": parent1,
            "parent2_gt": parent2,
            # Informativity flags
            "informative_d1": alt_dose(host_gt) != alt_dose(donor1_gt),
            "informative_d2": alt_dose(host_gt) != alt_dose(donor2_gt),
            "informative_any": (
                alt_dose(host_gt) != alt_dose(donor1_gt)
                or alt_dose(host_gt) != alt_dose(donor2_gt)
            ),
            # Can we distinguish donor1 from donor2?
            "donors_distinguishable": alt_dose(donor1_gt) != alt_dose(donor2_gt),
        })

    return markers


def _mendelian_child(
    parent1: tuple[int, int],
    parent2: tuple[int, int],
    rng: random.Random,
) -> tuple[int, int]:
    """Draw a child genotype by Mendelian segregation from two parents.

    Each parent transmits one allele (chosen uniformly at random).
    """
    a1 = parent1[rng.randint(0, 1)]
    a2 = parent2[rng.randint(0, 1)]
    return (min(a1, a2), max(a1, a2))
```

### Expected informativity for 3 siblings

With MAF range 0.2–0.5 and 100 markers:
- ~80 markers informative for at least one donor (host ≠ donor_i for some i)
- ~60 markers informative for donor1 specifically
- ~60 markers informative for donor2 specifically
- ~40 markers where both donors differ from host AND differ from each other (best markers for resolving individual contributions)

The key question is whether the loci where both donors differ from host *in the same direction* (e.g., host 0/0, donor1 0/1, donor2 0/1) can still help. Answer: yes — they contribute to total donor fraction estimation, even though they can't distinguish *which* donor contributed.

---

## Part 2: Multi-Donor Simulation (VCF Blending)

### New function: `blend_vcfs_multi`

Extend `blend_vcfs` to accept multiple donors and fractions.

```python
# In simulate.py

def expected_vaf_multi(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
) -> float:
    """Expected ALT VAF in a multi-donor chimeric mixture.

    VAF = ((1 - f1 - f2) * host_dose + f1 * donor1_dose + f2 * donor2_dose) / 2

    Args:
        host_gt: Host diploid genotype.
        donor_gts: List of donor diploid genotypes.
        donor_fractions: List of donor fractions (must sum to <= 1.0).

    Returns:
        Expected ALT allele frequency.
    """
    f_host = 1.0 - sum(donor_fractions)
    vaf = f_host * alt_dose(host_gt)
    for dgt, f in zip(donor_gts, donor_fractions):
        vaf += f * alt_dose(dgt)
    return vaf / 2.0


def blend_vcfs_multi(
    host_path: str | Path,
    donor_paths: list[str | Path],
    donor_fractions: list[float],
    target_depth: int | None = None,
    sample_name: str | None = None,
    seed: int | None = None,
    error_rate: float = 0.01,
    depth_cv: float = 0.0,
    realistic_biases: bool = False,
) -> BlendResult:
    """Blend host + multiple donor VCFs to create a synthetic chimeric VCF.

    Args:
        host_path: Path to host genotype VCF.
        donor_paths: List of paths to donor genotype VCFs.
        donor_fractions: List of donor DNA fractions (must sum to <= 1.0).
        target_depth: Fixed depth for all markers.
        sample_name: Sample name for output.
        seed: Random seed.
        error_rate: Per-read sequencing error rate.
        depth_cv: Depth coefficient of variation.
        realistic_biases: Use heavy-tailed bias distribution.

    Returns:
        BlendResult with synthetic chimeric VCF data.
    """
    if sum(donor_fractions) > 1.0 + 1e-9:
        raise ValueError(
            f"donor_fractions sum to {sum(donor_fractions):.4f}, must be <= 1.0"
        )
    if len(donor_paths) != len(donor_fractions):
        raise ValueError("donor_paths and donor_fractions must have same length")

    rng = random.Random(seed)
    host_header, host_records = parse_vcf(host_path)

    # Parse all donor VCFs and index by locus
    donor_records_list = []
    for dp in donor_paths:
        _, records = parse_vcf(dp)
        donor_records_list.append({rec.locus: rec for rec in records})

    # ... (join on shared loci, compute expected_vaf_multi, sample counts)
    # Implementation follows the same pattern as blend_vcfs but with
    # expected_vaf_multi() instead of expected_vaf()
```

### Alternative: generate from marker dicts directly

For the sibling trio test data, we can skip VCF round-tripping and generate directly from the marker dicts produced by `generate_sibling_trio_genotypes()`. This is simpler and avoids needing to write 3 separate genotype VCFs just to blend them.

```python
def blend_from_genotype_dicts(
    markers: list[dict],
    donor_fractions: list[float],
    target_depth: int = 1000,
    seed: int | None = None,
    error_rate: float = 0.01,
    depth_cv: float = 0.0,
) -> BlendResult:
    """Create synthetic chimeric VCF directly from genotype dicts.

    Designed for use with generate_sibling_trio_genotypes() output.

    Args:
        markers: List of marker dicts with host_gt, donor1_gt, donor2_gt.
        donor_fractions: [f_donor1, f_donor2].
        target_depth: Mean sequencing depth.
        seed: Random seed.
        error_rate: Sequencing error rate.
        depth_cv: Depth CV across markers.

    Returns:
        BlendResult with synthetic chimeric VCF data.
    """
    rng = random.Random(seed)
    n = len(markers)

    if depth_cv > 0:
        depths = sample_marker_depths(n, target_depth, depth_cv, rng)
    else:
        depths = [target_depth] * n

    header = _build_vcf_header("simulated")
    out_records = []
    n_informative = 0

    for i, m in enumerate(markers):
        host_gt = m["host_gt"]
        donor_gts = [m["donor1_gt"], m["donor2_gt"]]

        vaf = expected_vaf_multi(host_gt, donor_gts, donor_fractions)
        ref_count, alt_count = sample_allele_counts(vaf, depths[i], rng, error_rate)

        if m["informative_any"]:
            n_informative += 1

        # Build VCF record line ...
        out_records.append(_format_vcf_record(m, ref_count, alt_count))

    return BlendResult(
        header=header,
        records=out_records,
        num_markers=n,
        num_informative=n_informative,
    )
```

---

## Part 3: Core Algorithm — 2D MLE Estimation

### 3a. Extended likelihood function

The key mathematical change is extending the expected weight formula from 1 donor to 2:

**Single donor:**
$$w_i(f) = (1 - f) \cdot \frac{g_{h,i}}{2} + f \cdot \frac{g_{d,i}}{2}$$

**Two donors:**
$$w_i(f_1, f_2) = (1 - f_1 - f_2) \cdot \frac{g_{h,i}}{2} + f_1 \cdot \frac{g_{d1,i}}{2} + f_2 \cdot \frac{g_{d2,i}}{2}$$

The per-marker log-likelihood formula is unchanged — only the expected weight changes.

```python
# In chimerism.py

def expected_weight_multi(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
    bias: float = 0.0,
) -> float:
    """Expected reference allele weight for multi-donor chimerism.

    w = (1 - f1 - f2) * host_ref_dose/2 + f1 * d1_ref_dose/2 + f2 * d2_ref_dose/2

    Args:
        host_gt: Host diploid genotype.
        donor_gts: List of donor diploid genotypes.
        donor_fractions: List of donor fractions (sum <= 1.0).
        bias: Per-marker amplification bias.

    Returns:
        Expected reference allele weight (0.0 to 1.0).
    """
    host_ref_dose = 2 - (host_gt[0] + host_gt[1])
    f_host = 1.0 - sum(donor_fractions)
    w = f_host * host_ref_dose / 2.0
    for dgt, f in zip(donor_gts, donor_fractions):
        d_ref_dose = 2 - (dgt[0] + dgt[1])
        w += f * d_ref_dose / 2.0
    if bias != 0.0:
        w = max(1e-6, min(1.0 - 1e-6, w - bias))
    return w


def total_log_likelihood_multi(
    markers: list[InformativeMarker],
    donor_fractions: list[float],
    error_rate: float = 0.01,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """Total log-likelihood for multi-donor model.

    A marker contributes to the likelihood if it is informative for ANY
    donor (i.e., host differs from at least one donor). The expected
    weight uses ALL donor genotypes simultaneously.

    Args:
        markers: Informative markers (informative for at least one donor).
        donor_fractions: [f_donor1, f_donor2].
        error_rate: Sequencing error rate.
        marker_biases: Optional per-marker bias dict.

    Returns:
        Total log-likelihood.
    """
    ll = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)
        ll += log_likelihood_marker(m.admix_ad_ref, m.admix_ad_alt, w, error_rate)
    return ll
```

### 3b. Optimization: 2D grid search + Nelder-Mead

The constraint is *f₁ + f₂ ≤ 1* with both *f₁, f₂ ≥ 0*. This defines a triangular feasible region (a simplex).

**Strategy:**
1. **Triangular grid search**: Evaluate likelihood on a grid within the triangle *f₁ + f₂ ≤ 1*. At resolution *k*, this gives *k(k+1)/2* points. With *k=100*, that's 5,050 evaluations — fast enough.
2. **Nelder-Mead refinement**: Use `scipy.optimize.minimize` with Nelder-Mead from the grid maximum. Include the simplex constraint via a penalty or by clamping.

```python
def estimate_multi_donor(
    markers: list[InformativeMarker],
    n_donors: int = 2,
    error_rate: float = 0.01,
    grid_steps: int = 101,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> MultiDonorResult:
    """Estimate multi-donor chimerism fractions via maximum likelihood.

    Algorithm:
        1. Triangular grid search over (f1, f2) with f1 + f2 <= 1
        2. Nelder-Mead refinement from grid maximum
        3. Profile likelihood CI per donor
        4. Per-marker residuals and outlier flagging

    Args:
        markers: Informative markers (for at least one donor).
        n_donors: Number of donors (currently 2).
        error_rate: Sequencing error rate.
        grid_steps: Grid resolution per dimension.
        marker_biases: Optional per-marker bias dict.

    Returns:
        MultiDonorResult with per-donor fractions and CIs.
    """
    n_informative = len(markers)
    if n_informative == 0:
        return _empty_multi_result(n_donors, error_rate)

    # Step 1: Triangular grid search
    best_ll = -math.inf
    best_f = [0.0] * n_donors
    step = 1.0 / (grid_steps - 1)

    for i in range(grid_steps):
        f1 = i * step
        max_f2 = 1.0 - f1
        n_f2_steps = max(1, int(max_f2 / step) + 1)
        for j in range(n_f2_steps):
            f2 = j * step
            if f1 + f2 > 1.0 + 1e-9:
                break
            fracs = [f1, f2]
            ll = total_log_likelihood_multi(
                markers, fracs, error_rate, marker_biases
            )
            if ll > best_ll:
                best_ll = ll
                best_f = fracs[:]

    # Step 2: Nelder-Mead refinement
    from scipy.optimize import minimize

    def neg_ll(x):
        f1, f2 = x
        # Enforce constraints via penalty
        if f1 < 0 or f2 < 0 or f1 + f2 > 1.0:
            return 1e30
        return -total_log_likelihood_multi(
            markers, [f1, f2], error_rate, marker_biases
        )

    result = minimize(
        neg_ll,
        x0=best_f,
        method="Nelder-Mead",
        options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": 1000},
    )

    f_mle = [max(0.0, x) for x in result.x]
    # Clamp to simplex
    if sum(f_mle) > 1.0:
        scale = 1.0 / sum(f_mle)
        f_mle = [f * scale for f in f_mle]
    ll_max = -result.fun

    # Step 3: Profile likelihood CIs (see section 3c below)
    cis = _profile_likelihood_cis_multi(
        markers, f_mle, ll_max, error_rate, marker_biases
    )

    # Step 4: Per-marker residuals (see section 3d below)
    per_marker = _compute_per_marker_multi(markers, f_mle, error_rate, marker_biases)

    return MultiDonorResult(
        donor_fractions=f_mle,
        donor_fraction_cis=cis,
        host_fraction=1.0 - sum(f_mle),
        log_likelihood=ll_max,
        n_informative=n_informative,
        n_markers_used=sum(1 for m in per_marker if m.included),
        per_marker=per_marker,
        error_rate=error_rate,
    )
```

### 3c. Profile likelihood CIs for 2D case

For each donor, we construct a profile likelihood CI by fixing that donor's fraction and re-optimizing the other donor's fraction. This gives a 1D profile for each donor.

**Important:** The chi-squared threshold uses **df=1** (not df=2) because we are profiling over one parameter at a time. The df=2 threshold (5.99) would be used for a joint confidence region — but we want marginal CIs.

```python
def _profile_likelihood_cis_multi(
    markers: list[InformativeMarker],
    f_mle: list[float],
    ll_max: float,
    error_rate: float,
    marker_biases: dict | None,
) -> list[tuple[float, float]]:
    """Profile likelihood CIs for each donor fraction.

    For donor_i, scan f_i while optimizing f_j (j != i) at each point.
    The threshold is chi2(df=1) because we profile one parameter at a time.

    Returns:
        List of (ci_lo, ci_hi) tuples, one per donor.
    """
    threshold = chi2.ppf(0.95, df=1)  # ~3.84
    half_threshold = threshold / 2.0
    step = 0.001
    cis = []

    for donor_idx in range(len(f_mle)):

        def profile_ll(fi: float) -> float:
            """Max LL over other donors, with donor_idx fixed at fi."""
            # For 2 donors, the other donor's fraction is optimized in [0, 1-fi]
            other_idx = 1 - donor_idx  # works for 2 donors
            from scipy.optimize import minimize_scalar

            def neg_ll_other(fj):
                fracs = [0.0, 0.0]
                fracs[donor_idx] = fi
                fracs[other_idx] = fj
                return -total_log_likelihood_multi(
                    markers, fracs, error_rate, marker_biases
                )

            res = minimize_scalar(
                neg_ll_other,
                bounds=(0.0, max(0.0, 1.0 - fi)),
                method="bounded",
            )
            return -res.fun

        # Scan left from MLE
        f_lo = f_mle[donor_idx]
        while f_lo > 0.0:
            f_test = max(0.0, f_lo - step)
            ll_test = profile_ll(f_test)
            if (ll_max - ll_test) > half_threshold:
                f_lo = f_test
                break
            f_lo = f_test
            if f_test == 0.0:
                break

        # Scan right from MLE
        f_hi = f_mle[donor_idx]
        while f_hi < 1.0:
            f_test = min(1.0, f_hi + step)
            ll_test = profile_ll(f_test)
            if (ll_max - ll_test) > half_threshold:
                f_hi = f_test
                break
            f_hi = f_test
            if f_test == 1.0:
                break

        cis.append((f_lo, f_hi))

    return cis
```

### 3d. Result types

```python
@dataclass
class MultiDonorResult:
    """Result of multi-donor chimerism estimation."""

    donor_fractions: list[float]              # [f_donor1, f_donor2, ...]
    donor_fraction_cis: list[tuple[float, float]]  # [(lo, hi), ...] per donor
    host_fraction: float                       # 1 - sum(donor_fractions)
    log_likelihood: float
    n_informative: int
    n_markers_used: int
    per_marker: list[MarkerResult]
    error_rate: float
```

### 3e. Backward compatibility: single-donor as special case

The CLI will auto-detect the number of donors from `--donor` arguments:
- 1 donor → call existing `estimate_single_donor()` (unchanged)
- 2 donors → call new `estimate_multi_donor()`

This avoids any regression risk for the existing single-donor path.

---

## Part 4: Genotype Module Changes

### 4a. Per-donor marker informativity

Currently `classify_markers()` only computes `marker_type` for the first donor. For multi-donor, we need to track which donors each marker is informative for.

```python
# Changes to genotype.py

@dataclass
class InformativeMarker:
    """A marker where host and at least one donor have different genotypes."""

    chrom: str
    pos: int
    ref: str
    alt: str
    host_gt: tuple[int, int]
    donor_gts: list[tuple[int, int]]
    marker_type: int               # Vynck type for FIRST donor (backward compat)
    admix_ad_ref: int
    admix_ad_alt: int
    admix_dp: int
    # New fields for multi-donor:
    marker_types: list[int | None] | None = None  # Vynck type per donor, None if non-informative
    informative_for: list[bool] | None = None      # True per donor if informative
```

### 4b. Updated classify_markers

The informativity check in `classify_markers()` currently only considers the first donor. Change to: a marker is informative if the host differs from **any** donor.

```python
# In classify_markers():

# Current (line 258):
mtype = marker_type(h.gt, ds[0].gt)

# Change to:
mtypes = [marker_type(h.gt, d.gt) for d in ds]
any_informative = any(mt is not None for mt in mtypes)
mtype_first = mtypes[0]  # backward compat

if any_informative:
    informative.append(InformativeMarker(
        chrom=key[0],
        pos=key[1],
        ref=key[2],
        alt=key[3],
        host_gt=h.gt,
        donor_gts=donor_gts,
        marker_type=mtype_first if mtype_first is not None else mtypes[1],
        admix_ad_ref=a.ad_ref,
        admix_ad_alt=a.ad_alt,
        admix_dp=a.dp,
        marker_types=mtypes,
        informative_for=[mt is not None for mt in mtypes],
    ))
```

---

## Part 5: Report and CLI Changes

### 5a. Report output for multi-donor

**TSV format** — add per-donor columns:

```
sample  donor1_pct  donor1_ci_lo  donor1_ci_hi  donor2_pct  donor2_ci_lo  donor2_ci_hi  host_pct  n_informative  n_used  mean_depth  gof_pval  qc_pass
```

**JSON format** — structured per-donor:

```json
{
  "sample": "d30_admix",
  "host_pct": 60.0,
  "donors": [
    {
      "label": "donor1",
      "donor_pct": 25.23,
      "ci_lo": 23.10,
      "ci_hi": 27.36,
      "n_informative": 58
    },
    {
      "label": "donor2",
      "donor_pct": 14.77,
      "ci_lo": 12.50,
      "ci_hi": 17.04,
      "n_informative": 62
    }
  ],
  "total_donor_pct": 40.0,
  "n_informative": 72,
  "n_used": 71,
  "mean_depth": 1045.3,
  "gof_pval": 0.2341,
  "qc_pass": true
}
```

### 5b. CLI changes

```python
# In cli.py _run_single_sample():

if len(donor_paths) == 1:
    result = estimate_single_donor(genotypes.informative, ...)
else:
    result = estimate_multi_donor(genotypes.informative, n_donors=len(donor_paths), ...)
```

The report functions (`to_tsv`, `to_json`, `timeline_json`) need multi-donor variants that accept `MultiDonorResult`. Use duck typing or a union type to handle both result types.

### 5c. Timeline output for multi-donor

```json
{
  "timepoints": [
    {
      "sample": "d30",
      "host_pct": 60.0,
      "donors": [
        {"label": "donor1", "donor_pct": 25.2, "ci_lo": 23.1, "ci_hi": 27.4},
        {"label": "donor2", "donor_pct": 14.8, "ci_lo": 12.5, "ci_hi": 17.0}
      ]
    },
    {
      "sample": "d60",
      "host_pct": 35.0,
      "donors": [
        {"label": "donor1", "donor_pct": 42.1, "ci_lo": 39.8, "ci_hi": 44.4},
        {"label": "donor2", "donor_pct": 22.9, "ci_lo": 20.4, "ci_hi": 25.4}
      ]
    }
  ]
}
```

---

## Part 6: QC Changes

Per-donor informative marker counts need to be tracked. The QC report should include:

```python
@dataclass
class MultiDonorQCReport(QCReport):
    """Extended QC for multi-donor results."""
    per_donor_n_informative: list[int] | None = None  # informative markers per donor
```

A new warning condition: if any single donor has fewer than `min_informative` markers, warn that that donor's estimate may be unreliable.

---

## Part 7: Test Data Generation — Three Brothers Scenario

### Script: `scripts/generate_multidonor_test_data.py`

Generates synthetic VCFs for the 3-brother scenario at a grid of (f₁, f₂) mixture points.

```python
#!/usr/bin/env python3
"""Generate multi-donor test data: 3 brothers (host + 2 donors).

Creates:
    - host.vcf, donor1.vcf, donor2.vcf (sibling genotype VCFs)
    - Chimeric VCFs at a grid of (f1, f2) mixture fractions
    - truth_table.tsv with ground truth fractions

Usage:
    python scripts/generate_multidonor_test_data.py --outdir tests/test_data/multidonor
"""

import argparse
import csv
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix.simulate import (
    generate_sibling_trio_genotypes,
    write_genotype_vcf,
    blend_from_genotype_dicts,
    write_vcf,
)


# Mixture fraction grid points on the simplex f1 + f2 <= 1
# Includes edges (single-donor cases) and interior points
MIXTURE_GRID = [
    # Pure host
    (0.00, 0.00),
    # Single-donor edges (test backward compat)
    (0.05, 0.00), (0.20, 0.00), (0.50, 0.00),
    (0.00, 0.05), (0.00, 0.20), (0.00, 0.50),
    # Balanced two-donor
    (0.10, 0.10),
    (0.25, 0.25),
    (0.40, 0.40),
    # Asymmetric two-donor
    (0.30, 0.10),
    (0.10, 0.30),
    (0.50, 0.20),
    (0.20, 0.50),
    (0.05, 0.15),
    (0.15, 0.05),
    # High total donor (near relapse)
    (0.45, 0.45),
    (0.60, 0.30),
    # Low-fraction detection
    (0.02, 0.02),
    (0.01, 0.05),
    # Pure donor1
    (1.00, 0.00),
    # Pure donor2
    (0.00, 1.00),
]


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="tests/test_data/multidonor")
    parser.add_argument("--n-markers", type=int, default=100)
    parser.add_argument("--depth", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # Generate 3-sibling genotypes
    markers = generate_sibling_trio_genotypes(args.n_markers, rng)

    # Write genotype VCFs
    write_genotype_vcf(markers, outdir / "host.vcf", "HOST", key="host_gt")
    write_genotype_vcf(markers, outdir / "donor1.vcf", "DONOR1", key="donor1_gt")
    write_genotype_vcf(markers, outdir / "donor2.vcf", "DONOR2", key="donor2_gt")

    # Report informativity
    n_inf_d1 = sum(1 for m in markers if m["informative_d1"])
    n_inf_d2 = sum(1 for m in markers if m["informative_d2"])
    n_inf_any = sum(1 for m in markers if m["informative_any"])
    n_distinguishable = sum(1 for m in markers if m["donors_distinguishable"])
    print(f"Markers: {args.n_markers}", file=sys.stderr)
    print(f"  Informative for donor1: {n_inf_d1}", file=sys.stderr)
    print(f"  Informative for donor2: {n_inf_d2}", file=sys.stderr)
    print(f"  Informative for any: {n_inf_any}", file=sys.stderr)
    print(f"  Donors distinguishable: {n_distinguishable}", file=sys.stderr)

    # Generate chimeric VCFs at grid points
    truth_rows = []
    for f1, f2 in MIXTURE_GRID:
        name = f"host_{100-int((f1+f2)*100)}_d1_{int(f1*100)}_d2_{int(f2*100)}"
        result = blend_from_genotype_dicts(
            markers, [f1, f2],
            target_depth=args.depth,
            seed=rng.randint(0, 2**31),
            error_rate=0.01,
            depth_cv=0.43,
        )
        write_vcf(result, outdir / f"{name}.vcf")
        truth_rows.append({
            "sample_name": name,
            "true_donor1_fraction": f1,
            "true_donor2_fraction": f2,
            "true_host_fraction": 1.0 - f1 - f2,
            "num_markers": result.num_markers,
            "num_informative_any": result.num_informative,
        })

    # Write truth table
    with open(outdir / "truth_table.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=truth_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_rows)

    print(f"\nGenerated {len(MIXTURE_GRID)} chimeric VCFs in {outdir}/", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
```

---

## Part 8: Test Suite

### New test file: `tests/test_multidonor.py`

```python
"""Tests for multi-donor chimerism estimation."""

import math
import random
from pathlib import Path

import pytest

from allomix.chimerism import (
    expected_weight_multi,
    log_likelihood_marker,
    total_log_likelihood_multi,
    estimate_multi_donor,
    MultiDonorResult,
)
from allomix.genotype import InformativeMarker, classify_markers
from allomix.simulate import (
    generate_sibling_trio_genotypes,
    blend_from_genotype_dicts,
    write_genotype_vcf,
    write_vcf,
    expected_vaf_multi,
)


class TestExpectedWeightMulti:
    """Test the 2D expected weight function."""

    def test_pure_host(self):
        """f1=f2=0: weight depends only on host."""
        # Host 0/0 → ref_dose=2, weight=1.0
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.0, 0.0])
        assert w == pytest.approx(1.0)

    def test_pure_donor1(self):
        """f1=1.0, f2=0: weight depends only on donor1."""
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [1.0, 0.0])
        assert w == pytest.approx(0.0)  # donor1 is 1/1, ref_dose=0

    def test_pure_donor2(self):
        """f1=0, f2=1.0: weight depends only on donor2."""
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.0, 1.0])
        assert w == pytest.approx(0.5)  # donor2 is 0/1, ref_dose=1

    def test_equal_mix(self):
        """f1=f2=0.25: weighted average."""
        # Host 0/0 (ref_dose=2), D1 1/1 (ref_dose=0), D2 0/1 (ref_dose=1)
        # w = 0.5 * 2/2 + 0.25 * 0/2 + 0.25 * 1/2 = 0.5 + 0 + 0.125 = 0.625
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.25, 0.25])
        assert w == pytest.approx(0.625)

    def test_reduces_to_single_donor(self):
        """With f2=0, should match single-donor expected_weight."""
        from allomix.chimerism import expected_weight
        for f in [0.0, 0.1, 0.5, 1.0]:
            w_single = expected_weight((0, 0), (1, 1), f)
            w_multi = expected_weight_multi((0, 0), [(1, 1), (0, 0)], [f, 0.0])
            assert w_multi == pytest.approx(w_single)


class TestMultiDonorLikelihood:
    """Test 2D log-likelihood computation."""

    def test_likelihood_peaks_at_truth(self):
        """LL should be maximized near the true fractions."""
        rng = random.Random(42)
        markers_data = generate_sibling_trio_genotypes(100, rng)
        # ... (create InformativeMarker objects at known fractions)
        # ... (verify LL at truth > LL at wrong fractions)


class TestEstimateMultiDonor:
    """End-to-end multi-donor estimation tests."""

    @pytest.fixture
    def sibling_markers(self):
        """Generate 100 markers for 3 siblings."""
        return generate_sibling_trio_genotypes(100, random.Random(42))

    def test_pure_host(self, sibling_markers, tmp_path):
        """f1=f2=0: both donors should estimate near 0%."""
        # ... generate chimeric VCF at (0, 0), run estimation
        # assert result.donor_fractions[0] < 0.03
        # assert result.donor_fractions[1] < 0.03

    def test_single_donor_only(self, sibling_markers, tmp_path):
        """f1=0.20, f2=0: should recover donor1=20%, donor2≈0%."""
        # ... test that multi-donor correctly handles single-donor case

    def test_balanced_mix(self, sibling_markers, tmp_path):
        """f1=f2=0.25: should recover both donors near 25%."""
        # ... test balanced two-donor mixture

    def test_asymmetric_mix(self, sibling_markers, tmp_path):
        """f1=0.30, f2=0.10: should distinguish major vs minor donor."""
        # ... key test: can we tell which donor contributes more?

    def test_low_fraction_detection(self, sibling_markers, tmp_path):
        """f1=0.05, f2=0.02: can we detect small contributions?"""
        # ... sensitivity test at low fractions

    def test_ci_contains_truth(self, sibling_markers, tmp_path):
        """Profile likelihood CIs should contain the true fractions."""
        # ... test CI coverage across several mixture points

    def test_fractions_sum_le_one(self, sibling_markers, tmp_path):
        """Estimated fractions must satisfy f1 + f2 <= 1."""
        # ... constraint enforcement test
```

---

## Part 9: Paper Figure — Multi-Donor Estimation with Sibling Donors

### Figure concept: 3-panel figure

This figure demonstrates multi-donor chimerism estimation accuracy using the 3-sibling scenario, which is the clinically hardest case (maximum genotype sharing between related donors).

**Panel A — Ternary simplex plot** showing true vs estimated mixture compositions for the full grid of test points. Each point is a mixture of (host%, donor1%, donor2%). Arrows from true to estimated positions visualize the error vector. This is the most natural way to show 3-component mixtures and immediately communicates the constraint that fractions sum to 100%.

**Panel B — Per-donor accuracy scatter** with donor1 on the left and donor2 on the right. Each panel shows true fraction (x-axis) vs estimated fraction (y-axis) with identity line. Error bars show profile likelihood 95% CIs. This is the most quantitatively readable format.

**Panel C — 2D log-likelihood contour** for one representative mixture point (e.g., 50% host, 30% donor1, 20% donor2). Shows the log-likelihood surface as a heatmap/contour in (f₁, f₂) space with the MLE point and 95% CI contour overlaid. The triangular feasible region (f₁ + f₂ ≤ 1) is clearly visible.

### Figure generation script

```python
#!/usr/bin/env python3
"""Generate multi-donor validation figure for the allomix paper.

Produces a 3-panel figure:
  A: Ternary simplex showing true vs estimated compositions
  B: Per-donor accuracy scatter with CIs
  C: 2D likelihood contour for one example mixture

Usage:
    python scripts/generate_multidonor_figure.py \
        --truth tests/test_data/multidonor/truth_table.tsv \
        --results output/multidonor_validation/results.tsv \
        --output paper/figures/fig_multidonor.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


def ternary_to_cartesian(host, d1, d2):
    """Convert ternary coordinates to 2D Cartesian for plotting.

    Uses the standard equilateral triangle with:
    - Host at bottom-left (0, 0)
    - Donor1 at bottom-right (1, 0)
    - Donor2 at top (0.5, sqrt(3)/2)
    """
    x = d1 + d2 * 0.5
    y = d2 * np.sqrt(3) / 2
    return x, y


def plot_panel_a(ax, truth_data, estimated_data):
    """Ternary simplex: true vs estimated positions."""
    # Draw triangle outline
    triangle = plt.Polygon(
        [(0, 0), (1, 0), (0.5, np.sqrt(3)/2)],
        fill=False, edgecolor="black", linewidth=1.5
    )
    ax.add_patch(triangle)

    # Plot grid lines inside triangle
    for frac in [0.2, 0.4, 0.6, 0.8]:
        # ... draw iso-fraction lines

    # Plot true positions (open circles) and estimated (filled circles)
    for true, est in zip(truth_data, estimated_data):
        tx, ty = ternary_to_cartesian(true["host"], true["d1"], true["d2"])
        ex, ey = ternary_to_cartesian(est["host"], est["d1"], est["d2"])
        ax.plot(tx, ty, "o", color="steelblue", markersize=6,
                markerfacecolor="none", markeredgewidth=1.5)
        ax.plot(ex, ey, "o", color="firebrick", markersize=4)
        ax.annotate("", xy=(ex, ey), xytext=(tx, ty),
                     arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, np.sqrt(3)/2 + 0.05)
    ax.set_aspect("equal")
    ax.axis("off")

    # Corner labels
    ax.text(0, -0.04, "Host 100%", ha="center", fontsize=9)
    ax.text(1, -0.04, "Donor 1\n100%", ha="center", fontsize=9)
    ax.text(0.5, np.sqrt(3)/2 + 0.04, "Donor 2\n100%", ha="center", fontsize=9)
    ax.set_title("A", fontsize=14, fontweight="bold", loc="left")


def plot_panel_b(ax_left, ax_right, truth_data, estimated_data):
    """Per-donor accuracy scatter with CIs."""
    for ax, donor_idx, label in [(ax_left, 0, "Donor 1"), (ax_right, 1, "Donor 2")]:
        true_vals = [t[f"d{donor_idx+1}"] * 100 for t in truth_data]
        est_vals = [e[f"d{donor_idx+1}"] * 100 for e in estimated_data]
        ci_los = [e[f"d{donor_idx+1}_ci_lo"] * 100 for e in estimated_data]
        ci_his = [e[f"d{donor_idx+1}_ci_hi"] * 100 for e in estimated_data]

        ax.errorbar(true_vals, est_vals,
                     yerr=[np.array(est_vals) - np.array(ci_los),
                           np.array(ci_his) - np.array(est_vals)],
                     fmt="o", color="steelblue", markersize=5,
                     capsize=2, elinewidth=0.8)
        ax.plot([0, 100], [0, 100], "k--", alpha=0.4)
        ax.set_xlabel(f"True {label} %")
        ax.set_ylabel(f"Estimated {label} %")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)


def plot_panel_c(ax, markers, true_f1, true_f2):
    """2D log-likelihood contour for one mixture."""
    from allomix.chimerism import total_log_likelihood_multi

    # Evaluate LL on a grid
    n_grid = 200
    f1_range = np.linspace(0, 1, n_grid)
    f2_range = np.linspace(0, 1, n_grid)
    ll_grid = np.full((n_grid, n_grid), np.nan)

    for i, f1 in enumerate(f1_range):
        for j, f2 in enumerate(f2_range):
            if f1 + f2 <= 1.0:
                ll_grid[j, i] = total_log_likelihood_multi(
                    markers, [f1, f2], error_rate=0.01
                )

    ll_max = np.nanmax(ll_grid)
    # Plot as delta-LL from maximum
    delta_ll = ll_max - ll_grid

    # Contour levels: 1.92 (95% CI for 2df joint region)
    contour = ax.contour(
        f1_range * 100, f2_range * 100, delta_ll,
        levels=[0.5, 1.92, 4.0, 8.0],
        colors=["firebrick", "steelblue", "gray", "lightgray"],
    )
    ax.clabel(contour, fmt={0.5: "0.5", 1.92: "95%", 4.0: "4", 8.0: "8"})

    # Fill the infeasible region
    ax.fill_between(
        f1_range * 100, (1 - f1_range) * 100, 100,
        color="lightgray", alpha=0.3
    )

    # Mark true and estimated positions
    ax.plot(true_f1 * 100, true_f2 * 100, "k*", markersize=12, label="Truth")
    # ... add MLE point

    ax.set_xlabel("Donor 1 fraction (%)")
    ax.set_ylabel("Donor 2 fraction (%)")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal")
    ax.set_title("C", fontsize=14, fontweight="bold", loc="left")


def main():
    fig = plt.figure(figsize=(14, 5))

    # Panel A: ternary (left)
    ax_a = fig.add_subplot(1, 3, 1)
    # ... load data, call plot_panel_a

    # Panel B: per-donor scatter (center, split into two sub-panels)
    ax_b1 = fig.add_subplot(1, 3, 2)
    # ... use gridspec for two sub-panels

    # Panel C: likelihood contour (right)
    ax_c = fig.add_subplot(1, 3, 3)
    # ... load one example, call plot_panel_c

    fig.tight_layout()
    fig.savefig("paper/figures/fig_multidonor.png", dpi=300, bbox_inches="tight")
```

### Alternative: simpler 2-panel figure

If the ternary plot adds complexity without proportional insight (the reviewer may not be familiar with ternary diagrams), a simpler alternative:

**Panel A — Per-donor scatter** (donor1 and donor2 on same axes with different colors/shapes)
**Panel B — 2D likelihood contour** for one example mixture

This is more conventional and easier to read at a glance.

---

## Part 10: Implementation Order

| Step | File | Status |
|------|------|--------|
| 1 | **simulate.py** — `generate_sibling_trio_genotypes()`, `_mendelian_child()`, `expected_vaf_multi()`, `blend_from_genotype_dicts()` | ✅ Done |
| 2 | **genotype.py** — `marker_types`, `informative_for` fields; `classify_markers()` any-donor informativity | ✅ Done |
| 3 | **chimerism.py** — `MultiDonorResult`, `expected_weight_multi()`, `total_log_likelihood_multi()`, `estimate_multi_donor()`, profile CIs | ✅ Done |
| 4 | **qc.py** — `per_donor_n_informative`, per-donor CI/informativity warnings | ✅ Done |
| 5 | **report.py** — `_write_tsv_multi()`, multi-donor `to_json()`, `timeline_json()` | ✅ Done |
| 6 | **cli.py** — Auto-detect donor count in `_run_single_sample()` | ✅ Done |
| 7 | **scripts/generate_multidonor_test_data.py** — 3-brothers test data | ✅ Done |
| 8 | **tests/test_multidonor.py** — 48 tests (unit + integration + CLI) | ✅ Done |
| 9 | **paper/scripts/run_multidonor_validation.py** — Systematic validation + facts | 📝 TODO |
| 10 | **paper/scripts/generate_multidonor_figure.py** — Paper figure | 📝 TODO |
| 11 | **paper/** — Update methods, results, discussion, intro, abstract | 📝 TODO |
| 12 | **README.md** — Remove "not yet implemented" line | 📝 TODO |

### Test matrix (all passing)

| (f₁, f₂) | Description | Key assertion | Status |
|---|---|---|---|
| (0.00, 0.00) | Pure host | Both donors < 3% | ✅ |
| (0.20, 0.00) | Single donor only | d1 ≈ 20%, d2 < 8% | ✅ |
| (0.25, 0.25) | Balanced | Both ≈ 25% ± 10% | ✅ |
| (0.30, 0.10) | Asymmetric | d1 > d2 | ✅ |
| (0.10, 0.30) | Asymmetric flipped | d2 > d1 | ✅ |
| (0.25, 0.15) | CI coverage | CIs contain truth | ✅ |
| (0.40, 0.40) | Fractions constraint | f1+f2 ≤ 1 | ✅ |
| (1.00, 0.00) | Pure donor1 | d1 > 90% | ✅ |
| (0.00, 1.00) | Pure donor2 | d2 > 90% | ✅ |

### Performance considerations

- Grid search at 101 steps: 101 × 51 ≈ 5,151 evaluations (each sums ~61 markers) — <0.5s
- Nelder-Mead: typically converges in <100 iterations — <0.1s
- Profile likelihood CI: 2 × ~2000 evaluations (scanning + inner optimization) — <2s per donor
- **Total: ~5s per sample** — acceptable for clinical use

---

## References

- Crysup & Woerner (2022) — MLE likelihood framework for mixture deconvolution
- De Vynck et al. (2023) — Marker type classification and bias correction
- Wilks (1938) — Profile likelihood CI theory
- FABCASE `copy_num2_2donor_calculate.pl` — prior art for 2-donor likelihood (different implementation)
