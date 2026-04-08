# Plan: Calibrate Simulations from Empirical Panel Data

## Context

We measured panel characteristics from 210 joint-called VCFs (18,047 samples) on the 76-SNP IDT rhAmpSeq Sample ID panel. The empirical data is checked in at `paper/empirical_results/` and loaded as vibepaper facts from `output/facts/panel_empirical.csv`.

## Empirical Measurements

| Parameter | Value | Current simulation | Gap |
|---|---|---|---|
| Bias SD | 0.018 | 0.02 | Minor — already close |
| Locus dropout rate | 1.6% | 0% | **Missing** |
| Allele dropout (het deficit) | ~1.6% | 0% | Negligible at high depth |
| Mean depth | 1,732x | 50–1,000x | Sims cover lower depths |
| Depth CV across markers | 0.43 | 0 (uniform) | **Missing** |
| Depth CV within markers | 0.68 | 0 (uniform) | **Missing** |
| Markers with >5% no-call | 1/71 | 0 | Covered by locus dropout |
| Max single-marker bias | 0.10 | drawn from N(0,0.02) | Distribution is heavy-tailed |

## Changes Needed

### 1. Non-uniform depth per marker

**Problem**: Currently all markers get the same `target_depth`. In reality, depth varies substantially across markers (CV=0.43 within a sample, range 285x–2789x at mean 1,732x).

**Approach**: Draw per-marker depth from a log-normal distribution that matches the empirical mean and CV.

```python
def sample_marker_depths(
    n_markers: int,
    mean_depth: int,
    depth_cv: float,
    rng: random.Random,
) -> list[int]:
    """Draw per-marker depths from a log-normal matching empirical CV.
    
    The log-normal is parameterised so that:
      E[X] = mean_depth
      CV[X] = depth_cv
    """
    if depth_cv <= 0:
        return [mean_depth] * n_markers
    # Log-normal parameters from desired mean and CV
    sigma2 = math.log(1 + depth_cv ** 2)
    mu = math.log(mean_depth) - sigma2 / 2
    sigma = math.sqrt(sigma2)
    return [max(1, round(math.exp(rng.gauss(mu, sigma)))) for _ in range(n_markers)]
```

**Integration**: Add `depth_cv` parameter to `blend_vcfs()`. When set, draw per-marker depths instead of using a flat `target_depth`. Default to 0.0 (current behaviour).

```python
# In blend_vcfs():
if depth_cv > 0 and target_depth is not None:
    marker_depths = sample_marker_depths(n_shared, target_depth, depth_cv, rng)
else:
    marker_depths = [target_depth or extract_depth(host_rec) or 1000] * n_shared

# Then per marker:
depth = marker_depths[bias_idx]
```

### 2. Locus dropout

**Status**: Already implemented in `blend_vcfs()` via `locus_dropout_rate` parameter.

**Change needed**: Wire it into the validation scripts. Currently the depth validation and relatedness validation scripts don't pass it.

```python
# In run_depth_validation.py generate_and_run():
result = blend_vcfs(
    ...,
    locus_dropout_rate=0.016,  # empirical value
)
```

**Note**: Dropped markers will produce fewer informative markers, which is realistic. allomix already handles this gracefully (missing markers are just absent from the VCF).

### 3. Heavy-tailed bias distribution

**Problem**: The empirical bias distribution is heavy-tailed — median |bias| is 0.005 but 95th percentile is 0.041 and max is 0.10. A Gaussian with SD=0.018 underestimates the tails (Gaussian 95th pct would be ~0.035).

**Approach**: Use a mixture model or a t-distribution to capture the heavy tails.

```python
def generate_marker_biases_realistic(
    n_markers: int,
    rng: random.Random,
    sd: float = 0.018,
    outlier_frac: float = 0.05,
    outlier_sd: float = 0.08,
) -> list[float]:
    """Generate biases with a heavy-tailed distribution.
    
    95% of markers: N(0, sd)       — typical markers
    5% of markers:  N(0, outlier_sd) — outlier markers with extreme bias
    
    This matches the empirical observation of a few markers with
    |bias| > 0.04 while most are below 0.01.
    """
    biases = []
    for _ in range(n_markers):
        if rng.random() < outlier_frac:
            biases.append(rng.gauss(0, outlier_sd))
        else:
            biases.append(rng.gauss(0, sd))
    return biases
```

**Calibration**: From empirical data:
- 95% of markers: SD ≈ 0.012 (the bulk)
- ~5% of markers: SD ≈ 0.08 (the outliers, including the marker with 0.10 bias)
- Overall SD of the mixture ≈ 0.018 (matches observed)

### 4. Per-marker depth reproducibility (run-to-run)

**Problem**: The same marker tends to get similar depth across runs (it's a property of the primer/probe). Currently simulations treat each marker's depth as independent each time.

**Approach**: Pre-generate a per-marker "efficiency" that scales the target depth, then add per-run noise on top.

```python
def sample_marker_depths_with_efficiency(
    efficiencies: list[float],  # pre-generated, one per marker
    mean_depth: int,
    run_cv: float,  # within-run noise (smaller than between-marker CV)
    rng: random.Random,
) -> list[int]:
    """Depth = efficiency * mean_depth * (1 + noise).
    
    efficiencies are fixed per marker (property of the panel).
    run_cv adds per-run stochastic variation on top.
    """
    depths = []
    for eff in efficiencies:
        noise = rng.gauss(1.0, run_cv) if run_cv > 0 else 1.0
        depths.append(max(1, round(eff * mean_depth * noise)))
    return depths
```

**Priority**: Low. The simpler log-normal approach (#1) captures the main effect. This refinement matters mainly for serial monitoring (same patient, multiple timepoints) where the same markers consistently underperform.

### 5. Allele dropout at low depth

**Status**: Already implemented in `blend_vcfs()` via `allele_dropout_rate`.

**Empirical finding**: At 1,700x depth, ADO is negligible (het ratio 0.984). But at lower depths (50x, 100x) ADO would be more significant. The ADO rate should scale with depth.

**Approach**: Make ADO rate depth-dependent.

```python
def ado_rate_for_depth(depth: int, base_rate: float = 0.3) -> float:
    """Allele dropout probability decreases with depth.
    
    Model: ADO ~ base_rate / sqrt(depth)
    At depth=1: ADO ≈ base_rate (30%)
    At depth=100: ADO ≈ 3%
    At depth=1000: ADO ≈ 1%
    At depth=1700: ADO ≈ 0.7% (close to empirical ~1.6% which includes other effects)
    """
    if depth <= 0:
        return base_rate
    return min(base_rate, base_rate / math.sqrt(depth))
```

**Priority**: Low for this paper. Only matters for very low depth (<100x) simulations.

## Implementation Order

1. **Non-uniform depth** (#1) — biggest impact, straightforward
2. **Locus dropout in validation scripts** (#2) — trivial, just pass existing parameter
3. **Heavy-tailed bias** (#3) — moderate impact on CI coverage analysis
4. **Depth-dependent ADO** (#5) — only if we pursue low-depth validation
5. **Per-marker efficiency** (#4) — defer to serial monitoring work

## Validation

After implementation, re-run `paper/scripts/run_depth_validation.py` and `paper/scripts/run_relatedness_validation.py` and compare metrics vs current results. Key things to check:

- MAE should increase slightly (more realistic noise)
- CI coverage should change (non-uniform depth and heavy-tailed bias affect the model-data mismatch)
- Informative marker counts should decrease slightly (locus dropout removes some)
- Results should remain clinically acceptable (MAE < 2%, sufficient markers for siblings)

## Paper Impact

The methods section already references the empirical panel characterisation. After implementing these changes:

- Update the simulation description to mention non-uniform depth and locus dropout
- Re-run all validation scripts and rebuild the paper (`vibepaper build`)
- The facts CSVs will update automatically, flowing new numbers into the text and tables
- Figures will regenerate with the more realistic simulation
