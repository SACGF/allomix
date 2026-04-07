# Step 5: Detailed Implementation Plan

Synthesises Steps 2 (VCF-first), 3 (test data), and 4 (reference tools) into a file-level plan for allomix v1.

---

## Architecture Overview

```
VCF (host)  ─┐
VCF (donor) ─┤── genotype.py ──> MarkerGenotypes
VCF (admix) ─┘       │
                      ▼
               chimerism.py ──> ChimerismResult
                      │              │
                      ▼              ▼
                 qc.py          report.py ──> TSV / JSON
                                     │
                                     ▼
                                 cli.py
```

allomix is a pipeline of four stages:

1. **Parse & compare genotypes** — read VCFs, match markers, classify informativeness
2. **Estimate chimerism** — MLE on allele counts at informative markers
3. **QC** — flag outlier markers, assess confidence
4. **Report** — output results in TSV/JSON, optional timeline

---

## Decision Summary

| Decision | Choice | Reasoning |
|----------|--------|-----------|
| Primary input | VCF (GT, AD, DP required) | Pipeline is controlled; joint calling provides ref+alt AD at all sites (Step 2) |
| Core algorithm | MLE with known genotypes (Demixtify Formula 5 simplified) | Proper statistical framework; known genotypes eliminate marginalisation (Step 4) |
| Bias correction | Per-marker amplification bias (Vynck et al.) | Single biggest accuracy improvement; ~30% error reduction (Step 4) |
| Error model | Sequencing error rate `e` in likelihood; estimate empirically or default 0.01 | Standard approach from Demixtify/Conpair (Step 4) |
| Multi-donor | 2D grid search over (f_h, f_d1), Nelder-Mead refinement | Extends single-donor MLE naturally (Step 4) |
| Confidence intervals | Profile likelihood, chi-square approximation | Standard; used by Demixtify and Conpair (Step 4) |
| Test data | Synthetic VCF blending via binomial sampling | Already built in Step 3 (`simulate.py`) |

---

## Module Breakdown

### 1. `src/allomix/genotype.py` — Genotype Parsing & Comparison

**Responsibility:** Read VCFs, extract genotypes and allele counts at shared loci, classify markers.

**Key types:**

```python
@dataclass
class MarkerData:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int, int]      # (0,0), (0,1), or (1,1)
    ad_ref: int               # ref allele depth
    ad_alt: int               # alt allele depth
    dp: int                   # total depth
    gq: int | None            # genotype quality (optional)

@dataclass
class InformativeMarker:
    chrom: str
    pos: int
    ref: str
    alt: str
    host_gt: tuple[int, int]
    donor_gts: list[tuple[int, int]]   # 1 or 2 donors
    marker_type: int                    # Vynck classification (0,1,10,11,20,21)
    admix_ad_ref: int
    admix_ad_alt: int
    admix_dp: int
    bias: float | None                 # per-marker bias estimate (if available)

@dataclass
class MarkerGenotypes:
    informative: list[InformativeMarker]
    non_informative: list[MarkerData]   # for QC / bias estimation
    n_total: int                         # total markers in input
    n_shared: int                        # markers present in all VCFs
    n_filtered: int                      # markers excluded by QC
    sample_name: str
```

**Functions:**

```python
def parse_vcf(path: Path) -> list[MarkerData]:
    """Read a VCF and extract MarkerData at each record. Uses cyvcf2."""

def find_shared_markers(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admixture: list[MarkerData],
) -> list[tuple[MarkerData, ...]]:
    """Join markers by (chrom, pos, ref, alt). Return only shared loci."""

def classify_markers(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admixture: list[MarkerData],
    min_dp: int = 100,
    min_gq: int = 20,
) -> MarkerGenotypes:
    """Classify shared markers as informative/non-informative.
    Apply depth and quality filters. Assign Vynck marker types."""

def marker_type(host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> int | None:
    """Classify as Vynck type 0,1,10,11,20,21 or None if non-informative."""
```

**Design notes:**
- Uses cyvcf2 for VCF parsing (lightweight, fast for read-only).
- Markers joined on (chrom, pos, ref, alt) — exact match required.
- Multi-donor: a marker is informative for donor_i if host and donor_i differ. A marker can be informative for one donor but not the other.
- Filters: min depth (default 100x), min GQ (default 20), PASS filter only, no-call excluded.

---

### 2. `src/allomix/chimerism.py` — Core MLE Algorithm

**Responsibility:** Estimate chimerism fraction(s) from informative markers via maximum likelihood.

**Key types:**

```python
@dataclass
class ChimerismResult:
    donor_fraction: float              # MLE point estimate (0.0–1.0)
    donor_fraction_ci: tuple[float, float]  # 95% CI
    host_fraction: float               # 1 - sum(donor fractions)
    log_likelihood: float              # at MLE
    n_informative: int
    n_markers_used: int                # after outlier exclusion
    per_marker: list[MarkerResult]
    error_rate: float                  # estimated or default

@dataclass
class MultiDonorResult:
    donor_fractions: list[float]       # one per donor
    donor_fraction_cis: list[tuple[float, float]]
    host_fraction: float
    log_likelihood: float
    n_informative: list[int]           # per donor
    per_marker: list[MarkerResult]

@dataclass
class MarkerResult:
    chrom: str
    pos: int
    marker_type: int
    expected_vaf: float
    observed_vaf: float
    residual: float
    ad_ref: int
    ad_alt: int
    dp: int
    included: bool                     # False if outlier-excluded
    bias_corrected: bool
```

**Functions:**

```python
def expected_weight(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    f_donor: float,
) -> float:
    """Expected reference allele weight for a given chimerism fraction.
    w = (1-f) * host_ref_dose/2 + f * donor_ref_dose/2"""

def log_likelihood_marker(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
) -> float:
    """Per-marker log-likelihood (Demixtify Formula 5 with known genotypes).
    LL = n_ref * log(w*(1-e) + (1-w)*e/3) + n_alt * log((1-w)*(1-e) + w*e/3)"""

def total_log_likelihood(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
) -> float:
    """Sum of per-marker log-likelihoods."""

def estimate_single_donor(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
) -> ChimerismResult:
    """MLE for single-donor chimerism.
    1. Grid search over f in [0, 1] at 0.1% steps
    2. Brent refinement around grid maximum
    3. Profile likelihood CI (chi-sq 1df, alpha=0.05)
    4. Per-marker residuals and outlier flagging"""

def estimate_two_donors(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 201,
) -> MultiDonorResult:
    """MLE for two-donor chimerism.
    2D grid search over (f_d1, f_d2) with constraint f_d1+f_d2 <= 1.
    Nelder-Mead refinement. Profile likelihood CIs."""

def estimate_error_rate(markers: list[InformativeMarker]) -> float:
    """Empirical error rate from non-informative markers or
    from 'other' allele observations. Default fallback: 0.01."""
```

**MLE detail:**

For single donor, the likelihood is:

```
LL(f) = Σ_i [ n_ref_i * log(w_i(f)*(1-e) + (1-w_i(f))*e/3)
            + n_alt_i * log((1-w_i(f))*(1-e) + w_i(f)*e/3) ]
```

where `w_i(f) = (1-f) * host_ref_dose_i/2 + f * donor_ref_dose_i/2`.

Optimisation: Grid search at 0.1% resolution (1001 points) to find the approximate maximum, then `scipy.optimize.minimize_scalar(method='bounded')` (Brent) in a ±1% window for refinement.

CI: Find f_lo, f_hi where `2*(LL_max - LL(f)) = chi2.ppf(0.95, df=1) ≈ 3.84`.

---

### 3. `src/allomix/bias.py` — Per-Marker Bias Correction

**Responsibility:** Estimate and correct per-marker amplification bias (Vynck et al.).

```python
def estimate_bias(
    het_vafs: list[float],
) -> float:
    """Estimate bias for a single marker as median(VAF - 0.5)
    across heterozygous samples in a training set."""

def correct_chimerism(
    raw_host_pct: float,
    marker_type: int,
    bias: float,
) -> float:
    """Apply Vynck analytic correction for marker type and bias.
    E.g. type 0: corrected = (raw - 2*raw*b) / (-4*b*raw + 2*b + 1)"""

def load_bias_table(path: Path) -> dict[tuple[str, int], float]:
    """Load a pre-computed bias table (marker_id -> bias).
    TSV with columns: chrom, pos, bias"""

def compute_bias_table(
    training_vcfs: list[Path],
    markers: list[tuple[str, int]],
) -> dict[tuple[str, int], float]:
    """Compute per-marker bias from a set of training VCFs.
    For each marker, collect VAFs from het samples, compute median deviation."""
```

**Design notes:**
- Bias correction is optional — it requires a training set of pure-sample VCFs.
- If no bias table is provided, allomix runs without correction (acceptable for initial deployment).
- Bias estimation can be run as a standalone CLI subcommand: `allomix estimate-bias --vcfs *.vcf.gz --output bias_table.tsv`
- The MLE in `chimerism.py` can incorporate bias either:
  - (a) Pre-correct AD counts before MLE (simpler, what we'll do in v1)
  - (b) Incorporate bias into the likelihood function (more principled, future enhancement)

---

### 4. `src/allomix/qc.py` — Quality Control

**Responsibility:** Assess result quality, flag issues.

```python
@dataclass
class QCReport:
    n_total_markers: int
    n_shared_markers: int
    n_informative: int
    n_used: int                        # after depth/quality/outlier filtering
    n_excluded_depth: int
    n_excluded_quality: int
    n_excluded_outlier: int
    mean_depth: float
    median_depth: float
    min_depth: int
    goodness_of_fit: float             # chi-squared p-value
    warnings: list[str]
    pass_: bool

def assess_quality(
    result: ChimerismResult,
    markers: MarkerGenotypes,
    min_informative: int = 3,
    min_depth: int = 100,
) -> QCReport:
    """Compute QC metrics and raise warnings."""
```

**QC checks:**
- `n_informative < 3` → FAIL ("Insufficient informative markers")
- `mean_depth < min_depth` → WARN
- Goodness-of-fit: compare observed VAFs to model-predicted VAFs using chi-squared test. If p < 0.01, WARN ("Poor model fit — possible genotyping error, CNV, or sample issue")
- Per-marker residual > 3 SD → flag as outlier
- Host fraction < 0 or > 1 → FAIL ("Impossible result")
- CI width > 20% → WARN ("Wide confidence interval")

---

### 5. `src/allomix/report.py` — Output Formatting

**Responsibility:** Format results as TSV, JSON, or summary text.

```python
def to_tsv(result: ChimerismResult, qc: QCReport, output: Path | TextIO) -> None:
    """Write results as TSV. Two sections:
    1. Summary line: sample, donor_pct, ci_lo, ci_hi, n_markers, qc_pass
    2. Per-marker detail lines (if verbose)"""

def to_json(result: ChimerismResult, qc: QCReport) -> dict:
    """Structured JSON output for VariantGrid integration."""

def timeline_json(
    results: list[tuple[str, ChimerismResult, QCReport]],
) -> dict:
    """Timeline of chimerism across timepoints.
    List of {sample, date, donor_pct, ci_lo, ci_hi, qc_pass}."""
```

**TSV format (summary):**
```
sample	donor_pct	ci_lo	ci_hi	n_informative	n_used	mean_depth	gof_pval	qc_pass
day30	12.34	11.02	13.71	42	40	1850	0.45	PASS
```

**TSV format (per-marker detail):**
```
chrom	pos	marker_type	host_gt	donor_gt	ad_ref	ad_alt	observed_vaf	expected_vaf	residual	included
chr1	87923161	0	0/0	1/1	1752	248	0.124	0.123	0.001	True
```

---

### 6. `src/allomix/cli.py` — Command-Line Interface (update existing)

**Subcommands:**

```
allomix monitor   --host H --donor D [--donor D2] --sample S [-o out.tsv]
allomix timeline  --host H --donor D --sample S1 --sample S2 ... [-o out.json]
allomix estimate-bias --vcfs *.vcf.gz -o bias_table.tsv
allomix --version
```

**Common options:**
- `--min-dp N` — minimum depth filter (default: 100)
- `--min-gq N` — minimum genotype quality (default: 20)
- `--error-rate F` — sequencing error rate (default: auto-estimate, fallback 0.01)
- `--bias-table PATH` — per-marker bias correction table
- `--format {tsv,json}` — output format (default: tsv)
- `--verbose` — include per-marker detail in output

---

## Data Flow

```
allomix monitor --host h.vcf --donor d.vcf --sample s.vcf -o results.tsv

1. genotype.parse_vcf(h.vcf) -> host_markers
   genotype.parse_vcf(d.vcf) -> donor_markers
   genotype.parse_vcf(s.vcf) -> admix_markers

2. genotype.classify_markers(host, [donor], admix) -> MarkerGenotypes
     - join on (chrom, pos, ref, alt)
     - filter by depth, GQ, PASS
     - classify informative + marker type

3. [optional] bias.load_bias_table() -> adjust AD or w

4. chimerism.estimate_error_rate(markers) -> e
   chimerism.estimate_single_donor(markers, e) -> ChimerismResult
     - grid search LL(f) over f ∈ [0, 1]
     - Brent refinement
     - profile likelihood CI
     - per-marker residuals

5. qc.assess_quality(result, markers) -> QCReport

6. report.to_tsv(result, qc, output)
```

---

## Test Plan

### Unit tests (`tests/test_genotype.py`)
- Parse the example VCF (`data/idt_rhampseq_sid_example.vcf`), verify extracted fields
- Parse the joint-called VCF (`data/joint_called_example.vcf`), verify hom-ref sites have 2-element AD
- Marker type classification: all 9 GT×GT combinations
- Shared marker joining with mismatched loci
- Depth and GQ filtering

### Unit tests (`tests/test_chimerism.py`)
- `expected_weight`: all GT combinations at f=0, 0.5, 1.0
- `log_likelihood_marker`: known values (hand-computed)
- `total_log_likelihood`: sum matches individual markers
- `estimate_single_donor` with synthetic data from `simulate.py`:
  - f=0%: result ≈ 0% (pure host)
  - f=50%: result ≈ 50%
  - f=100%: result ≈ 100% (pure donor)
  - f=1%: result within CI of 1%
  - f=5%: result within CI of 5%
- CI width decreases with increasing depth
- CI width decreases with increasing number of markers

### Unit tests (`tests/test_bias.py`)
- `estimate_bias`: known VAFs → expected bias
- `correct_chimerism`: known corrections for each marker type
- Round-trip: bias + correction ≈ original

### Integration tests (`tests/test_chimerism.py`)
- Generate synthetic chimeric VCF at 13 fractions (using `simulate.py`)
- Run full pipeline (genotype → chimerism → QC → report)
- Verify: estimated fraction within CI of true fraction for all 13 levels
- Verify: CI narrows at higher depth
- Verify: QC passes for well-behaved data

### Unit tests (`tests/test_qc.py`)
- Insufficient markers → FAIL
- Low depth → WARN
- Wide CI → WARN
- Good data → PASS

### Unit tests (`tests/test_report.py`)
- TSV round-trip (write → read → verify)
- JSON schema validation
- Timeline with multiple timepoints

---

## Implementation Order

| Phase | Module | Depends on | Deliverable |
|-------|--------|-----------|-------------|
| 1 | `genotype.py` | cyvcf2 | Parse VCFs, classify markers |
| 2 | `chimerism.py` | genotype.py, scipy, numpy | MLE single-donor + CI |
| 3 | `qc.py` | chimerism.py | QC assessment |
| 4 | `report.py` | chimerism.py, qc.py | TSV/JSON output |
| 5 | `cli.py` (update) | all above | Wire up `monitor` subcommand |
| 6 | Integration test | simulate.py + all above | Full pipeline validation |
| 7 | `chimerism.py` (extend) | Phase 2 | Multi-donor (2D MLE) |
| 8 | `bias.py` | genotype.py | Bias estimation + correction |
| 9 | `cli.py` (extend) | all above | `timeline` + `estimate-bias` subcommands |

**Phase 1–6 = working single-donor tool with tests.**
**Phase 7–9 = multi-donor, bias correction, timeline.**

---

## Dependencies (confirmed in pyproject.toml)

| Package | Version | Used for |
|---------|---------|----------|
| cyvcf2 | >=0.30 | VCF parsing (genotype.py) |
| numpy | >=1.24 | Array operations |
| scipy | >=1.10 | `minimize_scalar` (Brent), `chi2.ppf`, stats |
| pytest | >=7.0 | Testing (dev) |
| ruff | >=0.4 | Linting (dev) |

No additional dependencies needed for v1.

---

## VariantGrid Integration Points

allomix's JSON output format is designed to be consumed by VariantGrid:

1. **Genotype database**: VG stores donor/host genotypes. allomix reads them as VCFs (VG exports VCFs).
2. **Results ingest**: allomix JSON output → VG API → chimerism results stored per patient.
3. **Timeline view**: VG queries chimerism results for a patient across timepoints → renders timeline chart.

Exact VG API integration is deferred to Step 10 of the overall plan. allomix's JSON schema should be agreed with the VG team before Phase 4.

---

## What This Plan Does NOT Cover

- BAM input (decided against in Step 2)
- Automated pipeline integration (Step 10)
- Real sample validation (Step 11)
- Publication (Step 12)
- GUI / web interface
- Automated donor-host matching (assumes user provides correct VCFs)
