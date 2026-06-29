# Step 21: Per-marker-type overdispersion (remove the sub-0.5% MLE floor)

Implementation plan for issue #33. **No code is written yet** — this is the design,
the snippets to write, and the test/benchmark/validation sequence. Ship as an
**opt-in** so the shared-rho default stays byte-identical (the paper's fixtures and
validated numbers are untouched until we deliberately flip the default).

## 1. Problem recap

The single-donor beta-binomial MLE carries a near-constant **positive offset of
~0.22 pp** at sub-0.5% host fractions. On the committed SRP434573 semi-synthetic
ladder (`paper/public_data/SRP434573/genotypes_synthetic`, 10 mixes x 5 reps) the
median host estimate is ~0.32% at 0.1% nominal, rising in step to ~0.74% at 0.5%.

Root cause (already diagnosed in the issue and `claude/allomix_overall_plan.md`
Step 28 #5 follow-up): a **single shared `rho`** fit jointly with `f` over all
markers lets the symmetric extra-binomial scatter at **donor-heterozygous** markers
(background VAF ~0.5, where amplification scatter is largest and symmetric) rectify
into a small positive host fraction. **Donor-homozygous** markers (VAF near 0/1)
carry one-sided host signal and are not rectified, which is why the host-presence
test (donor-hom only) stays specific at 0.1%.

It is a *dispersion* effect, not a *mean* effect: the donor-het VAF is already
centered on 0.5, a plain binomial returns 0 there, and per-site error / bias tables
(mean corrections) do not remove it. On the donor-het markers the floor scales with
overdispersion: 0.000 at rho=inf, 0.30 at rho=1000, 0.43 at rho=100, 0.48 at the
fitted rho~71.

## 2. Fix: separate rho per marker class

Fit `f` jointly with a **separate concentration `rho` per marker class** (donor-hom
vs donor-het), profiling each rho independently at every `f`:

```
total_LL(f) = max_{rho_hom} LL(hom_markers, f, rho_hom)
            + max_{rho_het} LL(het_markers, f, rho_het)
```

This down-weights the over-dispersed het markers by their own measured variance. At
low fraction they are near-noise, so their fitted rho is small, their effective
weight is ~zero, and the estimate is unbiased. At high fraction the host signal
pushes their VAF off 0.5, their residual scatter shrinks, their rho rises, and they
regain weight and contribute precision. This is why it beats a hard "use hom-only
below X%" switch: it keeps the het markers and weights them adaptively, avoiding the
CI penalty a hom-only switch pays above ~5%.

Marker class (single donor), from the donor genotype's ALT dose:

```python
donor_alt_dose = m.donor_gts[0][0] + m.donor_gts[0][1]
is_het = donor_alt_dose == 1        # background VAF ~0.5 -> over-dispersed class
is_hom = donor_alt_dose in (0, 2)   # background VAF ~0/1 -> one-sided class
```

(Mirrors the existing `donor_ref_dose = 2 - (m.donor_gts[0][0] + m.donor_gts[0][1])`
at `chimerism.py:801` / `:469`.)

## 3. Where the change lands

Everything is inside `src/allomix/chimerism.py`. The estimator profiles a single rho
in three places today, each becomes two rhos:

| Location | Today (single rho) | New (two rho) |
|---|---|---|
| Grid search over f (`:920-934`) | 1 `minimize_scalar` per f point | 2 (hom, het) per f point |
| Nelder-Mead refine (`:937-951`) | simplex over `(f, log_rho)` | simplex over `(f, log_rho_hom, log_rho_het)` |
| Profile-likelihood CI (`:962-989`) | profile 1 rho at each f | profile 2 rhos at each f |
| `detection_limit` / `fraction_se` (`:763-870`) | scalar rho | per-class rho (hom markers use rho_hom, het use rho_het) |

The plumbing already exists: `total_log_likelihood_bb` / `_total_ll_vec` /
`_ll_from_p_alt` / `_p_alt_for_f` all operate on whatever `_MarkerArrays` they are
handed, so each class's LL is just those functions called on a filtered
`_MarkerArrays`. We build **two** `_MarkerArrays` (one per class) with the existing
`_precompute_marker_arrays`, no new likelihood kernel needed.

## 4. Design decisions (resolved)

1. **Flag surface.** New keyword `marker_type_overdispersion: bool = False` on
   `_estimate_single_donor_bb_core` and `estimate_single_donor_bb`, threaded through
   `analyse_sample` and the CLI as `--marker-type-overdispersion` (default off). A
   boolean opt-in, matching the issue.

2. **Byte-identical default.** Branch at the top of the core: `if not
   marker_type_overdispersion:` runs the existing code path verbatim (untouched), so
   fixtures are guaranteed bit-identical. The new path is fully separate (chosen over
   making the two-rho path collapse to the one-rho path at runtime, which is too easy
   to perturb in the last ULPs and break fixtures).

3. **Sparse-class fallback -> QC warning, not stderr.** If either class has fewer than
   `MIN_CLASS_MARKERS` informative markers, its rho is not identifiable; fall back to
   the **shared-rho** path for the whole sample. Here ~348 hom / ~234 het, both ample.
   `MIN_CLASS_MARKERS = 30`. **The fallback must be visible in the structured output,
   not stderr** (the SRP runner and the paper pipeline shell out to `allomix monitor`
   as a subprocess, so stderr is lost). Record the fallback on the result and have the
   QC report surface it as a warning. See 5.3 / 5.8.

4. **Reported `rho`.** Confirmed: add optional fields `rho_hom: float | None` and
   `rho_het: float | None`, and set the headline `rho` to the **het-class rho** (the
   one that governs the floor and the low-fraction CI). `detection_limit` gets the
   per-class rhos (see 5.5).

5. **Single-donor first.** Ship single-donor only in the first commit; the issue's
   floor evidence is single-donor. Multi-donor (`_estimate_multi_donor_core`) is a
   later phase; only the class-partition helper is stubbed now (5.6).

## 5. Code snippets (illustrative, to be written)

### 5.1 Marker-class partition

```python
def _donor_het_mask(markers: list[InformativeMarker]) -> np.ndarray:
    """Boolean mask, True where the (single) donor genotype is heterozygous.

    Donor-het markers sit at background VAF ~0.5 (symmetric amplification
    scatter); donor-hom markers sit near 0/1 (one-sided host signal). The two
    classes carry different overdispersion, which is what the per-marker-type
    mode fits separately.
    """
    return np.fromiter(
        ((m.donor_gts[0][0] + m.donor_gts[0][1]) == 1 for m in markers),
        dtype=bool,
        count=len(markers),
    )
```

### 5.2 Two-rho profiled log-likelihood at fixed f

```python
# Module constant, near _RHO_MIN / _RHO_MAX.
MIN_CLASS_MARKERS = 30  # below this a class's rho is not identifiable -> shared-rho fallback


def _profile_rho_at_f(arr: _MarkerArrays, f: float, error_rate: float) -> tuple[float, float]:
    """Max LL over rho in [1, _RHO_MAX] at fixed f for one marker class.

    Returns (ll_max, rho_at_max). Mirrors the single-rho profile already used in
    the grid search and CI, so the per-class arithmetic is identical to today's.
    """
    p_alt = _p_alt_for_f(arr, f, error_rate)
    opt = minimize_scalar(
        lambda log_r: -_ll_from_p_alt(arr, p_alt, math.exp(log_r)),
        bounds=(math.log(1.0), math.log(_RHO_MAX)),
        method="bounded",
    )
    return -float(opt.fun), math.exp(float(opt.x))


def _two_rho_profile_ll(
    arr_hom: _MarkerArrays, arr_het: _MarkerArrays, f: float, error_rate: float
) -> float:
    """total_LL(f) = max_rho_hom LL(hom) + max_rho_het LL(het)."""
    ll_hom, _ = _profile_rho_at_f(arr_hom, f, error_rate)
    ll_het, _ = _profile_rho_at_f(arr_het, f, error_rate)
    return ll_hom + ll_het
```

### 5.3 New core branch

Add a parameter and branch at the top of `_estimate_single_donor_bb_core`:

```python
def _estimate_single_donor_bb_core(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 1001,
    calibration: PanelCalibration | None = None,
    marker_type_overdispersion: bool = False,
) -> ChimerismResult:
    cal = calibration or PanelCalibration()
    n_informative = len(markers)
    if n_informative == 0:
        return _empty_single_result(error_rate)  # existing empty-case block

    if marker_type_overdispersion:
        het_mask = _donor_het_mask(markers)
        n_het = int(het_mask.sum())
        n_hom = n_informative - n_het
        if n_het >= MIN_CLASS_MARKERS and n_hom >= MIN_CLASS_MARKERS:
            return _estimate_single_donor_two_rho(
                markers, het_mask, error_rate, grid_steps, cal
            )
        # Sparse class: rho not identifiable. Fall through to shared-rho, but
        # record it so QC can surface a warning (stderr is lost when the CLI is
        # driven as a subprocess by the runner / paper pipeline).
        fell_back = (
            f"per-marker-type overdispersion requested but a class is sparse "
            f"(hom={n_hom}, het={n_het}, min={MIN_CLASS_MARKERS}); used shared rho"
        )
    else:
        fell_back = None

    # ---- existing shared-rho code path, unchanged except the final result sets
    #      marker_type_overdispersion_fallback=fell_back ----
    arr = _precompute_marker_arrays(markers, cal)
    ...
    return ChimerismResult(..., marker_type_overdispersion_fallback=fell_back)
```

When the flag is off, `fell_back` is `None` and the result field stays `None`, so the
default path is still byte-identical (a `None` default field does not change existing
fixtures that never set it; confirm the fixture comparison ignores or matches it).

`_estimate_single_donor_two_rho` is a near-copy of the shared-rho body with the three
profile points doubled:

```python
def _estimate_single_donor_two_rho(
    markers, het_mask, error_rate, grid_steps, cal,
) -> ChimerismResult:
    hom_markers = [m for m, h in zip(markers, het_mask) if not h]
    het_markers = [m for m, h in zip(markers, het_mask) if h]
    arr_hom = _precompute_marker_arrays(hom_markers, cal)
    arr_het = _precompute_marker_arrays(het_markers, cal)

    # Step 1: grid over f, both rhos profiled out at each grid point.
    grid = np.linspace(0.0, 1.0, grid_steps)
    best_ll, best_f, best_rho_hom, best_rho_het = -math.inf, 0.0, 100.0, 100.0
    for f in grid:
        ll_hom, r_hom = _profile_rho_at_f(arr_hom, f, error_rate)
        ll_het, r_het = _profile_rho_at_f(arr_het, f, error_rate)
        ll = ll_hom + ll_het
        if ll > best_ll:
            best_ll, best_f, best_rho_hom, best_rho_het = ll, float(f), r_hom, r_het

    # Step 2: Nelder-Mead over (f, log_rho_hom, log_rho_het).
    def neg_ll_joint(x):
        f_val, lr_hom, lr_het = x
        if f_val < 0.0 or f_val > 1.0:
            return _INFEASIBLE_PENALTY
        r_hom, r_het = math.exp(lr_hom), math.exp(lr_het)
        if not (_RHO_MIN <= r_hom <= _RHO_MAX) or not (_RHO_MIN <= r_het <= _RHO_MAX):
            return _INFEASIBLE_PENALTY
        return -(_total_ll_vec(arr_hom, f_val, error_rate, r_hom)
                 + _total_ll_vec(arr_het, f_val, error_rate, r_het))

    opt = minimize(
        neg_ll_joint,
        x0=[best_f, math.log(best_rho_hom), math.log(best_rho_het)],
        method="Nelder-Mead",
        options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 5000},
    )
    f_mle = max(0.0, min(1.0, float(opt.x[0])))
    rho_hom_mle = math.exp(float(opt.x[1]))
    rho_het_mle = math.exp(float(opt.x[2]))

    # Step 3: profile-likelihood CI for f, both rhos profiled out at each f.
    threshold = chi2.ppf(CI_LEVEL, df=1)
    half = threshold / 2.0

    def profile_ll_f(f_val):
        return _two_rho_profile_ll(arr_hom, arr_het, f_val, error_rate)

    ll_max = profile_ll_f(f_mle)
    def ci_func(f_val):
        return ll_max - profile_ll_f(f_val) - half
    f_lo = 0.0 if (f_mle <= 0.0 or ci_func(0.0) <= 0.0) else brentq(ci_func, 0.0, f_mle, xtol=1e-5)
    f_hi = 1.0 if (f_mle >= 1.0 or ci_func(1.0) <= 0.0) else brentq(ci_func, f_mle, 1.0, xtol=1e-5)

    # Step 4: per-marker residuals (full set, at f_mle) -- unchanged helper.
    per_marker = _compute_per_marker_results(markers, f_mle, cal)
    n_markers_used = sum(1 for mr in per_marker if mr.included)

    # Step 5: detection limit under per-class rho.
    lob, lod = detection_limit(
        markers, error_rate, rho=rho_het_mle, calibration=cal,
        rho_hom=rho_hom_mle, rho_het=rho_het_mle,   # see 6
    )

    return ChimerismResult(
        donor_fraction=f_mle,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - f_mle,
        log_likelihood=ll_max,
        n_informative=len(markers),
        n_markers_used=n_markers_used,
        per_marker=per_marker,
        error_rate=error_rate,
        rho=rho_het_mle,           # headline rho = het class (governs the floor)
        rho_hom=rho_hom_mle,
        rho_het=rho_het_mle,
        lob_fraction=lob,
        lod_fraction=lod,
    )
```

Note Step 4 keeps `_compute_per_marker_results` over the **full** marker set at the
single `f_mle`, so per-marker output, residuals, and the robust-refit interaction are
unchanged in shape; only the `f` they are evaluated at moves.

### 5.4 `ChimerismResult` new fields

```python
@dataclass
class ChimerismResult:
    ...
    rho: float = float("inf")  # headline concentration; in two-rho mode = het-class rho
    # Per-marker-type overdispersion (issue #33). None unless that mode ran with
    # both classes above MIN_CLASS_MARKERS. rho_het governs the donor-het class
    # (background VAF ~0.5, the over-dispersed one); rho_hom the donor-hom class.
    rho_hom: float | None = None
    rho_het: float | None = None
    # Set to a human-readable reason when --marker-type-overdispersion was
    # requested but a class was too sparse to identify its rho, so the estimator
    # fell back to shared rho. None otherwise. Surfaced as a QC warning (5.8)
    # because stderr is lost when the CLI runs as a subprocess.
    marker_type_overdispersion_fallback: str | None = None
```

### 5.5 `detection_limit` / `fraction_se` per-class rho

Keep the scalar path byte-identical; add an optional per-class override:

```python
def fraction_se(
    markers, f_donor, error_rate=DEFAULT_ERROR_RATE, rho=float("inf"),
    calibration=None, rho_hom=None, rho_het=None,
) -> float:
    """... rho_hom/rho_het: when both given, each marker uses its class rho
    instead of the scalar `rho` (per-marker-type overdispersion, issue #33)."""
    ...
    for m in markers:
        ...
        if rho_hom is not None and rho_het is not None:
            m_rho = rho_het if (m.donor_gts[0][0] + m.donor_gts[0][1]) == 1 else rho_hom
        else:
            m_rho = rho
        overdispersion = 1.0 if math.isinf(m_rho) else 1.0 + (n - 1.0) / (m_rho + 1.0)
        ...
```

`detection_limit` just forwards `rho_hom`/`rho_het` to its two `fraction_se` calls.
With both `None` the arithmetic is identical to today.

### 5.6 Multi-donor stub (defer the wiring)

```python
def _donor_het_mask_multi(markers, n_donors) -> np.ndarray:
    """True where the combined donor background ALT balance is intermediate
    (not near 0 or 1), the multi-donor analogue of donor-het. Background ALT
    dose = mean over donors of donor_alt_dose/2; 'intermediate' = not within
    eps of {0, 1}. To be applied in _estimate_multi_donor_core in phase 2."""
    ...
```

### 5.7 Wrapper + CLI threading

- `estimate_single_donor_bb(..., marker_type_overdispersion: bool = False)` -> forwards
  to `core` via the existing `core(mk)` closure (`chimerism.py:1150`).
- `analyse_sample(..., marker_type_overdispersion: bool = False)` (`analysis.py:87`) ->
  forwards to `estimate_single_donor_bb` (`:154`) and, later, `estimate_multi_donor`.
- CLI `--marker-type-overdispersion` (store_true) in `cli.py`, threaded through
  `_analyse` (`cli.py:240`) the same way `robust` is (`:251`, `:280`).

### 5.8 QC warning for the sparse-class fallback

`qc.py` builds `QCReport.warnings: list[str]` (`qc.py:119`) from the result. Add one
line where the other result-driven warnings are appended (e.g. near the robust-drop
block at `qc.py:472-494`):

```python
# Per-marker-type overdispersion fell back to shared rho (sparse class). Inform
# the user in the structured output; stderr is lost when allomix runs as a
# subprocess. A plain warning, not a REVIEW promotion: the shared-rho estimate is
# the validated default, just not the mode that was asked for.
fb = getattr(result, "marker_type_overdispersion_fallback", None)
if fb:
    warnings.append(fb)
```

Keep it a warning (no status change). The shared-rho result it falls back to is the
validated default; the message only tells the user the requested mode did not engage.
This is the path the runner/paper pipeline actually sees, so it is how a silent
fallback becomes visible.

### 5.9 Fast grid estimator (`paper/scripts/fast_grid_estimator.py`)

The fast grid is the paper sweep's opt-in single-rho approximation
(`estimate_single_donor_bb_grid`). It profiles a single rho by taking the max over a
log-spaced rho-grid (`_ll_grid_over_rho` -> `(n_f, n_rho)`, then `.max(axis=1)`). The
two-rho version profiles each class's rho the same way and sums the per-class
profiled-over-rho curves before the f argmax. No new kernel: reuse `_p_alt_grid`,
`_ll_grid_over_rho`, and the core's `_donor_het_mask` (import it from
`allomix.chimerism`).

```python
def estimate_single_donor_bb_grid(
    markers, error_rate=DEFAULT_ERROR_RATE, calibration=None,
    n_f=201, n_rho=32, refine=True, marker_type_overdispersion=False,
) -> GridChimerismResult:
    cal = calibration or PanelCalibration()
    n_informative = len(markers)
    if n_informative == 0:
        return GridChimerismResult(...)  # unchanged empty case

    if marker_type_overdispersion:
        het_mask = _donor_het_mask(markers)
        n_het = int(het_mask.sum())
        n_hom = n_informative - n_het
        if n_het >= MIN_CLASS_MARKERS and n_hom >= MIN_CLASS_MARKERS:
            return _estimate_single_donor_bb_grid_two_rho(
                markers, het_mask, error_rate, cal, n_f, n_rho, refine
            )
        # sparse class -> fall through to the unchanged single-rho grid, recording
        # the fallback on the result (5.8 analogue for GridChimerismResult).

    # ---- existing single-rho grid path, unchanged ----
    ...
```

```python
def _estimate_single_donor_bb_grid_two_rho(
    markers, het_mask, error_rate, cal, n_f, n_rho, refine,
) -> GridChimerismResult:
    hom_markers = [m for m, h in zip(markers, het_mask) if not h]
    het_markers = [m for m, h in zip(markers, het_mask) if h]
    arr_hom = _precompute_marker_arrays(hom_markers, cal)
    arr_het = _precompute_marker_arrays(het_markers, cal)

    f_grid = np.linspace(0.0, 1.0, n_f)
    rho_grid = np.exp(np.linspace(math.log(_GRID_RHO_LO), math.log(_RHO_MAX), n_rho))

    # Profile rho out of each class on the grid, then sum the per-class profiles.
    p_hom = _p_alt_grid(arr_hom, f_grid, error_rate)            # (n_f, M_hom)
    p_het = _p_alt_grid(arr_het, f_grid, error_rate)            # (n_f, M_het)
    ll_hom = _ll_grid_over_rho(arr_hom, p_hom, rho_grid)        # (n_f, n_rho)
    ll_het = _ll_grid_over_rho(arr_het, p_het, rho_grid)        # (n_f, n_rho)
    prof_hom = ll_hom.max(axis=1)                               # (n_f,)
    prof_het = ll_het.max(axis=1)
    prof_total = prof_hom + prof_het                           # rho profiled, per f

    fi = int(np.argmax(prof_total))
    best_f = float(f_grid[fi])
    best_ll = float(prof_total[fi])
    best_rho_hom = float(rho_grid[int(np.argmax(ll_hom[fi]))])
    best_rho_het = float(rho_grid[int(np.argmax(ll_het[fi]))])

    def profile_two_rho_ll(f):  # exact per-class rho profile, mirrors the core
        return (_profile_rho_at_f_grid(arr_hom, f, error_rate)
                + _profile_rho_at_f_grid(arr_het, f, error_rate))

    if refine:
        step = 1.0 / (n_f - 1)
        lo, hi = max(0.0, best_f - 2 * step), min(1.0, best_f + 2 * step)
        if hi > lo:
            opt = minimize_scalar(lambda f: -profile_two_rho_ll(f),
                                  bounds=(lo, hi), method="bounded",
                                  options={"xatol": 1e-8})
            if -float(opt.fun) >= best_ll:
                best_f = max(0.0, min(1.0, float(opt.x)))
                best_ll = -float(opt.fun)
                best_rho_hom = _argmax_rho_at_f(arr_hom, best_f, error_rate, rho_grid)
                best_rho_het = _argmax_rho_at_f(arr_het, best_f, error_rate, rho_grid)

    # Coarse CI from the summed profiled-over-rho curve (same bracket as single-rho).
    half = float(chi2.ppf(CI_LEVEL, df=1)) / 2.0
    above = prof_total >= (best_ll - half)
    idx = np.nonzero(above)[0]
    f_lo, f_hi = (float(f_grid[idx[0]]), float(f_grid[idx[-1]])) if idx.size else (best_f, best_f)

    return GridChimerismResult(
        donor_fraction=best_f, donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - best_f, log_likelihood=best_ll,
        n_informative=n_informative, n_markers_used=n_informative,
        error_rate=error_rate, rho=best_rho_het,  # headline = het class
        rho_hom=best_rho_hom, rho_het=best_rho_het,   # new GridChimerismResult fields
    )
```

`_profile_rho_at_f_grid` is the 1-D bounded rho profile already inlined in the
single-rho `profile_ll_f` (lines `:245-253`); factor it out so both modes share it.
`_argmax_rho_at_f` is the same `minimize_scalar` returning the rho, mirroring
`:274-280`. Add `rho_hom`/`rho_het` (and optionally a fallback string) to
`GridChimerismResult` alongside `rho`.

**Validate the fast two-rho path against the exact two-rho estimator** in
`paper/scripts/validate_grid_estimator.py`, exactly as the single-rho grid is
validated against the exact single-rho MLE today (target: fraction agreement < 1e-3
over the LoD/ladder regime). Only after that passes is `fast_grid=1` trustworthy for
the new mode.

## 6. Tests

### 6.1 Unit tests (`tests/test_chimerism.py`, fast)

Reuse the existing `_make_marker` / `_make_markers_for_fraction` helpers. Add a
helper that builds a mixed hom+het marker set at a known low fraction with injected
symmetric het overdispersion (the floor's source):

```python
def _make_mixed_class_markers(f_host, n_hom, n_het, dp, het_overdisp_rho, seed):
    """Donor-hom markers: clean one-sided host signal at fraction f_host.
    Donor-het markers: background VAF 0.5 with beta-binomial scatter at
    het_overdisp_rho (the symmetric scatter that rectifies into a false +f)."""
```

Cases:

1. **`test_default_path_byte_identical`** — with `marker_type_overdispersion=False`,
   `estimate_single_donor_bb` returns exactly today's `donor_fraction`, `ci`,
   `log_likelihood`, `rho` on `_make_markers_overdispersed`. Assert equality (not
   `approx`) against a value captured from the current estimator. This is the
   fixture-safety guard.
2. **`test_partition_mask`** — `_donor_het_mask` is True only for donor `(0,1)` /
   `(1,0)` and False for `(0,0)` / `(1,1)`.
3. **`test_two_rho_removes_het_floor`** — on `_make_mixed_class_markers(f_host=0.0,
   het_overdisp_rho=70, ...)`, shared-rho mode returns a positive offset (~the floor)
   while `marker_type_overdispersion=True` returns ~0 (assert `< 0.1 * shared_estimate`
   and `< 0.0015`). The core regression test for the fix.
4. **`test_two_rho_recovers_signal_at_higher_fraction`** — at `f_host=0.05` the
   two-rho estimate is within tolerance of truth and its CI is **no wider** than the
   hom-only estimate (the precision-recovery claim vs a hom-only switch).
5. **`test_sparse_class_falls_back`** — with `n_het < MIN_CLASS_MARKERS`, the two-rho
   request returns the shared-rho result exactly (and the stderr note fires).
6. **`test_rho_fields_populated`** — two-rho result has non-None `rho_hom`/`rho_het`,
   `rho == rho_het`; shared-rho result has them `None`.
7. **`test_ci_monotone_and_contains_mle`** — `profile_ll_f` is maximized at `f_mle`
   and `f_lo <= f_mle <= f_hi`, `0 <= f_lo`, `f_hi <= 1` (CI well-formedness; the
   thing the issue's validation TODO flags as the real risk).
8. **`test_detection_limit_per_class_rho`** — `fraction_se` with `rho_hom=rho_het=R`
   equals `fraction_se` with scalar `rho=R` (consistency), and with `rho_het < rho_hom`
   the SE is larger than the all-`rho_hom` case (het markers down-weighted).
9. **`test_fallback_sets_qc_warning`** — request the mode with a sparse het class; the
   `ChimerismResult.marker_type_overdispersion_fallback` string is set and the built
   `QCReport.warnings` contains it (status unchanged). Guards the stderr-is-lost fix.

Fast grid (`tests/test_fast_grid_estimator.py`, alongside the existing single-rho
agreement tests):

10. **`test_grid_two_rho_matches_exact_two_rho`** — `estimate_single_donor_bb_grid(...,
    marker_type_overdispersion=True)` agrees with the exact two-rho MLE on
    `donor_fraction` to < 1e-3 on a panel-sized mixed-class set (the fast-grid
    accuracy contract, extended to the new mode).
11. **`test_grid_two_rho_default_unchanged`** — with the flag off, the grid result is
    identical to today's (the fast-grid byte-identical guard).

### 6.2 One-off scripts (`output/`, not committed, gitignored)

Standalone scripts run by hand to confirm the mechanism on real-shaped data before
the full sweep. These are the "have I got it right" checks, not pytest.

- **`scripts/oneoff_rho_vs_floor.py`** — on the committed synthetic ladder's donor-het
  markers, sweep a fixed rho over `{inf, 1000, 100, 71}` and print the recovered f at
  0.1% and 0.5% nominal. Reproduces the issue's table (binomial -> 0.000, fitted rho
  -> ~0.48). Confirms our class partition selects the same markers the investigation
  used.
- **`scripts/oneoff_two_rho_ladder.py`** — run `estimate_single_donor_bb` with the flag
  off vs on across the full synthetic ladder (10 mixes x 5 reps), print median host%
  per nominal level for both. Target reproduction of the issue's bottom row:
  `0.100 (+0.00)` at 0.1%, `0.450 (-0.05)` at 0.5%, and confirm flag-off reproduces
  `0.319 / 0.738`.
- **`scripts/oneoff_ci_coverage.py`** — Monte Carlo: simulate N>=200 replicate samples
  at each of f in {0, 0.001, 0.005, 0.02, 0.1, 0.5}, count how often the two-rho 95%
  CI contains truth. The whole point of rho is interval calibration, and only the
  point estimate was prototyped in the issue. Report coverage for flag-off vs flag-on.

### 6.3 Benchmarks (`scripts/bench_two_rho.py`)

The core is ~99% of paper build time (`chimerism.py:899-901`), and the two-rho path
roughly doubles the per-f-point rho profiling plus adds a simplex dimension. Measure
it so we know the cost of flipping the default and whether the fast-grid sweep needs
the new mode too.

```python
# Wall-clock per estimate, single donor, ~580 markers (panel-sized), grid_steps=1001:
#   - flag off (shared rho)         baseline
#   - flag on  (two rho, both ample)
#   - flag on  (sparse -> fallback) should equal baseline
# Report mean +/- sd over 20 calls each; assert flag-on is within ~3x of baseline
# (expected ~2x). Also time grid_steps=101 (quick) to size the LoD sweep impact.
```

Acceptance: if flag-on is <=~2.5x baseline at panel size we can consider it for the
default later; if it is much worse, the grid rho-profiling needs the vectorized
`_ll_grid_over_rho` treatment from the fast-grid estimator (note for a follow-up).

## 7. Validation sequence (gated)

The fastgrid two-rho path (5.9) is **not** built until the exact estimator is proven
on the ladder. fastgrid is an approximation, and its existing trust (0.0011 pp on the
LoD-summary cells) does not cover this regime: the het-class rho that controls the
floor sits near the grid's coarsest (lower) rho bound, and the effect we are measuring
is a ~0.22 pp floor at sub-0.5%, so a grid discretization artifact there could be a
real fraction of the signal in either direction. Building fastgrid before the exact
check risks both wasted rework (if the approach needs tuning) and validating the method
on a grid artifact rather than the method. So the order is gated, not parallel.

**Stage 0 -- env** (paper extras, per CLAUDE.md):
```bash
uv venv --python 3.13 && source .venv/bin/activate
uv pip install -e ".[dev,scripts,paper]"
```

**Stage 1 -- unit tests + one-offs (exact estimator only).** Section 6.1 unit tests
pass (esp. byte-identical default), and the section 6.2 one-offs reproduce the issue's
rho-vs-floor table and the flag-off ladder numbers (0.319 / 0.738). No fastgrid yet.

**Stage 2 -- GO/NO-GO: exact two-rho on the synthetic ladder.** This is the decisive
check and it uses **no fastgrid**. The SRP434573 ladder is small (10 mixes x 5 reps),
so it runs under the exact estimator via the CLI. `run_srp434573_allomix.py` shells out
to `allomix monitor` (`ALLOMIX = [..., "monitor"]`), so wiring
`--marker-type-overdispersion` into the CLI and adding it to that invocation is enough;
regenerate `output/srp434573_synthetic.tsv`.
  - **GO criteria:** median host% drops to ~`0.100` at 0.1% nominal and ~`0.450` at
    0.5% (the issue's target bottom row), flag-off still reproduces `0.319 / 0.738`,
    host-dose slope stays 0.92-0.99, and CI coverage (one-off `oneoff_ci_coverage.py`)
    is acceptable.
  - **If NO-GO** (floor not removed, slope degraded, or coverage broken): stop and
    re-tune the exact estimator. **Do not build fastgrid.** fastgrid would only be an
    approximation of a method that does not yet work.

**Stage 3 -- full-range exact check: extend the ladder past 0.5% (still no fastgrid).**
The synthetic ladder stops at 0.5%; the issue asks to confirm the fix holds across
0-100%. Regenerate a few higher-fraction semi-synthetic mixes (1%, 5%, 20%, 50%) and,
with the exact estimator, confirm two-rho tracks truth there and the het markers regain
weight (CI no wider than shared-rho at >=5%). This is a correctness gate, kept before
the expensive sweep deliberately: if precision recovery fails at higher fractions the
method loses its main advantage over a hom-only switch, and we want that before
committing the full build. Exact estimator, same CLI path as Stage 2.

**Stage 4 -- build + validate fastgrid two-rho (only after Stages 2-3 hold).** Implement
5.9, then validate it in `validate_grid_estimator.py` against the now-trusted exact
two-rho MLE. The agreement target is fraction < 1e-3, and the comparison **must include
the sub-0.5% regime**, not only the LoD-summary cells, because that is where the grid
is coarsest and the effect is smallest. If it cannot meet < 1e-3 sub-0.5%, raise
`n_rho` / refine bounds there or keep the LoD sweep on the exact estimator for the new
mode.

**Stage 5 -- full paper sweep with `--config fast_grid=1`** (tractability for the LoD
sweep, ~99% of build time):
```bash
snakemake -s paper/Snakefile --cores $(nproc) --config fast_grid=1
```
The runner's `--marker-type-overdispersion` toggle and the sweep's `fast_grid=1` are
independent: the SRP ladder uses the exact estimator via the CLI; the LoD sweep uses
the fast grid via `run_lod_validation.py`.

Report, per stage: flag-off vs flag-on median host% table, CI coverage, slope, and the
benchmark numbers.

## 8. Rollout (per davmlaw on #33)

1. Land opt-in (`--marker-type-overdispersion`, default off). Fixtures byte-identical.
2. Run + validate (sections 6-7): point estimate, **CI coverage**, full-range, timing.
3. If pure win: make it the default, wire into the paper Snakemake (the fast-grid
   two-rho path from 5.9 is already in place), and update **methods** and **results**
   (including the `paper/results.md` "semi-synthetic sub-0.5% ladder" paragraph and the
   quick MLE +0.22 edit that attributed the floor to overdispersion).
4. Multi-donor (`_estimate_multi_donor_core`) is a separate later phase.
5. Do **not** auto-close #33; reference with a plain `#33`.

## 9. Open risks

- **CI coverage is the real unknown.** Two classes each fit their own rho on fewer
  markers, so the profile-likelihood interval could be miscalibrated even when the
  point estimate is unbiased. Section 6.2 `oneoff_ci_coverage.py` and unit test 7 must
  clear before considering default.
- **Identifiability at the boundary.** When `f` is near 0 the het class is pure noise
  and its rho hits the lower bound; the profile there must stay finite (the
  `_RHO_MIN`/`_RHO_MAX` clamps and `_INFEASIBLE_PENALTY` already guard this, but verify
  no `brentq` sign error like the one the single-rho CI guards against at `:957`).
- **Robust-refit interaction.** `_robust_refit` calls `core_fn` on marker subsets; a
  trim can drop a class below `MIN_CLASS_MARKERS` mid-iteration and switch that refit
  to shared-rho. Acceptable (it is the documented fallback), and the
  `marker_type_overdispersion_fallback` string propagates through `_robust_refit`'s
  `replace(result, ...)` so it still reaches QC. Worth a test that the mode +
  `--robust auto` together do not crash and converge.
- **Multi-donor left for phase 2**; single-donor only in the first commit.
