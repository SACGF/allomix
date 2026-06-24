# Code review refactoring plan (June 24)

Agreed actions from the June 24 review of `dave_code_comments.md`. Two parts: a
structural reorganisation of `chimerism.py`, and a set of small local changes.

The structural part is purely organisational. `chimerism.py` is ~1510 lines and
mixes three concerns (result data types, the likelihood/weight model, and the
MLE estimators). Splitting them does not change behaviour or performance, it
just makes each file single-purpose and easier to read. The estimators stay
byte-identical, which matters because the single-donor output is fixture-locked
to ~0.1% at the limit of detection.

## Part 1: Split `chimerism.py` into three modules

Target layout:

- `results.py` (new): the output data types.
  - `MarkerResult`, `ChimerismResult`, `MultiDonorResult`.
  - These are the public data surface that `report.py`, `qc.py`, and
    `analysis.py` consume. None of those call an estimator, they just read
    result objects, so the data types belong on their own.
  - Imports `HostPresenceResult` (detect), `ContaminationResult`
    (contamination), `RelatednessResult` / `AdmixConsistencyResult`
    (relatedness), `RunUnitInfo` (runmeta). None of those import chimerism, so
    there is no cycle.

- `likelihood.py` (new): the likelihood and weight model, no optimisation.
  - `PanelCalibration` (input/config type; its `bias_for` / `error_for` are
    called by the likelihood functions, so it lives with them).
  - `W_EPS`, `apply_bias`, `inject_bias`.
  - `expected_weight`, `expected_weight_multi`.
  - `alt_read_probability`, `log_likelihood_marker_bb`.
  - `total_log_likelihood_bb`, `total_log_likelihood_multi_bb`.
  - The vectorised core: `_MarkerArrays`, `_precompute_marker_arrays`,
    `_total_ll_vec`, `_p_alt_for_f`, `_ll_from_p_alt`.
  - Imports only `constants`, `error_rates` (`MarkerErrorRates`), and
    `genotype` (`InformativeMarker`, `MarkerKey`). No cycle.

- `chimerism.py` (stays): the MLE / optimisation layer.
  - Robust-refit constants and `_robust_refit`, `_robust_trigger`,
    `_marker_key`.
  - `_compute_per_marker_results`, `_per_marker_results_multi`.
  - `fraction_se`, `detection_limit` (Fisher information / LoB / LoD).
  - `_estimate_single_donor_bb_core`, `estimate_single_donor_bb`.
  - `_estimate_multi_donor_core`, `estimate_multi_donor`,
    `_profile_likelihood_cis_multi`.

### Backward compatibility: re-export

~30 call sites across `src/`, `tests/`, `scripts/`, and `paper/scripts/` import
result types and functions from `allomix.chimerism`. To avoid touching all of
them, `chimerism.py` re-exports the moved names at the top:

```python
from allomix.results import ChimerismResult, MarkerResult, MultiDonorResult
from allomix.likelihood import (
    PanelCalibration,
    apply_bias,
    expected_weight,
    expected_weight_multi,
    inject_bias,
    total_log_likelihood_bb,
    total_log_likelihood_multi_bb,
    # ... and any other likelihood names current call sites import from chimerism
)
# The estimators (estimate_single_donor_bb, estimate_multi_donor) stay defined
# in chimerism.py, so they need no re-export.
```

Every existing `from allomix.chimerism import X` keeps working. This turns a
30-file edit into a 3-file change with no behaviour change.

### What this does not buy

It does not make `report.py` / `qc.py` cheaper to import. The result types
aggregate `HostPresenceResult` etc., and those modules pull in scipy regardless,
so `results.py` still transitively imports scipy. The benefit is readability:
three single-purpose files instead of one 1510-line file. That is the point.

## Part 2: Minor local changes

### 2.1 `total_log_likelihood_multi_bb` error-kwargs (now in `likelihood.py`)

Collapse the two `if entry is not None` checks into one, built as a 3-line
conditional with an indented body:

```python
entry = cal.error_for(m)
err_kwargs = {}
if entry is not None:
    err_kwargs = {"e_refalt": entry.e_refalt, "e_altref": entry.e_altref}
...
ll += log_likelihood_marker_bb(
    m.admix_ad_ref, m.admix_ad_alt, w, error_rate=error_rate, rho=rho, **err_kwargs
)
```

This is the multi-donor path, which is not the vectorised hot loop, so building
the per-marker dict is fine.

### 2.2 Disambiguate the two `2`s, and merge the weight helpers

There are two unrelated `2`s in this code:

- **Ploidy 2**: the `/2` and `2 - (gt[0] + gt[1])` dose math. "2 alleles per
  diploid genotype." Same everywhere.
- **Donor-count 2**: `n_donors=2`, the triangular `(f1, f2)` grid,
  `other_idx = 1 - donor_idx`. The multi-donor cap.

Actions:

- Add a `PLOIDY = 2` constant and use it for the ploidy `2`s, so they read as
  ploidy and are visibly distinct from the donor-count cap.
- Merge `expected_weight` into `expected_weight_multi`: they are algebraically
  identical (`expected_weight(h, d, f)` == `expected_weight_multi(h, [d], [f])`,
  bias branch included). `expected_weight` is not in the vectorised hot loop
  (that path is `_p_alt_for_f`), so make it a one-line delegator. Removes the
  duplicated dose math.
- In the multi-donor estimator, replace the bare `1 - donor_idx` and similar
  magic with a clearer construct plus a comment. The existing
  `if n_donors > 2: raise` already documents the cap.

Do **not** merge the two MLE estimators. They diverge for real reasons: the
single-donor path is vectorised, is ~99% of the paper build, and is
fixture-locked; the multi-donor path uses the per-marker Python loop and a
different optimiser (triangular grid, per-donor profile CIs). A `num_donors`
merge would either slow the validated hot path or keep the single-donor
specialisation anyway, and calling the multi path with `n_donors=1` would not
bit-match the fixtures.

### 2.3 Hoist the grid-loop bounds

In `_estimate_single_donor_bb_core`, the rho-profiling bounds
`(math.log(1.0), math.log(10000.0))` are recomputed on every grid point. Lift
them to a local computed once before the loop (a per-call local is fine, it does
not need to be a module constant; it just must be out of the hot loop).

While there: the grid loop uses an upper rho bound of `10000.0` but the
profile-CI search and Nelder-Mead use `_RHO_MAX = 50000.0`. Confirm whether that
asymmetry is intentional and add a one-line comment either way. Leave the values
unchanged for now (changing them would move fixtures).

### 2.4 report.py rounding constant

Dropped. The two conventions (`6` for raw fractions/VAFs, `4` for percentages,
`1` for depth) do not collapse to a single constant cleanly, and the literals
are clear enough in context.

### 2.5 Separate the paper tests from the main suite

`tests/test_semisynthetic_srp434573.py` and `tests/test_lod_validation.py` test
paper-runner glue in `paper/scripts/`, not the `src/allomix/` tool. We plan to
add CI when we go public and want the core gate fast and panel-focused.

- Move both files to `paper/tests/`.
- `testpaths = ["tests"]` already means the main `pytest` run stops collecting
  them, so no config change is needed for the core suite.
- Add a paper-test step to the Snakefile that runs `pytest paper/tests` before
  the paper build (either as a dependency of the `paper` target or a standalone
  `paper-test` target), so paper-glue regressions are still caught at build time.
- Do not delete them. They guard against shipping a wrong paper (e.g. a
  fraction-label parser regression mislabelling every LoD point).

## Sequencing and validation

1. Part 1 first (move + re-export), as three separate steps: `results.py`,
   then `likelihood.py`, then trim `chimerism.py` to re-exports. Run the
   relevant tests after each (`pytest tests/test_chimerism.py tests/test_multidonor.py tests/test_analysis.py -q`).
2. Then the 2.x local changes, each followed by the same targeted tests.
3. Confirm the single-donor fixtures are byte-identical (no estimator logic was
   touched). 2.2's `expected_weight` delegation is exact, so per-marker outputs
   should not move.
4. `ruff check` and `ruff format` on the new and edited files.
5. Commit to `main` per the repo convention.
</content>
</invoke>
