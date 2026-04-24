# Plan: Per-Marker Likelihood Context Refactor (Step 12)

Status: proposed, not implemented.

## TL;DR

Steps 14, 15, 16, and 17 each independently propose adding a new optional kwarg to `total_log_likelihood_bb`, `total_log_likelihood_multi_bb`, `estimate_single_donor_bb`, `estimate_multi_donor`, and `_profile_likelihood_cis_multi`, then threading it through every nested closure (`neg_ll_joint`, `profile_ll_f`, the per-donor inner optimisers, etc.). Each kwarg is a per-marker quantity (an empirical error rate, a dropout weight, a GQ weight, a base-quality-derived effective error). By the time all four land, those signatures will have grown 4 new optional parameters each, plus the existing `marker_biases`, and every call site has to forward all of them.

This step does that plumbing once, up front. The four downstream steps then mutate one field on a precomputed per-marker context list rather than threading another kwarg through every closure.

This is a pure refactor: no behavioural change, no new CLI flags, no new files for the user. The benefit is purely ergonomic for the next four steps.

## Goals

1. One signature change to the aggregator: `total_log_likelihood_bb(markers, f_donor, ctx, rho)` where `ctx` is a precomputed list aligned 1-to-1 with `markers`.
2. One precomputation pass per estimator call (in `estimate_single_donor_bb` / `estimate_multi_donor`) that builds the context list from whatever inputs are configured (today: `error_rate` + optional `marker_biases`).
3. The `_profile_likelihood_cis_multi` helper takes the same `ctx` and forwards it; closures capture it by reference.
4. `qc._compute_gof_pval` reads the same per-marker fields it needs (today: `error_rate`; later: per-site error and dropout) from `MarkerResult` rather than from a scalar argument.
5. Zero behavioural change, verified by the existing 261-test suite passing unchanged.

## Design

### `PerMarkerContext` dataclass

New small dataclass in `chimerism.py` (or a new `_likelihood_context.py` if it grows):

```python
@dataclass
class PerMarkerContext:
    """Per-marker quantities consumed by the likelihood, precomputed once
    per estimator call.

    Today this carries only the bias offset and the (scalar) error rate,
    which keeps the refactor a no-op. Future steps mutate fields on this
    object rather than adding kwargs to every aggregator and closure:

    - Step 14 (empirical error rates): add ``e_refalt`` / ``e_altref``.
    - Step 15 (dropout weighting): add ``dropout_weight``.
    - Step 16 (GQ weighting): add ``gq_weight``.
    - Step 17 (BQ-aware likelihood): replace scalar ``error_rate`` with
      a per-marker effective error.
    """

    bias: float = 0.0
    error_rate: float = 0.01
    # Future fields: e_refalt, e_altref, dropout_weight, gq_weight, ...
```

### Aggregator signature

```python
def total_log_likelihood_bb(
    markers: list[InformativeMarker],
    f_donor: float,
    ctx: list[PerMarkerContext],
    rho: float = 100.0,
) -> float:
    if len(ctx) != len(markers):
        raise ValueError(...)
    ll = 0.0
    for m, c in zip(markers, ctx):
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=c.bias)
        ll += log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w,
            error_rate=c.error_rate, rho=rho,
        )
    return ll
```

`total_log_likelihood_multi_bb` gets the same shape with the multi-donor `expected_weight_multi`.

### Estimator changes

`estimate_single_donor_bb` and `estimate_multi_donor` keep their public kwargs (`error_rate`, `marker_biases`) for backward compatibility. They build the `ctx` list once at the top:

```python
ctx = [
    PerMarkerContext(
        bias=marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
            if marker_biases is not None else 0.0,
        error_rate=error_rate,
    )
    for m in markers
]
```

Every `total_log_likelihood_bb(...)` call site inside the function (grid search, Nelder-Mead, profile-likelihood scan) passes `ctx` instead of `error_rate, marker_biases`. The closures (`neg_ll_joint`, `profile_ll_f`) capture `ctx` from the enclosing scope, same as today.

### `_profile_likelihood_cis_multi`

Add `ctx: list[PerMarkerContext]` to its signature; forward into the two inner optimisers.

### `_compute_per_marker_results` and `_per_marker_results_multi`

These currently take `error_rate` as a scalar. Change to take the same `ctx` so per-marker `error_rate` can be surfaced on `MarkerResult` for the QC step. Today `MarkerResult.effective_error_rate = ctx[i].error_rate` — identical to today's behaviour. When Step 17 lands the field starts varying per marker; `qc._compute_gof_pval` already plumbs `error_rate` through so this is a one-line change there.

### `qc.py`

`_compute_gof_pval` already accepts a scalar `error_rate`. Refactor it to read `m.effective_error_rate` from each `MarkerResult` (with a fallback to the scalar `error_rate` when absent for backward compat with old fixtures). No call-site change in `assess_quality`.

## Scope of the change

- `src/allomix/chimerism.py`: introduce `PerMarkerContext`; replace the `error_rate, marker_biases` kwargs on the two aggregators with a single `ctx` list; precompute `ctx` once at the top of each estimator; thread `ctx` through every closure and through `_profile_likelihood_cis_multi`; add `effective_error_rate: float | None = None` to `MarkerResult`.
- `src/allomix/qc.py`: `_compute_gof_pval` reads per-marker `effective_error_rate` if set, else falls back to the scalar arg.
- Tests: the existing 261 tests must pass unchanged. Any test that calls `total_log_likelihood_bb(...)` directly with `error_rate=...` or `marker_biases=...` will need a one-line update to pass a `ctx` list (or use a small `_make_ctx_for(markers, error_rate, biases)` helper added to the test module). Spot-check before committing — `grep -rn 'total_log_likelihood_bb\|total_log_likelihood_multi_bb' tests/` to find callers.

No CLI changes, no new files (other than possibly `_likelihood_context.py` if the dataclass grows during follow-up steps), no behavioural changes.

## Verification plan

1. `pytest -x -q` — full suite must stay green. This is the primary success criterion: a refactor that changes test results is a refactor that has slipped behaviour somewhere.
2. Numerical regression spot-check: run `allomix monitor` on `tests/test_data/multidonor/sample_h_d1_d2_f1_20_f2_10.vcf` before and after the refactor; donor fractions and CIs identical to ~12 decimals.
3. Real-data spot-check: re-run the April-24 validation batch (`output/validation_run_new_bias2/`) and diff `batch.tsv` against the pre-refactor copy. Bit-identical apart from any floating-point reordering inside the precomputed `ctx` list (none expected).

## File-by-file checklist

- [ ] `src/allomix/chimerism.py`: `PerMarkerContext` dataclass; `total_log_likelihood_bb` and `total_log_likelihood_multi_bb` signature change; `estimate_single_donor_bb` and `estimate_multi_donor` precompute `ctx` and thread it through; `_profile_likelihood_cis_multi` takes and forwards `ctx`; `_compute_per_marker_results` and `_per_marker_results_multi` populate `MarkerResult.effective_error_rate`; `MarkerResult` gains `effective_error_rate: float | None = None`.
- [ ] `src/allomix/qc.py`: `_compute_gof_pval` reads per-marker `effective_error_rate` with scalar fallback.
- [ ] `tests/test_chimerism.py`, `tests/test_multidonor.py`: any direct callers of the aggregators get a `_make_ctx_for` helper or inline `ctx` list. Existing assertions unchanged.
- [ ] Numerical regression spot-check on the multi-donor fixture and on the April-24 batch.

## Out of scope

- Any new behaviour. The point is to land the plumbing so Steps 14–17 are each a one-field addition to `PerMarkerContext` plus the math change in `log_likelihood_marker_bb` (or a multiplier in the aggregator), not a five-call-site signature growth.
- Renaming the public estimator kwargs (`error_rate`, `marker_biases`). Keep them; they are how `cli.py` and `scripts/run_xls_batch.py` already drive things, and there is no value in churning the surface area for a pure refactor.
