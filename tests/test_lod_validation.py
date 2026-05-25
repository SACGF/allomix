"""Tests for the LoB / logistic-fit helpers in paper/scripts/run_lod_validation.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the paper/scripts module importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "paper" / "scripts"))

import run_lod_validation as lod  # noqa: E402


def test_compute_lob_parametric() -> None:
    rng = np.random.default_rng(0)
    blanks = rng.normal(0.0, 0.001, size=200).tolist()
    arr = np.asarray(blanks)
    expected = float(arr.mean() + 1.645 * arr.std(ddof=1))
    assert lod.compute_lob(blanks) == pytest.approx(expected)


def test_compute_lob_approaches_95th_percentile_for_normal() -> None:
    # For Gaussian blanks, parametric LoB and empirical 95th percentile converge
    # in distribution. At large n they should agree to within a few percent.
    rng = np.random.default_rng(7)
    blanks = rng.normal(0.0, 0.001, size=10000).tolist()
    parametric = lod.compute_lob(blanks)
    empirical = float(np.quantile(np.asarray(blanks), 0.95))
    assert parametric == pytest.approx(empirical, rel=0.05)


def test_compute_lob_empty_is_nan() -> None:
    assert math.isnan(lod.compute_lob([]))


def test_detection_rate_basic() -> None:
    assert lod.detection_rate([0.0, 0.001, 0.005, 0.010], lob=0.002) == pytest.approx(0.5)
    assert lod.detection_rate([], lob=0.0) == 0.0


def test_logistic_fit_recovers_known_params() -> None:
    # Construct deterministic data from a known logistic curve, then check we
    # recover the f95 to within floating-point precision.
    a_true, b_true = 2.5, 3.0  # P(det) = sigmoid(2.5 + 3 log10 f)
    fractions = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    rates = [1.0 / (1.0 + math.exp(-(a_true + b_true * math.log10(f)))) for f in fractions]
    fit = lod.fit_lod(fractions, rates)
    assert fit is not None
    f95, a_hat, b_hat = fit
    expected_log10_f95 = (lod.LOGIT_95 - a_true) / b_true
    assert math.log10(f95) == pytest.approx(expected_log10_f95, rel=1e-3)
    assert a_hat == pytest.approx(a_true, rel=1e-3)
    assert b_hat == pytest.approx(b_true, rel=1e-3)


def test_logistic_fit_handles_too_few_points() -> None:
    assert lod.fit_lod([0.01], [0.5]) is None


def test_logistic_fit_handles_degenerate() -> None:
    # All rates equal: curve_fit may converge with b ~= 0 — guard rejects it.
    fit = lod.fit_lod([0.001, 0.002, 0.005, 0.01], [0.5, 0.5, 0.5, 0.5])
    assert fit is None or not math.isfinite(fit[0]) or fit[0] > 0


def test_interp_lod_brackets_target() -> None:
    # det rate 0.95 lies between (0.005, 0.5) and (0.01, 1.0) -> log10-interpolate
    fractions = [0.001, 0.005, 0.01, 0.02]
    rates = [0.0, 0.5, 1.0, 1.0]
    f95 = lod._interp_lod(fractions, rates, target=0.95)
    assert f95 is not None
    assert 0.005 < f95 < 0.01
    expected = 10 ** (math.log10(0.005) + 0.9 * (math.log10(0.01) - math.log10(0.005)))
    assert f95 == pytest.approx(expected)


def test_interp_lod_returns_none_when_never_crossing() -> None:
    # Detection never reaches target -> bracketing impossible, callers handle
    # via the LOD_ABOVE_RANGE sentinel.
    assert lod._interp_lod([0.001, 0.01], [0.0, 0.5]) is None
    # Always at or above target -> bracketing also impossible (below smallest
    # tested fraction). Callers use LOD_BELOW_RANGE.
    assert lod._interp_lod([0.001, 0.01], [1.0, 1.0], target=0.95) is None


def test_fit_lod_falls_back_to_interp_on_step_data() -> None:
    # Step-like detection 0->1 across one fraction -- logistic slope is
    # unidentifiable. Interp fallback should still yield an LoD.
    fractions = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    rates = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    fit = lod.fit_lod(fractions, rates)
    assert fit is not None
    f95 = fit[0]
    assert 0.005 < f95 < 0.01


def test_derive_seed_is_stable_across_invocations() -> None:
    # Regression: Python's hash() is randomised per-process for str (PEP 456),
    # so an older `hash(repr(parts))`-based derive_seed silently produced
    # different "deterministic" seeds in each run of this sweep. Pin a known
    # SHA-256-derived value so future refactors can't reintroduce the bug.
    assert lod.derive_seed("gt", "unrelated", 0, 42) == 3162746855


def test_derive_seed_distinct_inputs_distinct_outputs() -> None:
    a = lod.derive_seed("gt", "unrelated", 0, 42)
    b = lod.derive_seed("gt", "unrelated", 1, 42)
    c = lod.derive_seed("bias", "unrelated", 0, 42)
    assert len({a, b, c}) == 3


def test_fit_lod_rejects_negative_slope() -> None:
    # Ultra-easy corner: ~50% detection at smallest fraction, 100% above.
    # curve_fit can converge to a negative-slope solution that algebraically
    # inverts to f95 > 1 (the regression that produced LoD = 156% on
    # unrelated/2000x/400 markers). Negative slope must be rejected so the
    # interp fallback supplies the real LoD.
    fractions = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    rates = [0.467, 1.0, 1.0, 1.0, 1.0, 1.0]
    fit = lod.fit_lod(fractions, rates)
    assert fit is not None
    f95 = fit[0]
    assert 0 < f95 < 0.01, f"expected sub-1% LoD, got {f95}"


def test_to_pct_preserves_sentinels() -> None:
    assert lod._to_pct(lod.LOD_BELOW_RANGE) == -1.0
    assert lod._to_pct(lod.LOD_ABOVE_RANGE) == float("inf")
    assert math.isnan(lod._to_pct(float("nan")))
    assert lod._to_pct(0.0123) == pytest.approx(1.23)


def test_lod_for_aggregation_maps_sentinels() -> None:
    # all_detected -> probed floor, none_detected -> probed ceiling, NaN dropped.
    assert lod._lod_for_aggregation(lod.LOD_BELOW_RANGE) == lod.MIN_POS_FRACTION
    assert lod._lod_for_aggregation(lod.LOD_ABOVE_RANGE) == lod.MAX_FRACTION
    assert lod._lod_for_aggregation(float("nan")) is None
    assert lod._lod_for_aggregation(0.003) == 0.003


def test_summarise_cell_median_and_band() -> None:
    # Curve point = median across pairs; band = 10th-90th percentile across pairs.
    pairs = [
        {"lod": 0.001, "lob": 0.0005, "mean_n_informative": 40, "note": ""},
        {"lod": 0.002, "lob": 0.0007, "mean_n_informative": 42, "note": ""},
        {"lod": 0.004, "lob": 0.0009, "mean_n_informative": 41, "note": ""},
    ]
    s = lod.summarise_cell("sibling", 1000, 100, pairs, n_seq_reps=30)
    assert s["lod_pct"] == pytest.approx(0.2)  # median 0.002 -> 0.2%
    assert s["n_pairs"] == 3 and s["n_pairs_used"] == 3 and s["n_pairs_dropped"] == 0
    assert s["lod_pct_ci_lo"] <= s["lod_pct"] <= s["lod_pct_ci_hi"]


def test_summarise_cell_drops_failed_fits() -> None:
    pairs = [
        {"lod": 0.002, "lob": 0.0007, "mean_n_informative": 42, "note": ""},
        {"lod": float("nan"), "lob": float("nan"), "mean_n_informative": 5,
         "note": "fit_failed"},
    ]
    s = lod.summarise_cell("sibling", 1000, 50, pairs, n_seq_reps=30)
    assert s["n_pairs"] == 2
    assert s["n_pairs_used"] == 1
    assert s["n_pairs_dropped"] == 1
    assert "dropped" in s["note"]


def test_summarise_cell_clamps_all_detected_pair() -> None:
    # A pair that detects everything (LoD below probed range) is clamped to the
    # floor before the median, so it pulls the cell LoD down but stays finite.
    pairs = [
        {"lod": lod.LOD_BELOW_RANGE, "lob": 0.0002, "mean_n_informative": 120,
         "note": "all_detected"},
        {"lod": 0.002, "lob": 0.0007, "mean_n_informative": 122, "note": ""},
    ]
    s = lod.summarise_cell("unrelated", 2000, 200, pairs, n_seq_reps=30)
    # median of [0.001 (clamped floor), 0.002] = 0.0015 -> 0.15%
    assert s["lod_pct"] == pytest.approx(0.15)
    assert s["n_pairs_used"] == 2
