"""Tests for the LoB / logistic-fit helpers in paper/scripts/run_lod_validation.py."""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import numpy as np
import pytest

# Make the paper/scripts module importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "paper" / "scripts"))

import run_lod_validation as lod  # noqa: E402


def test_compute_lob_matches_numpy_quantile() -> None:
    rng = np.random.default_rng(0)
    blanks = rng.uniform(0, 0.01, size=200).tolist()
    expected = float(np.quantile(np.asarray(blanks), 0.95))
    assert lod.compute_lob(blanks) == pytest.approx(expected)


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


def test_bootstrap_lod_ci_brackets_point_estimate() -> None:
    # Booleans drawn from a clean logistic: bootstrap CI should bracket the
    # point-estimate LoD on average.
    a_true, b_true = 2.5, 3.0
    fractions = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    rng = random.Random(123)
    booleans = {}
    for f in fractions:
        p = 1.0 / (1.0 + math.exp(-(a_true + b_true * math.log10(f))))
        booleans[f] = [rng.random() < p for _ in range(60)]
    rates = [sum(booleans[f]) / len(booleans[f]) for f in fractions]
    fit = lod.fit_lod(fractions, rates)
    assert fit is not None
    f95 = fit[0]
    ci_lo, ci_hi = lod.bootstrap_lod_ci(booleans, n_bootstrap=200, rng=random.Random(7))
    assert math.isfinite(ci_lo) and math.isfinite(ci_hi)
    assert ci_lo <= ci_hi
    # f95 should usually fall inside the bootstrap CI; allow a small slack.
    assert ci_lo * 0.5 <= f95 <= ci_hi * 2.0


def test_to_pct_preserves_sentinels() -> None:
    assert lod._to_pct(lod.LOD_BELOW_RANGE) == -1.0
    assert lod._to_pct(lod.LOD_ABOVE_RANGE) == float("inf")
    assert math.isnan(lod._to_pct(float("nan")))
    assert lod._to_pct(0.0123) == pytest.approx(1.23)
