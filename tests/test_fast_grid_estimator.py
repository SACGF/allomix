"""Tests for the paper-only fast grid single-donor estimator.

The estimator lives in ``paper/scripts/fast_grid_estimator.py`` (outside the
``allomix`` package, since it is a paper-build speed optimisation, not part of the
clinical tool). These tests pin it to the exact estimator in ``allomix.chimerism``
so the LoD sweep stays trustworthy.
"""

import random
import sys
from pathlib import Path

import numpy as np
import pytest

from allomix.chimerism import estimate_single_donor_bb
from allomix.error_rates import MarkerErrorRates
from allomix.likelihood import (
    PanelCalibration,
    _p_alt_for_f,
    _precompute_marker_arrays,
)
from allomix.simulate import generate_marker_biases_realistic

# The fast grid estimator is a paper script, not an installed package module.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "paper" / "scripts"))

from fast_grid_estimator import (  # noqa: E402
    _p_alt_grid,
    estimate_single_donor_bb_grid,
)

# Marker-building helpers shared with the core chimerism tests.
from tests.test_chimerism import (  # noqa: E402
    _make_markers_for_fraction,
    _make_markers_overdispersed,
    _make_mixed_class_markers,
)


class TestGridEstimator:
    """The fast vectorized grid estimator must match the exact MLE fraction."""

    @pytest.mark.parametrize("f_true", [0.0, 0.005, 0.01, 0.05, 0.1, 0.5])
    @pytest.mark.parametrize("dp", [250, 1000])
    def test_grid_matches_exact_fraction(self, f_true: float, dp: int) -> None:
        markers = _make_markers_for_fraction(f_true, n_markers=100, dp=dp, seed=7)
        exact = estimate_single_donor_bb(markers, error_rate=0.01, grid_steps=201)
        grid = estimate_single_donor_bb_grid(markers, error_rate=0.01)
        # Match to < 0.01 percentage points (the required LoD tolerance).
        assert abs(grid.donor_fraction - exact.donor_fraction) < 1e-4

    def test_grid_with_bias(self) -> None:
        markers = _make_markers_overdispersed(0.02, n_markers=120, dp=1000, seed=3)
        rng = random.Random(11)
        biases = generate_marker_biases_realistic(len(markers), rng)
        cal = PanelCalibration(
            biases={(m.chrom, m.pos, m.ref, m.alt): b for m, b in zip(markers, biases)}
        )
        exact = estimate_single_donor_bb(markers, error_rate=0.01, grid_steps=201, calibration=cal)
        grid = estimate_single_donor_bb_grid(markers, error_rate=0.01, calibration=cal)
        assert abs(grid.donor_fraction - exact.donor_fraction) < 1e-4

    def test_grid_with_asymmetric_error_table(self) -> None:
        # Per-direction error table (PanelCalibration.errors): the grid must
        # match the exact estimator on the asymmetric REF/ALT-only likelihood
        # path, not just the symmetric error_rate fallback.
        markers = _make_markers_for_fraction(0.02, n_markers=120, dp=1000, seed=17)
        rng = random.Random(23)
        errors = {
            (m.chrom, m.pos, m.ref, m.alt): MarkerErrorRates(
                e_refalt=0.0005 + rng.random() * 0.002,
                e_altref=0.0005 + rng.random() * 0.002,
            )
            for m in markers
        }
        cal = PanelCalibration(errors=errors)
        exact = estimate_single_donor_bb(markers, error_rate=0.01, grid_steps=201, calibration=cal)
        grid = estimate_single_donor_bb_grid(markers, error_rate=0.01, calibration=cal)
        assert abs(grid.donor_fraction - exact.donor_fraction) < 1e-4

    def test_grid_empty_markers(self) -> None:
        res = estimate_single_donor_bb_grid([], error_rate=0.01)
        assert res.donor_fraction == 0.0
        assert res.n_informative == 0

    def test_grid_p_alt_matches_scalar(self) -> None:
        # _p_alt_grid row i must equal _p_alt_for_f(f_grid[i]) exactly.
        markers = _make_markers_overdispersed(0.03, n_markers=40, dp=500, seed=5)
        arr = _precompute_marker_arrays(markers, PanelCalibration())
        f_grid = np.linspace(0.0, 1.0, 11)
        grid_pa = _p_alt_grid(arr, f_grid, 0.01)
        for i, f in enumerate(f_grid):
            scalar = _p_alt_for_f(arr, float(f), 0.01)
            np.testing.assert_allclose(grid_pa[i], scalar, rtol=0, atol=1e-12)


class TestGridTwoRho:
    """The fast grid two-rho path must match the exact two-rho MLE (issue #33)."""

    @pytest.mark.parametrize("f_host", [0.0, 0.001, 0.005, 0.05])
    def test_grid_two_rho_matches_exact_two_rho(self, f_host: float) -> None:
        # Panel-sized mixed-class set with injected het overdispersion.
        markers = _make_mixed_class_markers(
            f_host=f_host, n_hom=150, n_het=100, dp=2000, het_overdisp_rho=71, seed=9
        )
        exact = estimate_single_donor_bb(
            markers, error_rate=0.01, grid_steps=201, marker_type_overdispersion=True
        )
        grid = estimate_single_donor_bb_grid(
            markers, error_rate=0.01, marker_type_overdispersion=True
        )
        # The two-rho path engaged (both classes ample) and matches to < 1e-3.
        assert exact.rho_hom is not None and grid.rho_hom is not None
        assert abs(grid.donor_fraction - exact.donor_fraction) < 1e-3

    def test_grid_shared_rho_path_byte_identical(self) -> None:
        # The opt-out single-rho grid is identical to the pre-#33 path. On this
        # hom-only fixture the default (two-rho on) also falls back, so it matches.
        markers = _make_markers_overdispersed(0.02, n_markers=120, dp=1000, seed=3)
        shared = estimate_single_donor_bb_grid(
            markers, error_rate=0.01, marker_type_overdispersion=False
        )
        default = estimate_single_donor_bb_grid(markers, error_rate=0.01)  # default on
        assert default.donor_fraction == shared.donor_fraction
        assert default.donor_fraction_ci == shared.donor_fraction_ci
        assert default.log_likelihood == shared.log_likelihood
        assert default.rho == shared.rho
        assert shared.rho_hom is None and shared.rho_het is None

    def test_grid_two_rho_sparse_falls_back(self) -> None:
        # Too few het markers: the grid falls back to the single-rho path and
        # records the reason, matching the single-rho grid fraction.
        markers = _make_mixed_class_markers(
            f_host=0.0, n_hom=120, n_het=10, dp=2000, het_overdisp_rho=71, seed=13
        )
        single = estimate_single_donor_bb_grid(
            markers, error_rate=0.01, marker_type_overdispersion=False
        )
        default = estimate_single_donor_bb_grid(markers, error_rate=0.01)  # default on
        assert default.donor_fraction == single.donor_fraction
        assert default.rho_hom is None and default.rho_het is None
        assert default.marker_type_overdispersion_fallback is not None
