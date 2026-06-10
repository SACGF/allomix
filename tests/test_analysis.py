"""Tests for allomix.analysis wiring, focused on the contamination LoD floor.

The contamination-floor-into-LoD logic (Obs 2 in
``claude/further_improvements.md``) is validated here against a synthetic
contamination scalar, which isolates the flooring rule from the in-data
contamination estimator (tested separately in ``test_contamination.py``).
"""

import math

from allomix.analysis import _floor_detection_limits
from allomix.chimerism import ChimerismResult, MultiDonorResult


def _single_result(lob: float, lod: float) -> ChimerismResult:
    """A minimal single-donor result carrying only the LoB/LoD under test."""
    return ChimerismResult(
        donor_fraction=0.99,
        donor_fraction_ci=(0.98, 1.0),
        host_fraction=0.01,
        log_likelihood=0.0,
        n_informative=50,
        n_markers_used=50,
        per_marker=[],
        error_rate=1e-3,
        lob_fraction=lob,
        lod_fraction=lod,
    )


class TestFloorDetectionLimits:
    def test_floor_raises_both_when_above(self):
        """A floor above the analytical LoD replaces both limits."""
        res = _single_result(lob=0.0008, lod=0.0015)
        _floor_detection_limits(res, 0.002)
        assert res.lob_fraction == 0.002
        assert res.lod_fraction == 0.002

    def test_floor_below_lod_leaves_lod(self):
        """A floor between LoB and LoD raises only LoB; LoD stays and remains the
        larger of the two.
        """
        res = _single_result(lob=0.0008, lod=0.005)
        _floor_detection_limits(res, 0.002)
        assert res.lob_fraction == 0.002
        assert res.lod_fraction == 0.005
        assert res.lod_fraction >= res.lob_fraction

    def test_floor_below_both_is_noop(self):
        res = _single_result(lob=0.003, lod=0.006)
        _floor_detection_limits(res, 0.001)
        assert res.lob_fraction == 0.003
        assert res.lod_fraction == 0.006

    def test_zero_floor_is_noop(self):
        res = _single_result(lob=0.003, lod=0.006)
        _floor_detection_limits(res, 0.0)
        assert res.lob_fraction == 0.003
        assert res.lod_fraction == 0.006

    def test_infinite_lod_stays_infinite_above_floor(self):
        """An uninformative sample (inf LoD) stays inf; a finite floor cannot
        lower it.
        """
        res = _single_result(lob=float("inf"), lod=float("inf"))
        _floor_detection_limits(res, 0.002)
        assert math.isinf(res.lob_fraction)
        assert math.isinf(res.lod_fraction)

    def test_multidonor_result_unchanged(self):
        """Multi-donor results carry no LoB/LoD fields, so flooring is a no-op
        and must not raise.
        """
        res = MultiDonorResult(
            donor_fractions=[0.5, 0.49],
            donor_fraction_cis=[(0.45, 0.55), (0.44, 0.54)],
            host_fraction=0.01,
            log_likelihood=0.0,
            n_informative=50,
            n_markers_used=50,
            per_marker=[],
            error_rate=1e-3,
        )
        _floor_detection_limits(res, 0.002)  # must not raise
        assert not hasattr(res, "lob_fraction")
