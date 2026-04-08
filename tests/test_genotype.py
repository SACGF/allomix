"""Tests for allomix.genotype — VCF parsing and marker classification."""

from pathlib import Path

import pytest

from allomix.genotype import (
    MarkerData,
    _alt_dose,
    classify_markers,
    marker_type,
    parse_vcf,
)

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
EXAMPLE_VCF = DATA_DIR / "idt_rhampseq_sid_example.vcf"
JOINT_VCF = DATA_DIR / "joint_called_example.vcf"


# ---------------------------------------------------------------------------
# marker_type classification
# ---------------------------------------------------------------------------


class TestMarkerType:
    """Test Vynck marker type classification."""

    def test_hom_ref_vs_hom_alt(self):
        assert marker_type((0, 0), (1, 1)) == 0

    def test_hom_alt_vs_hom_ref(self):
        assert marker_type((1, 1), (0, 0)) == 1

    def test_het_vs_hom_ref(self):
        assert marker_type((0, 1), (0, 0)) == 10

    def test_het_vs_hom_alt(self):
        assert marker_type((0, 1), (1, 1)) == 11

    def test_hom_ref_vs_het(self):
        assert marker_type((0, 0), (0, 1)) == 20

    def test_hom_alt_vs_het(self):
        assert marker_type((1, 1), (0, 1)) == 21

    def test_same_hom_ref(self):
        assert marker_type((0, 0), (0, 0)) is None

    def test_same_het(self):
        assert marker_type((0, 1), (0, 1)) is None

    def test_same_hom_alt(self):
        assert marker_type((1, 1), (1, 1)) is None


# ---------------------------------------------------------------------------
# _alt_dose
# ---------------------------------------------------------------------------


class TestAltDose:
    def test_hom_ref(self):
        assert _alt_dose((0, 0)) == 0

    def test_het(self):
        assert _alt_dose((0, 1)) == 1

    def test_hom_alt(self):
        assert _alt_dose((1, 1)) == 2


# ---------------------------------------------------------------------------
# parse_vcf — single-sample example
# ---------------------------------------------------------------------------


class TestParseVcfSingleSample:
    """Test parsing the de-identified single-sample example VCF."""

    @pytest.fixture
    def markers(self):
        return parse_vcf(EXAMPLE_VCF)

    def test_count(self, markers):
        # 15 data lines, but some are hom-ref with ALT="." — those should
        # still parse (AD has ref count only, ad_alt=0)
        assert len(markers) > 0

    def test_has_gt(self, markers):
        for m in markers:
            assert m.gt[0] >= 0 and m.gt[1] >= 0

    def test_has_ad(self, markers):
        for m in markers:
            assert m.ad_ref >= 0

    def test_het_has_alt_counts(self, markers):
        hets = [m for m in markers if m.gt == (0, 1)]
        assert len(hets) > 0
        for m in hets:
            assert m.ad_alt > 0

    def test_hom_alt_has_alt_counts(self, markers):
        hom_alts = [m for m in markers if m.gt == (1, 1)]
        assert len(hom_alts) > 0
        for m in hom_alts:
            assert m.ad_alt > 0

    def test_depths_are_high(self, markers):
        for m in markers:
            assert m.dp > 400

    def test_all_pass(self, markers):
        for m in markers:
            assert m.filter == "PASS"

    def test_gq_present(self, markers):
        for m in markers:
            assert m.gq is not None

    def test_min_dp_filter(self):
        all_markers = parse_vcf(EXAMPLE_VCF, min_dp=0)
        filtered = parse_vcf(EXAMPLE_VCF, min_dp=3000)
        assert len(filtered) < len(all_markers)


# ---------------------------------------------------------------------------
# classify_markers
# ---------------------------------------------------------------------------


class TestClassifyMarkers:
    """Test marker classification with synthetic MarkerData."""

    def _make_marker(self, chrom="chr1", pos=100, gt=(0, 0), ad_ref=1000, ad_alt=0, dp=1000):
        return MarkerData(
            chrom=chrom,
            pos=pos,
            ref="A",
            alt="G",
            gt=gt,
            ad_ref=ad_ref,
            ad_alt=ad_alt,
            dp=dp,
            gq=99,
        )

    def test_fully_informative_type_0(self):
        host = [self._make_marker(gt=(0, 0))]
        donor = [self._make_marker(gt=(1, 1))]
        admix = [self._make_marker(ad_ref=900, ad_alt=100, dp=1000)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert len(result.informative) == 1
        assert result.informative[0].marker_type == 0

    def test_fully_informative_type_1(self):
        host = [self._make_marker(gt=(1, 1))]
        donor = [self._make_marker(gt=(0, 0))]
        admix = [self._make_marker(ad_ref=100, ad_alt=900, dp=1000)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert len(result.informative) == 1
        assert result.informative[0].marker_type == 1

    def test_non_informative_same_gt(self):
        host = [self._make_marker(gt=(0, 1))]
        donor = [self._make_marker(gt=(0, 1))]
        admix = [self._make_marker(ad_ref=500, ad_alt=500, dp=1000)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert len(result.informative) == 0
        assert len(result.non_informative) == 1

    def test_depth_filter(self):
        host = [self._make_marker(gt=(0, 0))]
        donor = [self._make_marker(gt=(1, 1))]
        admix = [self._make_marker(dp=50)]
        result = classify_markers(host, [donor], admix, min_dp=100, min_gq=0)
        assert len(result.informative) == 0
        assert result.n_filtered == 1

    def test_gq_filter(self):
        host = [self._make_marker(gt=(0, 0))]
        host[0].gq = 10
        donor = [self._make_marker(gt=(1, 1))]
        admix = [self._make_marker(dp=1000)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=20)
        assert len(result.informative) == 0
        assert result.n_filtered == 1

    def test_shared_only(self):
        """Only markers present in all VCFs are considered."""
        host = [
            self._make_marker(pos=100, gt=(0, 0)),
            self._make_marker(pos=200, gt=(0, 0)),
        ]
        donor = [self._make_marker(pos=100, gt=(1, 1))]
        admix = [self._make_marker(pos=100, dp=1000)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert result.n_shared == 1
        assert len(result.informative) == 1

    def test_multiple_markers(self):
        """Mix of informative and non-informative markers."""
        host = [
            self._make_marker(pos=100, gt=(0, 0)),
            self._make_marker(pos=200, gt=(0, 1)),
            self._make_marker(pos=300, gt=(1, 1)),
        ]
        donor = [
            self._make_marker(pos=100, gt=(1, 1)),
            self._make_marker(pos=200, gt=(0, 1)),
            self._make_marker(pos=300, gt=(0, 0)),
        ]
        admix = [
            self._make_marker(pos=100, ad_ref=900, ad_alt=100, dp=1000),
            self._make_marker(pos=200, ad_ref=500, ad_alt=500, dp=1000),
            self._make_marker(pos=300, ad_ref=100, ad_alt=900, dp=1000),
        ]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert len(result.informative) == 2  # pos 100 and 300
        assert len(result.non_informative) == 1  # pos 200

    def test_two_donors(self):
        """Multi-donor: marker informative for first donor."""
        host = [self._make_marker(gt=(0, 0))]
        donor1 = [self._make_marker(gt=(1, 1))]
        donor2 = [self._make_marker(gt=(0, 0))]
        admix = [self._make_marker(dp=1000)]
        result = classify_markers(host, [donor1, donor2], admix, min_dp=0, min_gq=0)
        assert len(result.informative) == 1
        assert result.informative[0].donor_gts == [(1, 1), (0, 0)]

    def test_sample_name_not_chromosome(self):
        """classify_markers should not set sample_name to a chromosome."""
        host = [MarkerData("chr7", 100, "A", "T", (0, 0), 100, 0, 100, 99)]
        donor = [MarkerData("chr7", 100, "A", "T", (1, 1), 0, 100, 100, 99)]
        admix = [MarkerData("chr7", 100, "A", "T", (0, 1), 50, 50, 100, 99)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert result.sample_name != "chr7"
