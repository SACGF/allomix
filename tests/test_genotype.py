"""Tests for allomix.genotype — VCF parsing and marker classification."""

from pathlib import Path

import pytest

from allomix.genotype import (
    MarkerData,
    MarkerType,
    classify_markers,
    parse_vcf,
)

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"
# Synthetic single-sample VCF (made-up coordinates, fully synthetic counts).
# Deliberately not a real panel VCF: marker positions from a vendor panel are
# proprietary, so the suite must not depend on one.
EXAMPLE_VCF = TEST_DATA_DIR / "single_sample_example.vcf"
JOINT_VCF = TEST_DATA_DIR / "joint_single_donor.vcf"


# ---------------------------------------------------------------------------
# marker_type classification
# ---------------------------------------------------------------------------


class TestMarkerType:
    """Test Vynck marker type classification."""

    def test_hom_ref_vs_hom_alt(self):
        assert MarkerType.classify((0, 0), (1, 1)) == 0

    def test_hom_alt_vs_hom_ref(self):
        assert MarkerType.classify((1, 1), (0, 0)) == 1

    def test_het_vs_hom_ref(self):
        assert MarkerType.classify((0, 1), (0, 0)) == 10

    def test_het_vs_hom_alt(self):
        assert MarkerType.classify((0, 1), (1, 1)) == 11

    def test_hom_ref_vs_het(self):
        assert MarkerType.classify((0, 0), (0, 1)) == 20

    def test_hom_alt_vs_het(self):
        assert MarkerType.classify((1, 1), (0, 1)) == 21

    def test_same_hom_ref(self):
        assert MarkerType.classify((0, 0), (0, 0)) is None

    def test_same_het(self):
        assert MarkerType.classify((0, 1), (0, 1)) is None

    def test_same_hom_alt(self):
        assert MarkerType.classify((1, 1), (1, 1)) is None


# ---------------------------------------------------------------------------
# MarkerType metadata (dose + label carried on the enum members)
# ---------------------------------------------------------------------------


class TestMarkerTypeMetadata:
    def test_dose_pairs_round_trip_through_classify(self):
        # Every member's (host_dose, donor_dose) must classify back to itself.
        for mt in MarkerType:
            host_gt = (0, 0) if mt.host_dose == 0 else (0, 1) if mt.host_dose == 1 else (1, 1)
            donor_gt = (0, 0) if mt.donor_dose == 0 else (0, 1) if mt.donor_dose == 1 else (1, 1)
            assert MarkerType.classify(host_gt, donor_gt) is mt

    def test_label(self):
        assert MarkerType.HOST_HOMREF_DONOR_HOMALT.label == "host 0/0, donor 1/1"
        assert MarkerType.HOST_HOMALT_DONOR_HET.label == "host 1/1, donor 0/1"

    def test_label_for_bare_code(self):
        assert MarkerType.label_for(0) == "host 0/0, donor 1/1"
        assert MarkerType.label_for(999) == "999"  # unrecognised code passes through


# ---------------------------------------------------------------------------
# parse_vcf — single-sample example
# ---------------------------------------------------------------------------


class TestParseVcfSingleSample:
    """Test parsing the synthetic single-sample example VCF."""

    @pytest.fixture
    def markers(self):
        return parse_vcf(EXAMPLE_VCF)

    def test_count(self, markers):
        # 15 data lines, including hom-ref with ALT="." (AD has ref count only,
        # ad_alt=0); all parse.
        assert len(markers) == 15

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

    def test_fixture_integrity(self, markers):
        # Every fixture marker is a PASS call carrying a genotype quality.
        for m in markers:
            assert m.filter == "PASS"
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

    def test_sample_name_defaults_empty(self):
        """classify_markers leaves sample_name empty by default."""
        host = [MarkerData("chr7", 100, "A", "T", (0, 0), 100, 0, 100, 99)]
        donor = [MarkerData("chr7", 100, "A", "T", (1, 1), 0, 100, 100, 99)]
        admix = [MarkerData("chr7", 100, "A", "T", (0, 1), 50, 50, 100, 99)]
        result = classify_markers(host, [donor], admix, min_dp=0, min_gq=0)
        assert result.sample_name == ""


# ---------------------------------------------------------------------------
# parse_vcf — multi-sample joint VCF
# ---------------------------------------------------------------------------


class TestParseVcfMultiSample:
    """Test parsing a multi-sample joint-called VCF with sample selection."""

    @pytest.mark.parametrize("sample", ["HOST", "DONOR", "ADMIX_F0.50"])
    def test_parse_by_name(self, sample):
        markers = parse_vcf(JOINT_VCF, sample=sample)
        assert len(markers) > 0

    def test_different_samples_differ(self):
        host = parse_vcf(JOINT_VCF, sample="HOST")
        donor = parse_vcf(JOINT_VCF, sample="DONOR")
        # Host and donor should have different genotypes at most markers
        n_differ = sum(1 for h, d in zip(host, donor) if h.gt != d.gt)
        assert n_differ > 0

    def test_default_index_zero(self):
        """Default sample=0 should return the same as sample='HOST'."""
        by_idx = parse_vcf(JOINT_VCF, sample=0)
        by_name = parse_vcf(JOINT_VCF, sample="HOST")
        assert len(by_idx) == len(by_name)
        for a, b in zip(by_idx, by_name):
            assert a.gt == b.gt
            assert a.ad_ref == b.ad_ref
            assert a.ad_alt == b.ad_alt

    def test_invalid_sample_raises(self):
        with pytest.raises(ValueError, match="not found"):
            parse_vcf(JOINT_VCF, sample="NONEXISTENT")

    def test_admixture_at_zero_matches_host(self):
        """ADMIX_F0.00 should have genotypes very similar to HOST."""
        host = parse_vcf(JOINT_VCF, sample="HOST")
        admix = parse_vcf(JOINT_VCF, sample="ADMIX_F0.00")
        # At f=0, admixture GTs should mostly match host GTs
        n_match = sum(1 for h, a in zip(host, admix) if h.gt == a.gt)
        assert n_match > len(host) * 0.8
