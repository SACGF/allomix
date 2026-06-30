"""Tests for multi-donor chimerism estimation.

Covers:
- expected_weight_multi and total_log_likelihood_multi functions
- estimate_multi_donor end-to-end (grid search + Nelder-Mead + profile CIs)
- Sibling trio genotype generation
- Multi-donor VCF blending
- Multi-donor QC, report, and CLI integration
"""

import json
import math
import random
import tempfile
from io import StringIO
from pathlib import Path

import pytest

from allomix.cli import main
from allomix.estimate.chimerism import estimate_multi_donor
from allomix.estimate.likelihood import (
    expected_weight,
    expected_weight_multi,
    total_log_likelihood_multi_bb,
)
from allomix.genotype import InformativeMarker, MarkerType, classify_markers, parse_vcf
from allomix.qc.qc import assess_quality
from allomix.report.report import timeline_json, to_json, to_tsv
from allomix.results import MultiDonorResult
from allomix.simulate import (
    _mendelian_child,
    alt_dose,
    blend_from_genotype_dicts,
    build_joint_vcf,
    expected_vaf_multi,
    generate_sibling_trio_genotypes,
    write_genotype_vcf,
    write_joint_vcf,
    write_vcf,
)

# Path to pre-generated multi-donor test data
MULTIDONOR_DIR = Path(__file__).resolve().parent.parent / "tests" / "test_data" / "multidonor"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_marker(
    host_gt: tuple[int, int],
    donor1_gt: tuple[int, int],
    donor2_gt: tuple[int, int],
    ad_ref: int,
    ad_alt: int,
    pos: int = 100,
) -> InformativeMarker:
    """Create an InformativeMarker for testing."""
    mt1 = MarkerType.classify(host_gt, donor1_gt)
    mt2 = MarkerType.classify(host_gt, donor2_gt)
    mtype = mt1 if mt1 is not None else mt2
    return InformativeMarker(
        chrom="chr1",
        pos=pos,
        ref="A",
        alt="G",
        host_gt=host_gt,
        donor_gts=[donor1_gt, donor2_gt],
        marker_type=mtype if mtype is not None else 0,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
        marker_types=[mt1, mt2],
        informative_for=[mt1 is not None, mt2 is not None],
    )


def _make_sibling_markers_and_blend(
    f1: float,
    f2: float,
    n_markers: int = 100,
    depth: int = 2000,
    seed: int = 42,
) -> list[InformativeMarker]:
    """Generate sibling trio genotypes, blend, and return InformativeMarkers."""
    rng = random.Random(seed)
    markers = generate_sibling_trio_genotypes(n_markers, rng)

    blend = blend_from_genotype_dicts(
        markers,
        [f1, f2],
        target_depth=depth,
        seed=seed + 1,
        error_rate=0.01,
    )

    # Write to tmp files and parse with cyvcf2 via genotype module
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        write_genotype_vcf(markers, tmpdir / "host.vcf", "HOST", key="host_gt")
        write_genotype_vcf(markers, tmpdir / "d1.vcf", "D1", key="donor1_gt")
        write_genotype_vcf(markers, tmpdir / "d2.vcf", "D2", key="donor2_gt")
        write_vcf(blend, tmpdir / "admix.vcf")

        host = parse_vcf(tmpdir / "host.vcf", min_dp=0, min_gq=0)
        d1 = parse_vcf(tmpdir / "d1.vcf", min_dp=0, min_gq=0)
        d2 = parse_vcf(tmpdir / "d2.vcf", min_dp=0, min_gq=0)
        admix = parse_vcf(tmpdir / "admix.vcf", min_dp=0, min_gq=0)

        genotypes = classify_markers(host, [d1, d2], admix, min_dp=0, min_gq=0, pass_only=False)
        return genotypes.informative


# ---------------------------------------------------------------------------
# Simulate module tests
# ---------------------------------------------------------------------------


class TestMendelianChild:
    """Test Mendelian segregation helper."""

    def test_hom_parents_produce_hom_child(self):
        rng = random.Random(1)
        for _ in range(50):
            child = _mendelian_child((0, 0), (0, 0), rng)
            assert child == (0, 0)

    def test_het_het_produces_all_genotypes(self):
        rng = random.Random(42)
        gts = set()
        for _ in range(200):
            gts.add(_mendelian_child((0, 1), (0, 1), rng))
        assert (0, 0) in gts
        assert (0, 1) in gts
        assert (1, 1) in gts

    def test_sorted_output(self):
        rng = random.Random(99)
        for _ in range(100):
            gt = _mendelian_child((1, 0), (1, 0), rng)
            assert gt[0] <= gt[1]


class TestSiblingTrioGenotypes:
    """Test sibling trio genotype generation."""

    def test_correct_length(self):
        markers = generate_sibling_trio_genotypes(50, random.Random(1))
        assert len(markers) == 50

    def test_required_keys(self):
        markers = generate_sibling_trio_genotypes(10, random.Random(1))
        required = {
            "chrom",
            "pos",
            "ref",
            "alt",
            "host_gt",
            "donor1_gt",
            "donor2_gt",
            "p_alt",
            "informative_d1",
            "informative_d2",
            "informative_any",
            "donors_distinguishable",
        }
        assert required.issubset(markers[0].keys())

    def test_informativity_flags_correct(self):
        markers = generate_sibling_trio_genotypes(100, random.Random(42))
        for m in markers:
            assert m["informative_d1"] == (alt_dose(m["host_gt"]) != alt_dose(m["donor1_gt"]))
            assert m["informative_d2"] == (alt_dose(m["host_gt"]) != alt_dose(m["donor2_gt"]))
            assert m["informative_any"] == (m["informative_d1"] or m["informative_d2"])

    def test_siblings_have_some_shared_genotypes(self):
        """Siblings should share genotype at some loci (IBD=2 ~25%)."""
        markers = generate_sibling_trio_genotypes(200, random.Random(42))
        n_all_same = sum(1 for m in markers if m["host_gt"] == m["donor1_gt"] == m["donor2_gt"])
        # With 200 markers, expect ~50 all-same (25% IBD=2 for each pair,
        # plus population homozygosity). Should be at least 10.
        assert n_all_same > 10


class TestExpectedVafMulti:
    """Test multi-donor expected VAF function."""

    def test_single_donor_equivalent(self):
        """With one donor fraction=0, should match single-donor."""
        host = (0, 0)
        d1 = (1, 1)
        d2 = (0, 1)
        for f in [0.0, 0.1, 0.5, 1.0]:
            single = ((1.0 - f) * alt_dose(host) + f * alt_dose(d1)) / 2.0
            multi = expected_vaf_multi(host, [d1, d2], [f, 0.0])
            assert multi == pytest.approx(single)

    def test_pure_host(self):
        assert expected_vaf_multi((0, 0), [(1, 1), (0, 1)], [0.0, 0.0]) == 0.0

    def test_pure_donor1(self):
        assert expected_vaf_multi((0, 0), [(1, 1), (0, 1)], [1.0, 0.0]) == 1.0

    def test_pure_donor2(self):
        assert expected_vaf_multi((0, 0), [(1, 1), (0, 1)], [0.0, 1.0]) == 0.5

    def test_three_way_mix(self):
        # host 0/0 (dose=0), d1 1/1 (dose=2), d2 0/1 (dose=1)
        # f1=0.25, f2=0.25 -> host=0.5
        # vaf = (0.5*0 + 0.25*2 + 0.25*1) / 2 = 0.75/2 = 0.375
        assert expected_vaf_multi((0, 0), [(1, 1), (0, 1)], [0.25, 0.25]) == pytest.approx(0.375)


class TestBlendFromGenotypeDicts:
    """Test multi-donor VCF blending from genotype dicts."""

    def test_basic_blend(self):
        markers = generate_sibling_trio_genotypes(20, random.Random(42))
        result = blend_from_genotype_dicts(markers, [0.3, 0.2], target_depth=500, seed=1)
        assert result.num_markers == 20
        assert len(result.records) == 20

    def test_fractions_sum_exceeds_one_raises(self):
        markers = generate_sibling_trio_genotypes(5, random.Random(1))
        with pytest.raises(ValueError, match="must be <= 1.0"):
            blend_from_genotype_dicts(markers, [0.6, 0.5])


# ---------------------------------------------------------------------------
# Core likelihood tests
# ---------------------------------------------------------------------------


class TestExpectedWeightMulti:
    """Test multi-donor expected weight function."""

    def test_pure_host(self):
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.0, 0.0])
        assert w == pytest.approx(1.0)

    def test_pure_donor1(self):
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [1.0, 0.0])
        assert w == pytest.approx(0.0)

    def test_pure_donor2(self):
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.0, 1.0])
        assert w == pytest.approx(0.5)

    def test_equal_mix(self):
        # host 0/0 (ref=2), d1 1/1 (ref=0), d2 0/1 (ref=1)
        # f1=0.25, f2=0.25, fh=0.5
        # w = 0.5*2/2 + 0.25*0/2 + 0.25*1/2 = 0.5 + 0 + 0.125 = 0.625
        w = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.25, 0.25])
        assert w == pytest.approx(0.625)

    def test_reduces_to_single_donor(self):
        """With f2=0, should match single-donor expected_weight."""
        for f in [0.0, 0.1, 0.5, 1.0]:
            w_single = expected_weight((0, 0), (1, 1), f)
            w_multi = expected_weight_multi((0, 0), [(1, 1), (0, 0)], [f, 0.0])
            assert w_multi == pytest.approx(w_single)

    def test_bias_applied(self):
        # Bias is applied multiplicatively in logit space (issue #20), not as a
        # flat additive shift. In odds space the correction is exact:
        #   w' = odds(w) / (odds(w) + odds(0.5 + bias))
        # With w = 0.625 and bias = 0.05: odds(0.625)=5/3, odds(0.55)=11/9, so
        #   w' = (5/3) / (5/3 + 11/9) = 15/26.
        w_no_bias = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.25, 0.25])
        assert w_no_bias == pytest.approx(0.625)
        w_bias = expected_weight_multi((0, 0), [(1, 1), (0, 1)], [0.25, 0.25], bias=0.05)
        assert w_bias < w_no_bias
        assert w_bias == pytest.approx(15 / 26)


class TestMultiDonorLikelihood:
    """Test 2D log-likelihood computation."""

    def test_ll_is_finite(self):
        m = _make_marker((0, 0), (1, 1), (0, 1), ad_ref=800, ad_alt=200)
        ll = total_log_likelihood_multi_bb([m], [0.2, 0.1], error_rate=0.01)
        assert math.isfinite(ll)

    def test_ll_peaks_near_truth(self):
        """LL should be higher at the true fractions than far away."""
        markers = _make_sibling_markers_and_blend(0.30, 0.10, seed=42)
        ll_truth = total_log_likelihood_multi_bb(markers, [0.30, 0.10])
        ll_wrong = total_log_likelihood_multi_bb(markers, [0.10, 0.30])
        ll_zero = total_log_likelihood_multi_bb(markers, [0.0, 0.0])
        assert ll_truth > ll_wrong
        assert ll_truth > ll_zero


# ---------------------------------------------------------------------------
# Multi-donor estimation tests
# ---------------------------------------------------------------------------


class TestEstimateMultiDonor:
    """End-to-end multi-donor estimation tests."""

    def test_empty_markers(self):
        result = estimate_multi_donor([], n_donors=2)
        assert result.donor_fractions == [0.0, 0.0]
        assert result.host_fraction == 1.0
        assert result.n_informative == 0

    def test_pure_host(self):
        """f1=f2=0: both donors should estimate near 0%."""
        markers = _make_sibling_markers_and_blend(0.0, 0.0)
        result = estimate_multi_donor(markers)
        assert result.donor_fractions[0] < 0.03
        assert result.donor_fractions[1] < 0.03
        assert result.host_fraction > 0.94

    def test_single_donor_only(self):
        """f1=0.20, f2=0: should recover donor1~20%, donor2~0%."""
        markers = _make_sibling_markers_and_blend(0.20, 0.0, seed=100)
        result = estimate_multi_donor(markers)
        assert 0.10 < result.donor_fractions[0] < 0.35
        assert result.donor_fractions[1] < 0.08

    def test_balanced_mix(self):
        """f1=f2=0.25: should recover both donors near 25%."""
        markers = _make_sibling_markers_and_blend(0.25, 0.25, seed=200)
        result = estimate_multi_donor(markers)
        assert 0.15 < result.donor_fractions[0] < 0.35
        assert 0.15 < result.donor_fractions[1] < 0.35

    def test_asymmetric_mix(self):
        """f1=0.30, f2=0.10: should distinguish major vs minor donor."""
        markers = _make_sibling_markers_and_blend(0.30, 0.10, seed=300)
        result = estimate_multi_donor(markers)
        assert result.donor_fractions[0] > result.donor_fractions[1]

    def test_fractions_sum_le_one(self):
        """Estimated fractions must satisfy f1 + f2 <= 1."""
        markers = _make_sibling_markers_and_blend(0.40, 0.40, seed=400)
        result = estimate_multi_donor(markers)
        assert sum(result.donor_fractions) <= 1.0 + 1e-9

    def test_ci_structure(self):
        """CIs should be valid (lo <= mle <= hi for each donor)."""
        markers = _make_sibling_markers_and_blend(0.25, 0.25, seed=500)
        result = estimate_multi_donor(markers)
        for i in range(2):
            lo, hi = result.donor_fraction_cis[i]
            assert lo <= result.donor_fractions[i] + 1e-6
            assert hi >= result.donor_fractions[i] - 1e-6
            assert lo >= 0.0
            assert hi <= 1.0

    def test_ci_contains_truth(self):
        """Profile likelihood CIs should contain the true fractions."""
        f1_true, f2_true = 0.25, 0.15
        markers = _make_sibling_markers_and_blend(f1_true, f2_true, depth=3000, seed=600)
        result = estimate_multi_donor(markers)
        lo1, hi1 = result.donor_fraction_cis[0]
        lo2, hi2 = result.donor_fraction_cis[1]
        assert lo1 <= f1_true <= hi1, f"Donor1 CI [{lo1:.4f}, {hi1:.4f}] misses truth {f1_true}"
        assert lo2 <= f2_true <= hi2, f"Donor2 CI [{lo2:.4f}, {hi2:.4f}] misses truth {f2_true}"

    def test_per_marker_results_populated(self):
        markers = _make_sibling_markers_and_blend(0.20, 0.10, seed=700)
        result = estimate_multi_donor(markers)
        assert len(result.per_marker) == len(markers)
        for mr in result.per_marker:
            assert mr.dp > 0

    def test_per_donor_n_informative(self):
        markers = _make_sibling_markers_and_blend(0.20, 0.10, seed=800)
        result = estimate_multi_donor(markers)
        assert result.per_donor_n_informative is not None
        assert len(result.per_donor_n_informative) == 2
        assert all(n >= 0 for n in result.per_donor_n_informative)

    def test_result_is_multi_donor_type(self):
        markers = _make_sibling_markers_and_blend(0.20, 0.10, seed=900)
        result = estimate_multi_donor(markers)
        assert isinstance(result, MultiDonorResult)


# ---------------------------------------------------------------------------
# Integration tests with pre-generated VCFs
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not MULTIDONOR_DIR.exists(),
    reason="Multi-donor test data not generated",
)
class TestMultiDonorIntegration:
    """Integration tests using pre-generated 3-brothers VCFs."""

    @pytest.fixture(scope="class")
    @staticmethod
    def genotype_data():
        host = parse_vcf(MULTIDONOR_DIR / "host.vcf", min_dp=0, min_gq=0)
        d1 = parse_vcf(MULTIDONOR_DIR / "donor1.vcf", min_dp=0, min_gq=0)
        d2 = parse_vcf(MULTIDONOR_DIR / "donor2.vcf", min_dp=0, min_gq=0)
        return host, d1, d2

    def _run(self, genotype_data, name):
        host, d1, d2 = genotype_data
        admix = parse_vcf(MULTIDONOR_DIR / f"{name}.vcf", min_dp=0, min_gq=0)
        genotypes = classify_markers(
            host,
            [d1, d2],
            admix,
            min_dp=0,
            min_gq=0,
            pass_only=False,
        )
        genotypes.sample_name = name
        result = estimate_multi_donor(genotypes.informative, n_donors=2)
        qc = assess_quality(result, genotypes)
        return result, qc, genotypes

    def test_pure_host(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_100_d1_0_d2_0")
        assert result.donor_fractions[0] < 0.03
        assert result.donor_fractions[1] < 0.03

    def test_balanced_25_25(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_50_d1_25_d2_25")
        assert 0.15 < result.donor_fractions[0] < 0.35
        assert 0.15 < result.donor_fractions[1] < 0.35

    def test_asymmetric_30_10(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_60_d1_30_d2_10")
        assert result.donor_fractions[0] > result.donor_fractions[1]
        assert 0.20 < result.donor_fractions[0] < 0.40

    def test_asymmetric_10_30(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_60_d1_10_d2_30")
        assert result.donor_fractions[1] > result.donor_fractions[0]
        assert 0.20 < result.donor_fractions[1] < 0.40

    def test_pure_donor1(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_0_d1_100_d2_0")
        assert result.donor_fractions[0] > 0.90

    def test_pure_donor2(self, genotype_data):
        result, _, _ = self._run(genotype_data, "host_0_d1_0_d2_100")
        assert result.donor_fractions[1] > 0.90

    def test_qc_passes(self, genotype_data):
        _, qc, _ = self._run(genotype_data, "host_50_d1_25_d2_25")
        assert qc.pass_
        assert qc.per_donor_n_informative is not None

    def test_json_output(self, genotype_data):
        result, qc, genotypes = self._run(genotype_data, "host_50_d1_25_d2_25")
        data = to_json(result, qc, sample_name=genotypes.sample_name)
        assert "donors" in data
        assert len(data["donors"]) == 2
        assert "host_pct" in data
        # Should be serialisable
        json.dumps(data)

    def test_tsv_output(self, genotype_data):
        result, qc, _ = self._run(genotype_data, "host_50_d1_25_d2_25")
        buf = StringIO()
        to_tsv(result, qc, buf)
        content = buf.getvalue()
        assert "donor1_pct" in content
        assert "donor2_pct" in content
        assert "host_pct" in content

    def test_timeline_output(self, genotype_data):
        results = []
        for name in ["host_80_d1_10_d2_10", "host_50_d1_25_d2_25"]:
            result, qc, genotypes = self._run(genotype_data, name)
            results.append((genotypes.sample_name, result, qc))
        data = timeline_json(results)
        assert len(data["timepoints"]) == 2
        assert "donors" in data["timepoints"][0]
        json.dumps(data)


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not MULTIDONOR_DIR.exists(),
    reason="Multi-donor test data not generated",
)
class TestMultiDonorCLI:
    """Test CLI with multi-donor joint VCF inputs."""

    @pytest.fixture
    def joint_vcf(self, tmp_path):
        """Build a multi-donor joint VCF from the existing separate VCFs."""
        result = build_joint_vcf(
            host_path=str(MULTIDONOR_DIR / "host.vcf"),
            donor_paths=[
                str(MULTIDONOR_DIR / "donor1.vcf"),
                str(MULTIDONOR_DIR / "donor2.vcf"),
            ],
            admix_fractions=[0.25, 0.30, 0.10, 0.20],
            admix_sample_names=["D1_25_D2_25", "D1_30_D2_10", "D1_10_D2_10", "D1_20_D2_0"],
            host_sample_name="HOST",
            donor_sample_names=["DONOR1", "DONOR2"],
            target_depth=2000,
            seed=42,
        )
        path = tmp_path / "joint_multi.vcf"
        write_joint_vcf(result, path)
        return path

    def test_monitor_json(self, tmp_path, joint_vcf):
        out = tmp_path / "result.json"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(joint_vcf),
                "--admix-vcf",
                str(joint_vcf),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR1",
                "--donor-sample",
                "DONOR2",
                "--sample",
                "D1_25_D2_25",
                "--json",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())["analysis"]
        assert "donors" in data
        assert len(data["donors"]) == 2

    def test_monitor_tsv(self, tmp_path, joint_vcf):
        out = tmp_path / "result.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(joint_vcf),
                "--admix-vcf",
                str(joint_vcf),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR1",
                "--donor-sample",
                "DONOR2",
                "--sample",
                "D1_30_D2_10",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        content = out.read_text()
        assert "donor1_pct" in content

    def test_timeline(self, tmp_path, joint_vcf):
        out = tmp_path / "timeline.json"
        rc = main(
            [
                "timeline",
                "--genotype-vcf",
                str(joint_vcf),
                "--admix-vcf",
                str(joint_vcf),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR1",
                "--donor-sample",
                "DONOR2",
                "--sample",
                "D1_10_D2_10",
                "--sample",
                "D1_25_D2_25",
                "--json",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert len(data["timepoints"]) == 2

    def test_single_donor_still_works(self, tmp_path, joint_vcf):
        """Single --donor-sample should still use the single-donor estimator."""
        out = tmp_path / "single.json"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(joint_vcf),
                "--admix-vcf",
                str(joint_vcf),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR1",
                "--sample",
                "D1_20_D2_0",
                "--json",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())["analysis"]
        # Single-donor result should have donor_pct, not donors[]
        assert "donor_pct" in data
        assert "donors" not in data


# ---------------------------------------------------------------------------
# Genotype module multi-donor tests
# ---------------------------------------------------------------------------


class TestClassifyMarkersMultiDonor:
    """Test classify_markers with 2 donors."""

    def test_marker_informative_for_second_donor_only(self, tmp_path):
        """A marker where host==donor1 but host!=donor2 should be included."""
        markers = [
            {
                "chrom": "chr1",
                "pos": 1000,
                "ref": "A",
                "alt": "G",
                "host_gt": (0, 0),
                "donor1_gt": (0, 0),
                "donor2_gt": (1, 1),
                "p_alt": 0.3,
                "informative_d1": False,
                "informative_d2": True,
                "informative_any": True,
                "donors_distinguishable": True,
            },
        ]
        write_genotype_vcf(markers, tmp_path / "h.vcf", "H", key="host_gt")
        write_genotype_vcf(markers, tmp_path / "d1.vcf", "D1", key="donor1_gt")
        write_genotype_vcf(markers, tmp_path / "d2.vcf", "D2", key="donor2_gt")
        # Create admix VCF with some alt reads
        blend = blend_from_genotype_dicts(
            markers,
            [0.0, 0.3],
            target_depth=1000,
            seed=1,
        )
        write_vcf(blend, tmp_path / "admix.vcf")

        host = parse_vcf(tmp_path / "h.vcf", min_dp=0, min_gq=0)
        d1 = parse_vcf(tmp_path / "d1.vcf", min_dp=0, min_gq=0)
        d2 = parse_vcf(tmp_path / "d2.vcf", min_dp=0, min_gq=0)
        admix = parse_vcf(tmp_path / "admix.vcf", min_dp=0, min_gq=0)

        genotypes = classify_markers(
            host,
            [d1, d2],
            admix,
            min_dp=0,
            min_gq=0,
            pass_only=False,
        )
        assert len(genotypes.informative) == 1
        m = genotypes.informative[0]
        assert m.informative_for == [False, True]
        assert m.marker_types[0] is None  # not informative for donor1
        assert m.marker_types[1] is not None  # informative for donor2


class TestThreeDonorValidation:
    """estimate_multi_donor should reject n_donors > 2 with a clear error."""

    def test_three_donors_raises_value_error(self):
        rng = random.Random(42)
        markers = []
        for i in range(30):
            ad_alt = sum(1 for _ in range(2000) if rng.random() < 0.30)
            markers.append(
                InformativeMarker(
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                    ref="A",
                    alt="T",
                    host_gt=(0, 0),
                    donor_gts=[(1, 1), (0, 1), (0, 1)],
                    marker_type=0,
                    admix_ad_ref=2000 - ad_alt,
                    admix_ad_alt=ad_alt,
                    admix_dp=2000,
                    informative_for=[True, True, True],
                )
            )
        with pytest.raises(ValueError, match="not supported"):
            estimate_multi_donor(markers, n_donors=3)
