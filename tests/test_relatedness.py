"""Tests for allomix.relatedness — relatedness estimation and swap detection."""

import random

import pytest

from allomix.analysis import analyse_sample
from allomix.cli import main
from allomix.genotype import MarkerData, parse_vcf
from allomix.relatedness import (
    DEGREE_FIRST,
    DEGREE_IDENTICAL,
    DEGREE_THIRD,
    DEGREE_UNRELATED,
    MIN_CONSENSUS,
    admix_consistency,
    evaluate_expected,
    relatedness_coefficient,
)
from allomix.simulate import (
    blend_vcfs,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

SEEDS = [1, 2, 3, 4, 5]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _md(
    pos: int,
    gt: tuple[int, int],
    ad_ref: int = 0,
    ad_alt: int = 0,
    dp: int = 0,
    chrom: str = "chr1",
) -> MarkerData:
    """Build a MarkerData with the genotype and depths the tests need."""
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        gt=gt,
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=dp,
    )


def _markers_to_md(markers: list[dict], key: str) -> list[MarkerData]:
    """Convert generate_related_genotypes() dicts to a MarkerData list."""
    return [
        MarkerData(
            chrom=m["chrom"],
            pos=m["pos"],
            ref=m["ref"],
            alt=m["alt"],
            gt=tuple(m[key]),
            ad_ref=0,
            ad_alt=0,
            dp=0,
        )
        for m in markers
    ]


def _build_pair(n_shared_het: int, n_ibs0: int) -> tuple[list[MarkerData], list[MarkerData]]:
    """Two MarkerData lists with a known het-share / opposite-hom structure.

    Yields ``het_a == het_b == n_shared_het`` and ``ibs0 == n_ibs0``, so the
    coefficient is ``(n_shared_het - 2 * n_ibs0) / n_shared_het`` exactly and
    ``min(het_a, het_b) == n_shared_het`` drives confidence and the CI width.
    """
    a: list[MarkerData] = []
    b: list[MarkerData] = []
    pos = 0
    for _ in range(n_shared_het):
        pos += 1
        a.append(_md(pos, (0, 1)))
        b.append(_md(pos, (0, 1)))
    for _ in range(n_ibs0):
        pos += 1
        a.append(_md(pos, (0, 0)))
        b.append(_md(pos, (1, 1)))
    return a, b


# ---------------------------------------------------------------------------
# Coefficient over simulated relatedness
# ---------------------------------------------------------------------------


class TestRelatednessCoefficient:
    def test_identical_is_one(self):
        for seed in SEEDS:
            rng = random.Random(seed)
            markers = generate_related_genotypes(400, "unrelated", rng)
            host = _markers_to_md(markers, "host_gt")
            donor = _markers_to_md(markers, "host_gt")  # identical to host
            res = relatedness_coefficient(host, donor, "host", "donor")
            assert res.coefficient == 1.0
            assert res.degree == DEGREE_IDENTICAL
            assert res.relationship.startswith("identical")
            assert res.ibs0 == 0

    def test_unrelated_near_zero(self):
        coefs = []
        for seed in SEEDS:
            rng = random.Random(seed)
            markers = generate_related_genotypes(500, "unrelated", rng)
            host = _markers_to_md(markers, "host_gt")
            donor = _markers_to_md(markers, "donor_gt")
            res = relatedness_coefficient(host, donor, "host", "donor")
            assert res.coefficient is not None
            assert res.coefficient < 0.2  # well clear of the first-degree band
            coefs.append(res.coefficient)
        assert sum(coefs) / len(coefs) < 0.1

    def test_sibling_near_half(self):
        for seed in SEEDS:
            rng = random.Random(seed)
            markers = generate_related_genotypes(600, "sibling", rng)
            host = _markers_to_md(markers, "host_gt")
            donor = _markers_to_md(markers, "donor_gt")
            res = relatedness_coefficient(host, donor, "host", "donor")
            assert res.coefficient is not None
            assert 0.35 < res.coefficient < 0.7
            assert res.degree == DEGREE_FIRST

    def test_too_few_het_sites_gives_no_coefficient(self):
        # Only two shared het sites -> below MIN_HET_SITES.
        host = [_md(1, (0, 1)), _md(2, (0, 1)), _md(3, (0, 0))]
        donor = [_md(1, (0, 1)), _md(2, (0, 1)), _md(3, (0, 0))]
        res = relatedness_coefficient(host, donor, "host", "donor")
        assert res.coefficient is None
        assert res.degree is None
        assert res.relationship == "undetermined"

    def test_sex_chroms_excluded(self):
        # Het sites on chrX must not count toward the autosomal estimate.
        host = [_md(i, (0, 1)) for i in range(10)]
        donor = [_md(i, (0, 1)) for i in range(10)]
        x_host = [_md(i, (0, 1), chrom="chrX") for i in range(10)]
        x_donor = [_md(i, (0, 1), chrom="chrX") for i in range(10)]
        res = relatedness_coefficient(host + x_host, donor + x_donor, "host", "donor")
        assert res.n_sites == 10


class TestConfidenceAndCI:
    def test_confidence_steps_with_site_count(self):
        # min(het) drives the categorical confidence.
        assert relatedness_coefficient(*_build_pair(10, 0), "a", "b").confidence == "low"
        assert relatedness_coefficient(*_build_pair(25, 0), "a", "b").confidence == "med"
        assert relatedness_coefficient(*_build_pair(60, 0), "a", "b").confidence == "high"

    def test_ci_narrows_as_sites_grow(self):
        # Same coefficient (0.5), more sites -> tighter interval.
        small = relatedness_coefficient(*_build_pair(20, 5), "a", "b")
        large = relatedness_coefficient(*_build_pair(80, 20), "a", "b")
        assert abs(small.coefficient - 0.5) < 1e-9
        assert abs(large.coefficient - 0.5) < 1e-9
        small_w = small.ci_high - small.ci_low
        large_w = large.ci_high - large.ci_low
        assert large_w < small_w


# ---------------------------------------------------------------------------
# Expected-vs-detected comparison
# ---------------------------------------------------------------------------


class TestEvaluateExpected:
    def _result_with_degree(self, n_shared_het: int, n_ibs0: int):
        return relatedness_coefficient(*_build_pair(n_shared_het, n_ibs0), "host", "donor")

    def test_no_declaration_returns_none(self):
        res = self._result_with_degree(40, 0)
        assert evaluate_expected(res, None) is None
        assert evaluate_expected(res, "NA") is None
        assert evaluate_expected(res, "na") is None
        assert evaluate_expected(res, "") is None

    def test_unrecognised_declaration_raises(self):
        # A typo must error, not silently skip the check.
        res = self._result_with_degree(100, 25)
        with pytest.raises(ValueError, match="unrecognised expected relatedness"):
            evaluate_expected(res, "frist-degree")
        # "identical" is not an accepted declaration either.
        with pytest.raises(ValueError):
            evaluate_expected(res, "identical")

    def test_declared_related_detected_unrelated_fails(self):
        res = self._result_with_degree(40, 20)  # coef 0.0 -> unrelated
        assert res.degree == DEGREE_UNRELATED
        v = evaluate_expected(res, "first-degree")
        assert v.status == "FAIL"
        assert "first-degree" in v.message and "unrelated" in v.message

    def test_declared_unrelated_detected_identical_fails(self):
        # Identical where two distinct individuals were declared: sample reuse.
        res = self._result_with_degree(40, 0)  # coef 1.0 -> identical
        assert res.degree == DEGREE_IDENTICAL
        v = evaluate_expected(res, "unrelated")
        assert v.status == "FAIL"

    def test_declared_unrelated_detected_first_degree_reviews(self):
        # Gaining a (non-identical) relationship is not a random-swap signature,
        # so it is a REVIEW, not a FAIL.
        res = self._result_with_degree(100, 25)  # coef 0.5 -> first-degree
        assert res.degree == DEGREE_FIRST
        v = evaluate_expected(res, "unrelated")
        assert v.status == "REVIEW"

    def test_declared_first_degree_detected_identical_fails(self):
        # Expected a sibling, got a duplicate: still a reuse FAIL.
        res = self._result_with_degree(40, 0)  # identical
        v = evaluate_expected(res, "first-degree")
        assert v.status == "FAIL"

    def test_within_tolerance_passes(self):
        # coef 0.12 -> third-degree; declared second-degree, distance 1.
        res = self._result_with_degree(100, 44)
        assert res.degree == DEGREE_THIRD
        v = evaluate_expected(res, "second-degree", tolerance=1)
        assert v.status == "PASS"

    def test_two_level_gap_reviews(self):
        # coef 0.12 -> third-degree; declared first-degree, distance 2.
        res = self._result_with_degree(100, 44)
        assert res.degree == DEGREE_THIRD
        v = evaluate_expected(res, "first-degree", tolerance=1)
        assert v.status == "REVIEW"

    def test_related_catch_all_matches_first_degree(self):
        res = self._result_with_degree(100, 25)  # coef 0.5 -> first-degree
        assert res.degree == DEGREE_FIRST
        v = evaluate_expected(res, "related")
        assert v.status == "PASS"

    def test_related_catch_all_on_unrelated_reviews_not_fails(self):
        # "related" expected but none detected: worth a look, not a hard fail.
        res = self._result_with_degree(40, 20)  # unrelated
        v = evaluate_expected(res, "related")
        assert v.status == "REVIEW"

    def test_third_degree_declared_unrelated_detected_reviews(self):
        # Cousin (third-degree) declared but estimate dips to unrelated: the
        # third <-> unrelated boundary is noisy, so REVIEW rather than FAIL.
        res = self._result_with_degree(40, 20)  # unrelated
        v = evaluate_expected(res, "third-degree")
        assert v.status == "REVIEW"

    def test_second_degree_declared_unrelated_detected_fails(self):
        # A close relationship (second-degree) declared but none detected is the
        # swap/mislabel signal, so it stays a FAIL.
        res = self._result_with_degree(40, 20)  # unrelated
        v = evaluate_expected(res, "second-degree")
        assert v.status == "FAIL"

    def test_undetermined_coefficient_reviews(self):
        host = [_md(1, (0, 1)), _md(2, (0, 1))]
        donor = [_md(1, (0, 1)), _md(2, (0, 1))]
        res = relatedness_coefficient(host, donor, "host", "donor")
        assert res.coefficient is None
        v = evaluate_expected(res, "first-degree")
        assert v.status == "REVIEW"


# ---------------------------------------------------------------------------
# Admixture consensus-homozygote swap check
# ---------------------------------------------------------------------------


class TestAdmixConsistency:
    def _consensus_host_donor(self, n: int):
        """n consensus-hom-ref markers shared by host and one donor."""
        host = [_md(i, (0, 0)) for i in range(n)]
        donor = [_md(i, (0, 0)) for i in range(n)]
        return host, [donor]

    def test_clean_mixture_high_swap_pval(self):
        n = 30
        host, donors = self._consensus_host_donor(n)
        # Admix shows the hom-ref consensus with only error-level ALT reads.
        admix = [_md(i, (0, 0), ad_ref=995, ad_alt=5, dp=1000) for i in range(n)]
        res = admix_consistency(host, donors, admix, error_rate=0.01)
        assert res.n_consensus_hom == n
        assert res.n_discordant == 0
        assert res.swap_pval > 0.5

    def test_third_genome_low_swap_pval(self):
        n = 30
        host, donors = self._consensus_host_donor(n)
        # Admix carries the absent ALT allele at ~40% across consensus sites.
        admix = [_md(i, (0, 1), ad_ref=600, ad_alt=400, dp=1000) for i in range(n)]
        res = admix_consistency(host, donors, admix, error_rate=0.01)
        assert res.n_consensus_hom == n
        assert res.n_discordant == n
        assert res.discordant_fraction == 1.0
        assert res.swap_pval < 1e-6

    def test_hom_alt_consensus_uses_ref_as_minor(self):
        n = MIN_CONSENSUS
        host = [_md(i, (1, 1)) for i in range(n)]
        donor = [_md(i, (1, 1)) for i in range(n)]
        # Excess REF reads at a hom-alt consensus signal a third genome.
        admix = [_md(i, (0, 1), ad_ref=400, ad_alt=600, dp=1000) for i in range(n)]
        res = admix_consistency(host, [donor], admix, error_rate=0.01)
        assert res.n_consensus_hom == n
        assert res.n_discordant == n

    def test_informative_sites_not_counted(self):
        # Host hom-ref, donor hom-alt -> informative, never a consensus site.
        host = [_md(i, (0, 0)) for i in range(20)]
        donor = [_md(i, (1, 1)) for i in range(20)]
        admix = [_md(i, (0, 1), ad_ref=500, ad_alt=500, dp=1000) for i in range(20)]
        res = admix_consistency(host, [donor], admix, error_rate=0.01)
        assert res.n_consensus_hom == 0
        assert res.swap_pval == 1.0

    def test_low_depth_sites_skipped(self):
        n = 25
        host, donors = self._consensus_host_donor(n)
        admix = [_md(i, (0, 0), ad_ref=10, ad_alt=0, dp=10) for i in range(n)]
        res = admix_consistency(host, donors, admix, error_rate=0.01, min_dp=100)
        assert res.n_consensus_hom == 0


# ---------------------------------------------------------------------------
# End-to-end: synthetic related VCFs through the real pipeline
# ---------------------------------------------------------------------------


def _write_panel(tmp_path, markers: list[dict], host_name: str, donor_name: str):
    """Write a host VCF and a donor VCF from generate_related_genotypes markers.

    Returns (host_path, donor_path). Both sit at the same marker positions, so
    they line up when re-parsed and blended.
    """
    host_path = tmp_path / f"{host_name}.vcf"
    donor_path = tmp_path / f"{donor_name}.vcf"
    write_genotype_vcf(markers, host_path, host_name, key="host_gt", depth=1000)
    write_genotype_vcf(markers, donor_path, donor_name, key="donor_gt", depth=1000)
    return host_path, donor_path


def _analyse_related(tmp_path, relatedness: str, *, seed: int = 7, n: int = 300, **kwargs):
    """Generate a related host/donor pair, blend an admix, run analyse_sample.

    The admix is a clean host+donor mixture (donor_fraction 0.3). Returns the
    SampleAnalysis. Extra kwargs are forwarded to analyse_sample (e.g.
    expected_relatedness).
    """
    rng = random.Random(seed)
    markers = generate_related_genotypes(n, relatedness, rng, maf_range=(0.2, 0.5))
    host_path, donor_path = _write_panel(tmp_path, markers, "HOST", "DONOR")
    blended = blend_vcfs(
        host_path, donor_path, donor_fraction=0.3, target_depth=1000,
        sample_name="ADMIX", seed=1, error_rate=0.01,
    )
    admix_path = tmp_path / "admix.vcf"
    write_vcf(blended, admix_path)

    host = parse_vcf(host_path, sample="HOST", min_gq=0)
    donor = parse_vcf(donor_path, sample="DONOR", min_gq=0)
    admix = parse_vcf(admix_path, sample="ADMIX", min_dp=0)
    return analyse_sample(
        host, [donor], admix, min_dp=50, min_gq=0, error_rate=0.01, **kwargs
    )


class TestEndToEndSyntheticRelated:
    """Synthetic related VCFs parsed and run through analyse_sample."""

    def test_sibling_pair_detected_first_degree(self, tmp_path):
        analysis = _analyse_related(tmp_path, "sibling", seed=7)
        rel = analysis.result.relatedness[0]
        assert rel.a_name == "host" and rel.b_name == "donor"
        assert rel.degree == DEGREE_FIRST
        assert rel.confidence == "high"

    def test_unrelated_pair_detected_unrelated(self, tmp_path):
        analysis = _analyse_related(tmp_path, "unrelated", seed=11)
        rel = analysis.result.relatedness[0]
        assert rel.degree == DEGREE_UNRELATED

    def test_declared_boundary_mismatch_fails(self, tmp_path):
        # An unrelated pair declared first-degree crosses the boundary -> FAIL.
        analysis = _analyse_related(
            tmp_path, "unrelated", seed=11, expected_relatedness=["first-degree"]
        )
        assert analysis.qc.status == "FAIL"
        assert any("relatedness check FAIL" in w for w in analysis.qc.warnings)

    def test_correct_declaration_passes(self, tmp_path):
        analysis = _analyse_related(
            tmp_path, "sibling", seed=7, expected_relatedness=["first-degree"]
        )
        # No relatedness FAIL/REVIEW contributed by the identity check.
        assert not any("relatedness check FAIL" in w for w in analysis.qc.warnings)
        assert any("relatedness check PASS" in w for w in analysis.qc.warnings)

    def test_clean_mixture_no_swap_flag(self, tmp_path):
        analysis = _analyse_related(tmp_path, "unrelated", seed=11)
        ac = analysis.result.admix_consistency
        assert ac.n_consensus_hom >= MIN_CONSENSUS
        assert ac.n_discordant == 0
        assert ac.swap_pval > 0.5
        assert not any("sample swap" in w for w in analysis.qc.warnings)

    def test_duplicate_host_donor_flagged(self, tmp_path):
        # Host and donor written from the same genotypes (sample reuse). Flagged
        # as identical and FAILed even with no declared expectation.
        rng = random.Random(5)
        markers = generate_related_genotypes(300, "unrelated", rng)
        host_path = tmp_path / "HOST.vcf"
        donor_path = tmp_path / "DONOR.vcf"
        write_genotype_vcf(markers, host_path, "HOST", key="host_gt", depth=1000)
        write_genotype_vcf(markers, donor_path, "DONOR", key="host_gt", depth=1000)
        blended = blend_vcfs(
            host_path, donor_path, donor_fraction=0.3, target_depth=1000,
            sample_name="ADMIX", seed=1,
        )
        admix_path = tmp_path / "admix.vcf"
        write_vcf(blended, admix_path)
        host = parse_vcf(host_path, sample="HOST", min_gq=0)
        donor = parse_vcf(donor_path, sample="DONOR", min_gq=0)
        admix = parse_vcf(admix_path, sample="ADMIX", min_dp=0)
        analysis = analyse_sample(host, [donor], admix, min_dp=50, min_gq=0, error_rate=0.01)

        assert analysis.result.relatedness[0].degree == DEGREE_IDENTICAL
        assert analysis.qc.status == "FAIL"
        assert any("Identical reference samples" in w for w in analysis.qc.warnings)

    def test_cousin_declared_related_never_fails(self, tmp_path):
        # Cousins (true kinship ~0.125) sit right on the third-degree/unrelated
        # boundary, so the estimate scatters across both bands. A correct
        # "related" declaration must never hard-FAIL on that sampling noise; at
        # worst it is REVIEW.
        for seed in range(8):
            analysis = _analyse_related(
                tmp_path, "cousin", seed=seed, n=300,
                expected_relatedness=["related"],
            )
            fails = [w for w in analysis.qc.warnings if "relatedness check FAIL" in w]
            assert not fails, f"seed {seed} hard-failed a cousin declared related: {fails}"

    def test_cousin_declared_third_degree_never_fails(self, tmp_path):
        # Same guard for the explicit third-degree declaration.
        statuses = set()
        for seed in range(8):
            analysis = _analyse_related(
                tmp_path, "cousin", seed=seed, n=300,
                expected_relatedness=["third-degree"],
            )
            verdicts = [w for w in analysis.qc.warnings if "relatedness check" in w]
            assert not any("FAIL" in w for w in verdicts), f"seed {seed}: {verdicts}"
            statuses.update(w.split()[2] for w in verdicts)  # PASS / REVIEW token
        # Across seeds we should see the lenient outcomes, never FAIL.
        assert statuses <= {"PASS", "REVIEW"}

    def test_third_genome_admix_flags_swap(self, tmp_path):
        # Declared host/donor are one unrelated pair; the admix is blended from a
        # different pair at the same positions, so it carries alleles in neither.
        n = 300
        m_decl = generate_related_genotypes(n, "unrelated", random.Random(11))
        host_path, donor_path = _write_panel(tmp_path, m_decl, "HOST", "DONOR")
        m_wrong = generate_related_genotypes(n, "unrelated", random.Random(999))
        wrong_host, wrong_donor = _write_panel(tmp_path, m_wrong, "WH", "WD")
        blended = blend_vcfs(
            wrong_host, wrong_donor, donor_fraction=0.4, target_depth=1000,
            sample_name="ADMIX", seed=3, error_rate=0.01,
        )
        admix_path = tmp_path / "swap_admix.vcf"
        write_vcf(blended, admix_path)

        host = parse_vcf(host_path, sample="HOST", min_gq=0)
        donor = parse_vcf(donor_path, sample="DONOR", min_gq=0)
        admix = parse_vcf(admix_path, sample="ADMIX", min_dp=0)
        analysis = analyse_sample(host, [donor], admix, min_dp=50, min_gq=0, error_rate=0.01)

        ac = analysis.result.admix_consistency
        assert ac.n_consensus_hom >= MIN_CONSENSUS
        assert ac.discordant_fraction > 0.5
        assert ac.swap_pval < 1e-6
        assert analysis.qc.status == "REVIEW"
        assert any("sample swap" in w for w in analysis.qc.warnings)


# ---------------------------------------------------------------------------
# CLI input validation
# ---------------------------------------------------------------------------


class TestCliRejectsIdentical:
    """--expected-relatedness must refuse 'identical' with a clear reason."""

    def _base_argv(self, rel: str) -> list[str]:
        return [
            "monitor",
            "--panel-vcf", "x.vcf",
            "--admix-vcf", "x.vcf",
            "--host-sample", "H",
            "--donor-sample", "D",
            "--sample", "S",
            "--expected-relatedness", rel,
        ]

    def test_identical_is_rejected_at_parse(self, capsys):
        # Errors during argument parsing, before any VCF is opened.
        with pytest.raises(SystemExit):
            main(self._base_argv("identical"))
        err = capsys.readouterr().err.lower()
        assert "identical" in err
        assert "syngeneic" in err or "twin" in err

    def test_garbage_value_is_rejected(self, capsys):
        with pytest.raises(SystemExit):
            main(self._base_argv("cousin-ish"))
        err = capsys.readouterr().err.lower()
        assert "invalid expected relatedness" in err

    def test_count_mismatch_errors(self):
        # Two donors but one declaration: clear error, before any VCF is opened.
        argv = [
            "monitor",
            "--panel-vcf", "x.vcf",
            "--admix-vcf", "x.vcf",
            "--host-sample", "H",
            "--donor-sample", "D1",
            "--donor-sample", "D2",
            "--sample", "S",
            "--expected-relatedness", "unrelated",
        ]
        with pytest.raises(SystemExit) as exc:
            main(argv)
        assert "donor" in str(exc.value).lower()
