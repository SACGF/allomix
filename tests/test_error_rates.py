"""Tests for allomix.calibration.error_rates — per-site empirical error rate estimation."""

import tempfile
from pathlib import Path

import pytest

from allomix.calibration.error_rates import (
    DEFAULT_ERROR_FLOOR,
    MarkerError,
    MarkerErrorRates,
    errors_to_rates,
    estimate_error_rates,
    load_error_table,
    save_error_table,
)
from allomix.cli import main
from allomix.estimate.likelihood import (
    PanelCalibration,
    log_likelihood_marker_bb,
    total_log_likelihood_bb,
)
from allomix.genotype import InformativeMarker, MarkerData

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _homref(pos: int, depth: int, n_alt: int, chrom: str = "chr1") -> MarkerData:
    """A hom-ref MarkerData with the given depth and observed ALT count."""
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="T",
        gt=(0, 0),
        ad_ref=depth - n_alt,
        ad_alt=n_alt,
        dp=depth,
    )


def _homalt(pos: int, depth: int, n_ref: int, chrom: str = "chr1") -> MarkerData:
    """A hom-alt MarkerData with the given depth and observed REF count."""
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="T",
        gt=(1, 1),
        ad_ref=n_ref,
        ad_alt=depth - n_ref,
        dp=depth,
    )


def _het(pos: int, depth: int, n_alt: int, chrom: str = "chr1") -> MarkerData:
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="T",
        gt=(0, 1),
        ad_ref=depth - n_alt,
        ad_alt=n_alt,
        dp=depth,
    )


def _informative(
    pos: int,
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    ad_ref: int,
    ad_alt: int,
    chrom: str = "chr1",
) -> InformativeMarker:
    return InformativeMarker(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="T",
        host_gt=host_gt,
        donor_gts=[donor_gt],
        marker_type=0,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
    )


# ---------------------------------------------------------------------------
# estimate_error_rates
# ---------------------------------------------------------------------------


class TestEstimateErrorRates:
    def test_max_vaf_homref_drops_miscalled_het(self) -> None:
        """A hom-ref observation with vaf > max_vaf_homref is excluded entirely."""
        # Clean baseline + one obviously-miscalled-het sample at the same site.
        clean = [[_homref(100, 2000, 5)] for _ in range(10)]  # 5/2000 ~ 0.0025
        miscalled_het = [[_homref(100, 2000, 800)]]  # vaf 0.4 -- a het, not hom-ref
        out = estimate_error_rates(clean + miscalled_het, min_reads=1000)
        key = ("chr1", 100, "A", "T")
        est = out[key].e_refalt
        # Without the filter the rate would be inflated by the 800 ALT-read
        # outlier: (50 + 800) / (11 * 2000) ~ 0.039. With the filter it stays
        # near the clean baseline 50 / (10 * 2000) = 0.0025.
        assert est is not None and est < 0.005

    def test_min_vaf_homalt_drops_miscalled_het(self) -> None:
        """Same as above for the hom-alt direction."""
        clean = [[_homalt(100, 2000, 5)] for _ in range(10)]
        miscalled_het = [[_homalt(100, 2000, 800)]]  # vaf 0.6 != hom-alt
        out = estimate_error_rates(clean + miscalled_het, min_reads=1000)
        key = ("chr1", 100, "A", "T")
        est = out[key].e_altref
        assert est is not None and est < 0.005

    def test_min_reads_filter_per_direction(self) -> None:
        """A site below min_reads in one direction returns None for that direction."""
        # Site A: enough hom-ref reads (>1000), no hom-alt
        # Site B: enough hom-alt reads (>1000), no hom-ref
        samples = [
            [_homref(100, 2000, 6), _homalt(200, 500, 1)],  # 500 < 1000 for site B
        ]
        out = estimate_error_rates(samples, min_reads=1000)
        a = out[("chr1", 100, "A", "T")]
        assert a.e_refalt is not None
        assert a.e_altref is None  # no hom-alt obs
        assert ("chr1", 200, "A", "T") not in out  # both directions below floor

    def test_no_usable_observations_omits_site(self) -> None:
        """A site with only-het observations is omitted from the output."""
        samples = [[_het(100, 1000, 500)] for _ in range(50)]
        out = estimate_error_rates(samples, min_reads=1)
        assert out == {}

    def test_pooling_weights_by_depth(self) -> None:
        """Two samples at the same site pool proportionally to depth.

        Sample with depth 100, 1 ALT (rate 0.01) and sample with depth 10000,
        10 ALT (rate 0.001) should pool to (1+10)/(100+10000) ~ 0.00109, much
        closer to the deep sample's rate than the shallow one's.
        """
        samples = [
            [_homref(100, 100, 1)],
            [_homref(100, 10000, 10)],
        ]
        out = estimate_error_rates(samples, min_reads=100)
        est = out[("chr1", 100, "A", "T")].e_refalt
        assert est is not None
        # Read-pooled MLE: 11 / 10100 ~ 0.001089
        assert abs(est - 11.0 / 10100.0) < 1e-9


# ---------------------------------------------------------------------------
# save / load round-trip
# ---------------------------------------------------------------------------


class TestSaveLoadRoundtrip:
    def test_roundtrip_with_na(self) -> None:
        """NA entries round-trip to None; floor is applied at load."""
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1",
                100,
                "A",
                "T",
                e_refalt=1e-3,
                e_altref=None,
                n_reads_homref=5000,
                n_reads_homalt=0,
            ),
            ("chr1", 200, "C", "G"): MarkerError(
                "chr1",
                200,
                "C",
                "G",
                e_refalt=0.0,
                e_altref=2e-4,  # zero rate hits the floor
                n_reads_homref=5000,
                n_reads_homalt=5000,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "errors.tsv"
            save_error_table(errors, path)
            loaded = load_error_table(path)

        assert loaded[("chr1", 100, "A", "T")] == MarkerErrorRates(1e-3, None)
        # Floor applied to the 0.0 entry.
        rates = loaded[("chr1", 200, "C", "G")]
        assert rates.e_refalt == DEFAULT_ERROR_FLOOR
        assert rates.e_altref == 2e-4

    def test_loader_disable_floor(self) -> None:
        """Setting error_floor=0 disables the floor."""
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1",
                100,
                "A",
                "T",
                e_refalt=0.0,
                e_altref=0.0,
                n_reads_homref=5000,
                n_reads_homalt=5000,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "errors.tsv"
            save_error_table(errors, path)
            loaded = load_error_table(path, error_floor=0.0)
        assert loaded[("chr1", 100, "A", "T")] == MarkerErrorRates(0.0, 0.0)


# ---------------------------------------------------------------------------
# errors_to_rates
# ---------------------------------------------------------------------------


class TestErrorsToRates:
    def test_floor_applied_in_helper(self) -> None:
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1",
                100,
                "A",
                "T",
                e_refalt=1e-7,  # below default floor
                e_altref=1e-3,
                n_reads_homref=5000,
                n_reads_homalt=5000,
            ),
        }
        rates = errors_to_rates(errors)
        assert rates[("chr1", 100, "A", "T")] == MarkerErrorRates(DEFAULT_ERROR_FLOOR, 1e-3)


# ---------------------------------------------------------------------------
# Likelihood integration
# ---------------------------------------------------------------------------


class TestLikelihoodIntegration:
    def test_symmetric_per_marker_matches_global(self) -> None:
        """Per-marker (e/3, e/3) matches the legacy 4-state path.

        The 4-state model renormalises onto {REF, ALT}, adding an O(e^2) gap in
        general. At w=0.5 both paths collapse to p_alt=0.5 by symmetry, so the
        gap vanishes exactly; a tight bound guards that identity (the loose 0.1
        tolerance would have hidden a real regression).
        """
        m = _informative(100, host_gt=(0, 0), donor_gt=(1, 1), ad_ref=900, ad_alt=100)
        ll_symmetric = log_likelihood_marker_bb(
            m.admix_ad_ref,
            m.admix_ad_alt,
            w=0.5,
            error_rate=0.01,
            rho=100.0,
        )
        ll_asymmetric = log_likelihood_marker_bb(
            m.admix_ad_ref,
            m.admix_ad_alt,
            w=0.5,
            error_rate=0.01,
            rho=100.0,
            e_refalt=0.01 / 3.0,
            e_altref=0.01 / 3.0,
        )
        assert ll_symmetric == pytest.approx(ll_asymmetric, abs=1e-9)

    def test_asymmetric_shifts_likelihood(self) -> None:
        """Switching to a much larger e_refalt at a hom-ref-like marker (w=1)
        increases the likelihood of seeing ALT reads.
        """
        # ad_alt = 10 at depth 1000 is unlikely under e_refalt = 1e-4 but
        # easy under e_refalt = 1e-2.
        ll_low = log_likelihood_marker_bb(
            ad_ref=990,
            ad_alt=10,
            w=1.0,
            rho=100.0,
            e_refalt=1e-4,
            e_altref=1e-4,
        )
        ll_high = log_likelihood_marker_bb(
            ad_ref=990,
            ad_alt=10,
            w=1.0,
            rho=100.0,
            e_refalt=1e-2,
            e_altref=1e-4,
        )
        assert ll_high > ll_low

    def test_marker_errors_missing_falls_back(self) -> None:
        """A marker absent from the error table uses the global 4-state path
        (i.e. matches the same call with no calibration).
        """
        m = _informative(100, host_gt=(0, 1), donor_gt=(1, 1), ad_ref=400, ad_alt=600)
        ll_with_empty = total_log_likelihood_bb(
            [m],
            f_donor=0.5,
            error_rate=0.01,
            rho=100.0,
            calibration=PanelCalibration(),
        )
        ll_no_table = total_log_likelihood_bb(
            [m],
            f_donor=0.5,
            error_rate=0.01,
            rho=100.0,
        )
        assert ll_with_empty == pytest.approx(ll_no_table)

    def test_marker_errors_one_direction_missing_falls_back(self) -> None:
        """If only one per-direction rate is known for a marker, the marker
        falls back to the symmetric path (cannot specify p_alt fully without
        both).
        """
        m = _informative(100, host_gt=(0, 1), donor_gt=(1, 1), ad_ref=400, ad_alt=600)
        partial = {("chr1", 100, "A", "T"): MarkerErrorRates(1e-3, None)}
        ll_partial = total_log_likelihood_bb(
            [m],
            f_donor=0.5,
            error_rate=0.01,
            rho=100.0,
            calibration=PanelCalibration(errors=partial),
        )
        ll_baseline = total_log_likelihood_bb(
            [m],
            f_donor=0.5,
            error_rate=0.01,
            rho=100.0,
        )
        assert ll_partial == pytest.approx(ll_baseline)


# ---------------------------------------------------------------------------
# estimate-errors CLI: hom-ref background VCF folding (--homref-vcf)
# ---------------------------------------------------------------------------

_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def _write_vcf(path: Path, rows: list[str]) -> None:
    path.write_text(_VCF_HEADER + "".join(r + "\n" for r in rows))


class TestEstimateErrorsHomrefVcf:
    """`allomix estimate-errors --homref-vcf` folds hom-ref background sites in."""

    def _panel_vcf(self, path: Path) -> None:
        # One both-hom-alt panel site: REF reads there measure alt->ref.
        # No hom-ref-with-reads site, so ref->alt is unmeasurable from this VCF
        # alone (mirrors the variant-only joint call).
        _write_vcf(
            path,
            ["chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:DP\t1/1:3,1497:1500\t1/1:2,1498:1500"],
        )

    def _homref_vcf(self, path: Path) -> None:
        # A forced hom-ref background site (e.g. an amplicon midpoint): the few
        # ALT reads at a 0/0 call measure ref->alt.
        _write_vcf(
            path,
            ["chr1\t200\t.\tC\tT\t.\t.\t.\tGT:AD:DP\t0/0:1497,3:1500\t0/0:1498,2:1500"],
        )

    def test_panel_alone_has_no_refalt(self, tmp_path):
        panel = tmp_path / "panel.vcf"
        out = tmp_path / "err.tsv"
        self._panel_vcf(panel)
        main(
            [
                "estimate-errors",
                str(panel),
                "--sample",
                "S1",
                "--sample",
                "S2",
                "--min-reads",
                "100",
                "-o",
                str(out),
            ]
        )
        table = load_error_table(out, error_floor=0.0)
        # alt->ref recovered at the hom-alt site; no ref->alt anywhere.
        assert (("chr1", 100, "A", "G")) in table
        assert table[("chr1", 100, "A", "G")].e_altref is not None
        assert all(r.e_refalt is None for r in table.values())

    def test_homref_vcf_adds_refalt_direction(self, tmp_path):
        panel = tmp_path / "panel.vcf"
        homref = tmp_path / "homref.vcf"
        out = tmp_path / "err.tsv"
        self._panel_vcf(panel)
        self._homref_vcf(homref)
        main(
            [
                "estimate-errors",
                str(panel),
                "--homref-vcf",
                str(homref),
                "--sample",
                "S1",
                "--sample",
                "S2",
                "--min-reads",
                "100",
                "-o",
                str(out),
            ]
        )
        table = load_error_table(out, error_floor=0.0)
        # Background site now supplies the ref->alt rate: 5 ALT / 3000 DP.
        homref_entry = table[("chr1", 200, "C", "T")]
        assert homref_entry.e_refalt == pytest.approx(5 / 3000, rel=1e-6)
        assert homref_entry.e_altref is None
        # The panel hom-alt site's alt->ref is still present and unchanged.
        assert table[("chr1", 100, "A", "G")].e_altref is not None

    def test_homref_vcf_requires_samples(self, tmp_path):
        homref = tmp_path / "homref.vcf"
        self._homref_vcf(homref)
        with pytest.raises(SystemExit):
            main(
                [
                    "estimate-errors",
                    str(homref),
                    "--homref-vcf",
                    str(homref),
                    "-o",
                    str(tmp_path / "err.tsv"),
                ]
            )
