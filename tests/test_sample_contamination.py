"""Tests for allomix.qc.sample_contamination — in-data third-party contamination estimate."""

from allomix.genotype import MarkerData
from allomix.qc.sample_contamination import estimate_contamination


def _md(
    pos: int,
    gt: tuple[int, int],
    ad_ref: int = 0,
    ad_alt: int = 0,
    dp: int | None = None,
    chrom: str = "chr1",
) -> MarkerData:
    """MarkerData with the genotype and read counts a test needs."""
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        gt=gt,
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=ad_ref + ad_alt if dp is None else dp,
    )


def _consensus_set(minor_fracs: list[float], dp: int = 2000, hom_alt: bool = False):
    """Build host/donor/admix lists of consensus-hom markers.

    Host and every donor are homozygous for the same allele (hom-ref unless
    ``hom_alt``); each admix marker carries ``minor_fracs[i]`` of the minor
    (absent) allele. Returns ``(host, [donor], admix)``.
    """
    host, donor, admix = [], [], []
    gt = (1, 1) if hom_alt else (0, 0)
    for i, frac in enumerate(minor_fracs):
        pos = (i + 1) * 100
        host.append(_md(pos, gt, ad_ref=0 if hom_alt else dp, ad_alt=dp if hom_alt else 0, dp=dp))
        donor.append(_md(pos, gt, ad_ref=0 if hom_alt else dp, ad_alt=dp if hom_alt else 0, dp=dp))
        minor = round(frac * dp)
        if hom_alt:
            # minor allele is REF
            admix.append(_md(pos, gt, ad_ref=minor, ad_alt=dp - minor, dp=dp))
        else:
            # minor allele is ALT
            admix.append(_md(pos, gt, ad_ref=dp - minor, ad_alt=minor, dp=dp))
    return host, [donor], admix


class TestCleanSample:
    """No contamination: the minor fraction sits at the error floor everywhere."""

    def test_zero_minor_reads(self):
        host, donors, admix = _consensus_set([0.0] * 50)
        res = estimate_contamination(host, donors, admix)
        assert res.n_markers == 50
        assert res.contamination_fraction == 0.0
        assert res.p_value == 1.0  # Y == 0

    def test_uniform_low_error_not_called_contamination(self):
        # Every site at a flat 0.05% — uniform error, not heterogeneous leakage.
        host, donors, admix = _consensus_set([0.0005] * 50)
        res = estimate_contamination(host, donors, admix)
        # median ~ floor (both at 0.0005), so the excess is ~0.
        assert res.contamination_fraction < 1e-3


class TestContaminatedSample:
    """A heterogeneous minor-allele signal above the error floor."""

    def test_heterogeneous_signal_detected(self):
        # 40 clean error-only sites (p25 lands here), then 60 carrier sites at
        # 0.6% (the median lands here): a clear heterogeneous excess.
        fracs = [0.0] * 40 + [0.006] * 60
        host, donors, admix = _consensus_set(fracs)
        res = estimate_contamination(host, donors, admix)
        assert res.floor_empirical is True
        assert res.error_floor < 1e-9  # p25 is a clean site
        assert res.contamination_fraction > 0.004  # median(0.006) - floor(~0)
        assert res.p_value < 1e-6

    def test_hom_alt_consensus_minor_is_ref(self):
        # Same signal but the consensus genotype is hom-alt, so the minor allele
        # is REF. The estimate must be identical to the hom-ref case.
        fracs = [0.0] * 40 + [0.006] * 60
        h_ref, d_ref, a_ref = _consensus_set(fracs, hom_alt=False)
        h_alt, d_alt, a_alt = _consensus_set(fracs, hom_alt=True)
        r_ref = estimate_contamination(h_ref, d_ref, a_ref)
        r_alt = estimate_contamination(h_alt, d_alt, a_alt)
        assert r_alt.contamination_fraction == r_ref.contamination_fraction
        assert r_alt.n_markers == r_ref.n_markers


class TestHighFractionExclusion:
    """Genotype-miscall / swap sites above the cap must not inflate the estimate."""

    def test_high_sites_excluded_and_counted(self):
        # 40 clean + 10 at 50% (miscalls). Without the cap the median would jump.
        fracs = [0.0005] * 40 + [0.5] * 10
        host, donors, admix = _consensus_set(fracs)
        res = estimate_contamination(host, donors, admix)
        assert res.n_excluded_high == 10
        assert res.n_markers == 40
        assert res.contamination_fraction < 1e-3  # the 50% sites are gone

    def test_cap_is_configurable(self):
        fracs = [0.0005] * 40 + [0.2] * 10
        host, donors, admix = _consensus_set(fracs)
        # Default cap 0.10 drops the 0.2 sites...
        assert estimate_contamination(host, donors, admix).n_excluded_high == 10
        # ...a higher cap keeps them.
        res = estimate_contamination(host, donors, admix, max_site_frac=0.5)
        assert res.n_excluded_high == 0
        assert res.n_markers == 50


class TestMarkerSelection:
    """Only consensus-homozygous, adequately-covered autosomal markers qualify."""

    def test_no_consensus_markers(self):
        # Host hom-ref, donor hom-alt everywhere: informative, never consensus.
        host = [_md(i * 100, (0, 0), ad_ref=2000, dp=2000) for i in range(30)]
        donor = [_md(i * 100, (1, 1), ad_alt=2000, dp=2000) for i in range(30)]
        admix = [_md(i * 100, (0, 0), ad_ref=1000, ad_alt=1000, dp=2000) for i in range(30)]
        res = estimate_contamination(host, [donor], admix)
        assert res.n_markers == 0
        assert res.contamination_fraction == 0.0
        assert res.p_value == 1.0

    def test_min_dp_filter(self):
        host, donors, admix = _consensus_set([0.0] * 30, dp=50)
        # All admix dp = 50; require 100 -> nothing qualifies.
        res = estimate_contamination(host, donors, admix, min_dp=100)
        assert res.n_markers == 0

    def test_het_consensus_excluded(self):
        # Host and donor both het: non-informative but NOT homozygous, so it is
        # not a usable consensus-hom marker (both carry both alleles).
        host = [_md(i * 100, (0, 1), ad_ref=1000, ad_alt=1000, dp=2000) for i in range(30)]
        donor = [_md(i * 100, (0, 1), ad_ref=1000, ad_alt=1000, dp=2000) for i in range(30)]
        admix = [_md(i * 100, (0, 1), ad_ref=1000, ad_alt=1000, dp=2000) for i in range(30)]
        res = estimate_contamination(host, [donor], admix)
        assert res.n_markers == 0

    def test_sex_chroms_skipped(self):
        host, donors, admix = _consensus_set([0.0] * 30)
        for m in host + donors[0] + admix:
            m.chrom = "chrX"
        res = estimate_contamination(host, donors, admix)
        assert res.n_markers == 0


class TestErrorFloorFallback:
    """Too few markers for a stable percentile falls back to the error rate."""

    def test_small_panel_uses_error_rate_floor(self):
        host, donors, admix = _consensus_set([0.0] * 5)
        res = estimate_contamination(host, donors, admix, error_rate=0.012)
        assert res.floor_empirical is False
        # Fallback floor is the per-direction error rate (error_rate / 3).
        assert abs(res.error_floor - 0.012 / 3) < 1e-6

    def test_large_panel_uses_empirical_floor(self):
        host, donors, admix = _consensus_set([0.0] * 50)
        res = estimate_contamination(host, donors, admix)
        assert res.floor_empirical is True


class TestMultiDonor:
    """Every donor must share the consensus homozygote."""

    def test_one_discordant_donor_drops_marker(self):
        host = [_md(i * 100, (0, 0), ad_ref=2000, dp=2000) for i in range(30)]
        d1 = [_md(i * 100, (0, 0), ad_ref=2000, dp=2000) for i in range(30)]
        # d2 hom-alt at every site -> no consensus.
        d2 = [_md(i * 100, (1, 1), ad_alt=2000, dp=2000) for i in range(30)]
        admix = [_md(i * 100, (0, 0), ad_ref=1980, ad_alt=20, dp=2000) for i in range(30)]
        res = estimate_contamination(host, [d1, d2], admix)
        assert res.n_markers == 0
