"""In-data contamination estimate at consensus-homozygous markers.

Measures the third-party (neither-host-nor-donor) read signal directly from the
admixture VCF, independent of sequencing-run metadata. Where host and every donor
are homozygous for the *same* allele the admixture is a mixture of identical
homozygotes, so the other (minor) allele cannot come from either contributor: its
reads are sequencing error plus any DNA leaking in from co-loaded samples (index
hopping, or physical / library cross-contamination). The excess of the
minor-allele fraction over the per-site sequencing-error floor estimates that
contamination.

This is the "is contamination actually elevated" measurement of issue #12,
separate from the index-hopping provenance flag (shared sequencing run between
host and admix), which is metadata the joint-calling pipeline carries and allomix
never derives from BAMs. Complementary: this says whether a third genome is in the
reads; the run flag says whether index hopping is a plausible mechanism.

Relationship to the existing checks:

  - ``allomix.relatedness.admix_consistency`` is a gross-swap detector on the same
    marker set: it counts sites where the minor allele is *individually*
    significant (a whole third genome near 50%). It does not fire on a ~0.2% floor
    spread across every site, which this estimator targets.
  - ``allomix.detect.host_presence_test`` works on donor-homozygous markers where
    the *host* carries the minor allele, so it measures host, not third parties.

The headline estimate is the background-subtracted *median* per-site minor
fraction, not a read-weighted pooled mean: on real panel data a few
genotype-miscall or mapping-artifact sites sit at 40-100% minor allele and would
dominate a pooled mean (validated on SRP434573, issue #12). Those are dropped by
an upper fraction cap (``max_site_frac``) and the median is robust to any that
remain. A gross swap puts many sites above the cap and is left to
``admix_consistency``; this estimator stays on the low-level floor.
"""

import statistics
from dataclasses import dataclass

from scipy.stats import poisson

from allomix.constants import DEFAULT_ERROR_RATE, N_OTHER_BASES
from allomix.detect import ErrorRateSource
from allomix.error_rates import MarkerErrorRates
from allomix.genotype import MarkerData, MarkerKey, is_sex_chrom, marker_key

# Consensus-hom markers whose admix minor-allele fraction exceeds this are not
# low-level contamination: genotype miscalls, mapping artifacts, or (en masse) a
# gross swap. Dropped from the estimate and tallied in ``n_excluded_high``. Well
# above the realistic contamination range (SRP434573 ~0.2% typical, ~1.5% at p95)
# and well below a miscall (~50%).
DEFAULT_MAX_SITE_FRAC = 0.10

# Sequencing-error floor read off the data: at consensus-hom sites with no
# co-loaded carrier of the minor allele the only minor reads are error, so a low
# percentile of per-site minor fractions estimates the floor without trusting a
# global rate. Contamination is the excess of the median over this floor; uniform
# error lifts the floor too and is correctly not called contamination. The 10th
# percentile (not a quartile, which would already sit in the contaminated range on
# a dense pool) holds on no-carrier sites even when carriers dominate. Needs
# enough markers; below that the per-site/global error rate is used.
CONTAMINATION_FLOOR_PCTL = 0.10
MIN_MARKERS_FOR_EMPIRICAL_FLOOR = 20


@dataclass
class ContaminationResult:
    """In-data third-party contamination estimate at consensus-hom markers.

    Attributes:
        n_markers: Consensus-homozygous markers used (after the fraction cap).
        contamination_fraction: Headline estimate ``max(0, median_minor_frac -
            error_floor)``: typical per-site minor-allele fraction above the
            sequencing-error background.
        median_minor_frac: Median per-site minor fraction (raw, pre-subtraction).
        error_floor: Sequencing-error background subtracted from the median. The
            low percentile of per-site minor fractions when there are enough
            markers, else the per-site / global error rate. See ``floor_empirical``.
        floor_empirical: True when ``error_floor`` is the data percentile, False
            when it fell back to the supplied error rate.
        pooled_minor_frac: Read-weighted pooled minor fraction. Diagnostic only;
            inflated by residual high-fraction sites, hence the median headline.
        n_minor_reads: Pooled minor-allele reads ``Y`` over the used markers.
        total_depth: Pooled depth ``N`` over the used markers.
        p_value: One-sided pooled-Poisson presence-test p-value for the minor
            reads exceeding the error-floor background. Significant for any real
            excess at high depth; use ``contamination_fraction`` for clinical
            magnitude.
        n_excluded_high: Consensus-hom markers dropped above ``max_site_frac``.
        used_per_site_error: True when at least one marker used a per-site rate.
        error_rate_source: "per-site", "global-fallback", "mixed", or "none", as
            in ``allomix.detect``. The rates used for the test background, not the
            subtracted floor (see ``floor_empirical``).
    """

    n_markers: int
    contamination_fraction: float
    median_minor_frac: float
    error_floor: float
    floor_empirical: bool
    pooled_minor_frac: float
    n_minor_reads: int
    total_depth: int
    p_value: float
    n_excluded_high: int
    used_per_site_error: bool
    error_rate_source: ErrorRateSource


@dataclass
class _ConsensusHomRecord:
    """Per-marker record at a consensus-homozygous site."""

    minor_reads: int
    dp: int
    e: float  # per-site sequencing-error background for the minor allele


def _select_consensus_hom(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    marker_errors: dict[MarkerKey, MarkerErrorRates] | None,
    fallback_e: float,
    error_floor: float,
    min_dp: int,
) -> tuple[list[_ConsensusHomRecord], int, int]:
    """Collect consensus-homozygous markers with minor-allele counts and rates.

    A marker qualifies when the host and every donor are homozygous for the same
    allele and the admixture has depth ``>= min_dp``. The minor allele is the one
    absent from that homozygote (ALT for a hom-ref consensus, REF for hom-alt).
    Sex / mitochondrial contigs are skipped.

    Returns ``(records, n_per_site, n_fallback)`` where the tallies count how many
    markers drew a per-site error rate from the table vs the global fallback.
    """
    donor_maps = [{marker_key(m): m for m in d} for d in donors]
    admix_map = {marker_key(m): m for m in admix}

    records: list[_ConsensusHomRecord] = []
    n_per_site = 0
    n_fallback = 0
    for mh in host:
        if is_sex_chrom(mh.chrom):
            continue
        dose_h = mh.gt[0] + mh.gt[1]
        if dose_h not in (0, 2):
            continue  # host must be homozygous
        key = marker_key(mh)
        consensus = True
        for dm in donor_maps:
            md = dm.get(key)
            if md is None or (md.gt[0] + md.gt[1]) != dose_h:
                consensus = False
                break
        if not consensus:
            continue
        ma = admix_map.get(key)
        if ma is None or ma.dp < min_dp or ma.dp <= 0:
            continue

        # Minor allele is the one absent from the consensus homozygote, and its
        # per-site error direction follows: hom-ref consensus -> minor is ALT
        # (REF->ALT miscalls), hom-alt -> minor is REF (ALT->REF miscalls).
        entry = marker_errors.get(key) if marker_errors else None
        if dose_h == 0:
            minor_reads = ma.ad_alt
            e_dir = entry.e_refalt if entry is not None else None
        else:
            minor_reads = ma.ad_ref
            e_dir = entry.e_altref if entry is not None else None

        if e_dir is not None:
            e = max(e_dir, error_floor)
            n_per_site += 1
        else:
            e = fallback_e
            n_fallback += 1

        records.append(_ConsensusHomRecord(minor_reads=minor_reads, dp=ma.dp, e=e))
    return records, n_per_site, n_fallback


def estimate_contamination(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    *,
    marker_errors: dict[MarkerKey, MarkerErrorRates] | None = None,
    error_rate: float = DEFAULT_ERROR_RATE,
    error_floor: float = 1e-5,
    min_dp: int = 1,
    max_site_frac: float = DEFAULT_MAX_SITE_FRAC,
) -> ContaminationResult:
    """Estimate third-party contamination from consensus-homozygous markers.

    Args:
        host: Parsed host markers.
        donors: One parsed marker list per donor.
        admix: Parsed admixture markers (raw allele depths).
        marker_errors: Optional per-site, per-direction error table (see
            ``allomix.error_rates.load_error_table``). Each marker uses its
            per-direction rate; markers missing from the table fall back to
            ``error_rate / N_OTHER_BASES``.
        error_rate: Global symmetric sequencing error rate; the per-direction
            fallback is ``error_rate / N_OTHER_BASES``.
        error_floor: Lower bound applied to every per-marker background rate.
        min_dp: Minimum admixture depth for a marker to be used.
        max_site_frac: Drop consensus-hom markers whose admix minor-allele
            fraction exceeds this (genotype miscall / mapping artifact / swap
            site); they are not low-level contamination. See
            ``DEFAULT_MAX_SITE_FRAC``.

    Returns:
        A ``ContaminationResult``. ``n_markers == 0`` (with zero fractions and
        ``p_value`` 1.0) when no consensus-hom markers survive.
    """
    fallback_e = max(error_rate / N_OTHER_BASES, error_floor)
    records, n_per_site, n_fallback = _select_consensus_hom(
        host, donors, admix, marker_errors, fallback_e, error_floor, min_dp
    )

    # Drop high-fraction sites (not low-level contamination) before estimating.
    kept: list[_ConsensusHomRecord] = []
    n_excluded_high = 0
    for r in records:
        if r.dp > 0 and (r.minor_reads / r.dp) > max_site_frac:
            n_excluded_high += 1
        else:
            kept.append(r)

    if n_per_site and n_fallback:
        source: ErrorRateSource = "mixed"
    elif n_per_site:
        source = "per-site"
    elif n_fallback:
        source = "global-fallback"
    else:
        source = "none"

    if not kept:
        return ContaminationResult(
            n_markers=0,
            contamination_fraction=0.0,
            median_minor_frac=0.0,
            error_floor=0.0,
            floor_empirical=False,
            pooled_minor_frac=0.0,
            n_minor_reads=0,
            total_depth=0,
            p_value=1.0,
            n_excluded_high=n_excluded_high,
            used_per_site_error=n_per_site > 0,
            error_rate_source=source,
        )

    site_fracs = sorted(r.minor_reads / r.dp for r in kept)
    median_minor_frac = float(statistics.median(site_fracs))

    # Error floor: empirical low percentile when there are enough markers (the
    # no-carrier / error sites), else the median per-marker error rate. The
    # Poisson background rate below is additionally floored at ``error_floor`` so
    # it stays positive even when the empirical floor is exactly zero.
    if len(site_fracs) >= MIN_MARKERS_FOR_EMPIRICAL_FLOOR:
        idx = min(len(site_fracs) - 1, int(len(site_fracs) * CONTAMINATION_FLOOR_PCTL))
        floor = float(site_fracs[idx])
        floor_empirical = True
    else:
        floor = float(statistics.median([r.e for r in kept]))
        floor_empirical = False
    contamination_fraction = max(0.0, median_minor_frac - floor)

    Y = sum(r.minor_reads for r in kept)
    N = sum(r.dp for r in kept)
    pooled_minor_frac = Y / N if N else 0.0
    # Presence test: pooled minor reads against the error-floor background.
    background_rate = max(floor, error_floor)
    Lam = N * background_rate
    if Lam <= 0.0 or Y == 0:
        p_value = 1.0
    else:
        p_value = float(poisson.sf(Y - 1, Lam))

    return ContaminationResult(
        n_markers=len(kept),
        contamination_fraction=contamination_fraction,
        median_minor_frac=median_minor_frac,
        error_floor=floor,
        floor_empirical=floor_empirical,
        pooled_minor_frac=pooled_minor_frac,
        n_minor_reads=Y,
        total_depth=N,
        p_value=p_value,
        n_excluded_high=n_excluded_high,
        used_per_site_error=n_per_site > 0,
        error_rate_source=source,
    )


__all__ = ["ContaminationResult", "estimate_contamination"]
