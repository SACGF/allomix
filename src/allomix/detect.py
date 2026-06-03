"""Host-presence detection at donor-homozygous markers.

Tests whether a minor contributor (the host, re-occurring post-HSCT) is present
at all, using the markers where every donor is homozygous for the same allele
and the host carries the donor-absent allele. At those markers the
donor-absent allele is expected at the per-site sequencing-error background
in a pure-donor sample, so its read counts give a one-sided detection test
against that background.

This is complementary to the fraction MLE in ``chimerism``: the MLE estimates
the magnitude, this test guards the low end and answers "is host present?"
directly. See ``claude/20_host_presence_detection_plan.md`` for the full
rationale; the calibration evidence under realistic overdispersion lives in
``paper/scripts/run_presence_lod.py``.

The two statistics computed here are:

  - Pooled Poisson:  ``Y = sum y_i``, ``Lam = sum n_i * e_i``,
                     ``p_pois = P(Poisson(Lam) >= Y)``. Transparent and robust.
  - LRT:             ``q_i(f_h) = e_i + (h_i / 2) f_h``; bounded-MLE LRT with
                     a chi-bar-square one-sided p-value
                     (``0.5 * P(chi2_1 >= D)`` for ``D > 0``, else 1).

The LRT also returns ``f_host_mle`` and a profile-likelihood 95% CI.
"""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2, poisson

from allomix.constants import CI_LEVEL, DEFAULT_ERROR_RATE, N_OTHER_BASES
from allomix.error_rates import MarkerErrorRates
from allomix.genotype import InformativeMarker, MarkerKey  # MarkerKey re-exported below

# Direction of the per-site error rate we use at a donor-homozygous marker:
# "ref->alt" picks ``e_refalt`` (donor is hom-ref, donor-absent allele is ALT);
# "alt->ref" picks ``e_altref`` (donor is hom-alt, donor-absent allele is REF).
Direction = Literal["ref->alt", "alt->ref"]

# Tracks where each marker's background rate came from. "none" is reserved for
# the degenerate case of zero usable markers.
ErrorRateSource = Literal["per-site", "global-fallback", "mixed", "none"]


@dataclass
class HostPresenceResult:
    """Result of a host-presence detection test.

    Attributes:
        n_markers: Number of donor-homozygous markers actually used.
        n_donor_absent_reads: ``Y = sum y_i``, the pooled donor-absent count.
        expected_background: ``Lam = sum n_i * e_i`` under H0 (no host).
        poisson_pval: One-sided pooled-Poisson p-value.
        lrt_pval: Chi-bar-square one-sided LRT p-value.
        f_host_mle: Bounded MLE of the host fraction (>= 0).
        f_host_ci: Profile-likelihood 95% CI for the host fraction.
        used_per_site_error: True when at least one marker used a per-site rate.
        error_rate_source: "per-site" (all markers per-site), "global-fallback"
            (all markers fell back to ``error_rate / 3``), "mixed" (some of
            each), or "none" (no usable markers, the test is degenerate).
        n_artifact_filtered: Donor-homozygous markers dropped by the read-level
            artifact filter before testing (0 when the filter is off or the
            admix VCF carried no bias annotations).
    """

    n_markers: int
    n_donor_absent_reads: int
    expected_background: float
    poisson_pval: float
    lrt_pval: float
    f_host_mle: float
    f_host_ci: tuple[float, float]
    used_per_site_error: bool
    error_rate_source: ErrorRateSource
    n_artifact_filtered: int = 0


@dataclass
class _DonorAbsentMarker:
    """Internal per-marker record threaded through the detector.

    Kept private because the public surface is just ``host_presence_test``;
    exposing the rows directly would lock in a layout that may need to change
    once route B (the unified two-component likelihood) lands.
    """

    key: MarkerKey
    y: int  # donor-absent allele read count
    n: int  # admixture depth
    h: int  # host dose of the donor-absent allele (1 or 2)
    direction: Direction
    # Read-level bias from the admix mpileup, for the artifact filter. Strand
    # counts are for the donor-absent allele specifically; z-scores are
    # site-level. All None when the admix VCF carried no bias annotations.
    da_fwd: int | None = None  # donor-absent allele forward-strand reads
    da_rev: int | None = None  # donor-absent allele reverse-strand reads
    rpbz: float | None = None
    scbz: float | None = None


def _all_donors_uniform_hom(
    donor_gts: list[tuple[int, int]],
) -> tuple[int, int] | None:
    """Return the shared homozygous genotype if every donor is the same hom,
    else ``None``. The detector requires the donor-absent allele to be absent
    from *every* donor (plan section "Marker set").
    """
    if not donor_gts:
        return None
    first = donor_gts[0]
    if first not in ((0, 0), (1, 1)):
        return None
    for d in donor_gts[1:]:
        if d != first:
            return None
    return first


def select_donor_hom_markers(
    informative_markers: list[InformativeMarker],
) -> list[_DonorAbsentMarker]:
    """Pick markers where all donors are homozygous for the same allele.

    Restricts to Vynck types 0, 1, 10, 11 (donor homozygous, host carries the
    donor-absent allele). Types 20 and 21 are excluded because the donor is
    heterozygous and carries both alleles, so there is no donor-absent allele
    to count against a clean background.

    For multi-donor inputs every donor must be homozygous for the same allele;
    a heterozygous or differently-homozygous second donor disqualifies the
    marker.
    """
    out: list[_DonorAbsentMarker] = []
    for m in informative_markers:
        shared = _all_donors_uniform_hom(m.donor_gts)
        if shared is None:
            continue

        # marker_type is computed against the first donor; under the uniform-hom
        # filter above it is also the type for every donor.
        mt = m.marker_type
        if shared == (0, 0):
            # donor-absent allele = ALT; count ALT reads.
            if mt == 1:
                h = 2  # host 1/1
            elif mt == 10:
                h = 1  # host 0/1
            else:
                continue
            direction: Direction = "ref->alt"
            y = m.admix_ad_alt
        else:  # shared == (1, 1)
            # donor-absent allele = REF; count REF reads.
            if mt == 0:
                h = 2  # host 0/0
            elif mt == 11:
                h = 1  # host 0/1
            else:
                continue
            direction = "alt->ref"
            y = m.admix_ad_ref

        # Donor-absent allele strand counts from DP4 (ref_fwd, ref_rev,
        # alt_fwd, alt_rev): the donor-absent allele is ALT for "ref->alt"
        # and REF for "alt->ref".
        da_fwd = da_rev = None
        if m.admix_dp4 is not None:
            rf, rr, af, ar = m.admix_dp4
            if direction == "ref->alt":
                da_fwd, da_rev = af, ar
            else:
                da_fwd, da_rev = rf, rr

        out.append(
            _DonorAbsentMarker(
                key=(m.chrom, m.pos, m.ref, m.alt),
                y=y,
                n=m.admix_dp,
                h=h,
                direction=direction,
                da_fwd=da_fwd,
                da_rev=da_rev,
                rpbz=m.admix_rpbz,
                scbz=m.admix_scbz,
            )
        )
    return out


@dataclass
class ArtifactThresholds:
    """Read-level thresholds for the host-presence artifact filter.

    A donor-homozygous marker is dropped if the donor-absent allele's reads
    look like an alignment artifact rather than a real minor allele. The
    defaults are deliberately loose: genuine artifacts (e.g. the TP53 intron-3
    indel site chr17:7676483) sit far past them (minor-strand fraction <6%,
    |SCBZ| 5-11), while real host signal is strand-balanced (minor strand
    typically >20%) with near-zero bias z-scores, so the filter has a wide
    safety margin.

    Strand bias is judged by effect size (minor-strand fraction), not by a
    significance test: at high depth a real allele's mild 55:45 capture-strand
    skew is highly significant yet harmless, whereas an artifact is extreme
    (~95:5) regardless of depth.

    Attributes:
        min_strand_reads: Minimum donor-absent reads before the strand test
            is applied (below this the strand split is too noisy to judge).
        max_strand_minor_frac: Flag if the lesser strand holds less than this
            fraction of the donor-absent reads (one-sided strand artifact).
        max_abs_scbz: Maximum |soft-clip-length bias z|; above this the
            donor-absent reads are systematically soft-clipped (misaligned).
        max_abs_rpbz: Maximum |read-position bias z|; above this the allele is
            clustered at read positions (ends), a classic artifact signature.
    """

    min_strand_reads: int = 20
    max_strand_minor_frac: float = 0.10
    max_abs_scbz: float = 3.0
    max_abs_rpbz: float = 6.0


def _is_artifact_marker(r: _DonorAbsentMarker, thr: ArtifactThresholds) -> bool:
    """True if a donor-homozygous marker's donor-absent reads look artifactual.

    Three independent signatures, any of which flags the marker:
      1. Soft-clip-length bias (|SCBZ| large): donor-absent reads are
         soft-clipped, i.e. mismapped through a nearby indel/repeat.
      2. Read-position bias (|RPBZ| large): donor-absent reads cluster at read
         ends rather than spanning the site like a real allele.
      3. Strand bias: the donor-absent allele reads come almost entirely from
         one strand, which a genuine amplified/sequenced allele does not.

    Markers with no admix bias annotations (all fields None) never flag, so a
    panel processed without these tags is simply left unfiltered.
    """
    if r.scbz is not None and abs(r.scbz) > thr.max_abs_scbz:
        return True
    if r.rpbz is not None and abs(r.rpbz) > thr.max_abs_rpbz:
        return True
    if r.da_fwd is not None and r.da_rev is not None:
        total = r.da_fwd + r.da_rev
        if total >= thr.min_strand_reads:
            minor_frac = min(r.da_fwd, r.da_rev) / total
            if minor_frac < thr.max_strand_minor_frac:
                return True
    return False


@dataclass(frozen=True)
class DonorHomMarker:
    """Public per-marker view of a donor-homozygous host-presence marker.

    Returned by :func:`donor_hom_markers`. This is the stable surface the
    diagnostic plots consume, so they neither re-run the selection logic nor
    reach into the private ``_DonorAbsentMarker`` / ``_is_artifact_marker``
    internals. ``artifact`` is the verdict of the same filter
    ``host_presence_test`` applies, so a plot can show every marker while
    excluding ``artifact=True`` ones from pooled lines and counts exactly as
    the detector does.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    y: int  # donor-absent allele read count
    n: int  # admix depth
    h: int  # host dose of the donor-absent allele (1 or 2)
    direction: Direction
    artifact: bool

    @property
    def key(self) -> MarkerKey:
        """Marker key ``(chrom, pos, ref, alt)``, for joining to tables."""
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def coef(self) -> float:
        """Per-read host-fraction coefficient ``h / 2`` (``q = e + coef * f_h``)."""
        return self.h / 2.0


def donor_hom_markers(
    informative_markers: list[InformativeMarker],
    artifact_thresholds: ArtifactThresholds | None = None,
) -> list[DonorHomMarker]:
    """Select donor-homozygous markers and flag alignment artifacts.

    Combines :func:`select_donor_hom_markers` with the artifact filter that
    ``host_presence_test`` applies, returning one stable public record per
    marker (artifacts flagged, not dropped). The diagnostic plots consume
    this so the kept/filtered split they draw matches the detector's exactly.

    Args:
        informative_markers: Informative markers from ``classify_markers``.
        artifact_thresholds: Override the artifact-filter thresholds. Defaults
            to ``ArtifactThresholds()``. Markers lacking admix bias
            annotations never flag.

    Returns:
        One ``DonorHomMarker`` per donor-homozygous marker, in selection order.
    """
    thr = artifact_thresholds or ArtifactThresholds()
    out: list[DonorHomMarker] = []
    for r in select_donor_hom_markers(informative_markers):
        chrom, pos, ref, alt = r.key
        out.append(
            DonorHomMarker(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                y=r.y,
                n=r.n,
                h=r.h,
                direction=r.direction,
                artifact=_is_artifact_marker(r, thr),
            )
        )
    return out


def _resolve_e_per_marker(
    rows: list[DonorHomMarker],
    marker_errors: dict[MarkerKey, MarkerErrorRates] | None,
    fallback_e: float,
    error_floor: float,
) -> tuple[list[float], int, int]:
    """Assemble per-marker background rates and tally per-site vs fallback use.

    Returns ``(e_per_marker, n_per_site, n_fallback)``. A per-direction value
    of ``None`` in the table (or a missing key, or no table at all) counts as
    a fallback. The ``error_floor`` is applied uniformly: both the per-site
    loader (in ``allomix.error_rates``) and the fallback go through ``max(.,
    error_floor)`` so a zero rate cannot make a single stray read produce
    -inf log-likelihood.
    """
    fb = max(fallback_e, error_floor)
    e_per_marker: list[float] = []
    n_per_site = 0
    n_fallback = 0
    for r in rows:
        e_i: float | None = None
        if marker_errors is not None:
            entry = marker_errors.get(r.key)
            if entry is not None:
                if r.direction == "ref->alt" and entry.e_refalt is not None:
                    e_i = max(entry.e_refalt, error_floor)
                elif r.direction == "alt->ref" and entry.e_altref is not None:
                    e_i = max(entry.e_altref, error_floor)
        if e_i is None:
            e_i = fb
            n_fallback += 1
        else:
            n_per_site += 1
        e_per_marker.append(e_i)
    return e_per_marker, n_per_site, n_fallback


def _loglik(
    ys: np.ndarray,
    ns: np.ndarray,
    coef: np.ndarray,
    e: np.ndarray,
    f_h: float,
) -> float:
    """Binomial log-likelihood at f_h for the donor-absent allele counts.

    ``q_i(f_h) = e_i + (h_i / 2) * f_h``. Clipped away from 0 and 1 so log()
    stays finite at the boundary; the clip floor sits well below any
    realistic per-site rate.
    """
    q = e + coef * f_h
    q = np.clip(q, 1e-15, 1.0 - 1e-12)
    return float(np.sum(ys * np.log(q) + (ns - ys) * np.log1p(-q)))


def _profile_ci_for_f(
    ys: np.ndarray,
    ns: np.ndarray,
    coef: np.ndarray,
    e: np.ndarray,
    f_hat: float,
    ll_hat: float,
) -> tuple[float, float]:
    """Profile-likelihood 95% CI for the bounded host-fraction MLE.

    Drops 1.92 in log-likelihood from ``ll_hat`` (half the 0.95 quantile of
    chi-sq df=1). Lower bound is clipped at 0 since f_h >= 0 is the parameter
    constraint; upper bound is bracketed in (f_hat, 1].
    """
    drop = chi2.ppf(CI_LEVEL, df=1) / 2.0

    def gap(f_val: float) -> float:
        return ll_hat - _loglik(ys, ns, coef, e, f_val) - drop

    # Lower CI: scan from 0 upward. If ll at 0 is already within the drop,
    # the lower bound is 0 (the boundary).
    if f_hat <= 0.0 or gap(0.0) <= 0.0:
        f_lo = 0.0
    else:
        f_lo = float(brentq(gap, 0.0, f_hat, xtol=1e-9))

    # Upper CI: search upward until LL drops by `drop`. Cap at 1.
    if gap(1.0) <= 0.0:
        f_hi = 1.0
    else:
        f_hi = float(brentq(gap, f_hat, 1.0, xtol=1e-9))

    return f_lo, f_hi


def host_presence_test(
    informative_markers: list[InformativeMarker],
    marker_errors: dict[MarkerKey, MarkerErrorRates] | None = None,
    error_rate: float = DEFAULT_ERROR_RATE,
    error_floor: float = 1e-5,
    artifact_filter: bool = True,
    artifact_thresholds: ArtifactThresholds | None = None,
) -> HostPresenceResult:
    """Run the host-presence detection test.

    Args:
        informative_markers: Informative markers from ``classify_markers``.
            The detector internally selects the donor-homozygous subset.
        marker_errors: Optional per-site, per-direction error table from
            ``allomix.error_rates.load_error_table``. When provided, each
            marker uses its per-direction rate; when missing or ``None`` for
            the relevant direction the detector falls back to
            ``error_rate / 3`` (the per-direction floor implied by the
            symmetric 4-state error model).
        error_rate: Global symmetric sequencing error rate. Used as the
            per-direction fallback ``error_rate / 3`` for markers missing
            from the table.
        error_floor: Per-direction lower bound applied to every per-marker
            background rate. Prevents a zero rate from producing -inf
            log-likelihood on a single stray read.
        artifact_filter: When True (default), drop donor-homozygous markers
            whose donor-absent reads show alignment-artifact signatures
            (strand bias, soft-clip bias, read-position bias) before testing.
            See ``ArtifactThresholds``. Requires admix bias annotations
            (DP4/SCBZ/RPBZ); a no-op on VCFs lacking them.
        artifact_thresholds: Thresholds for the artifact filter. Defaults to
            ``ArtifactThresholds()`` when None.

    Returns:
        A ``HostPresenceResult`` summarising both statistics, the MLE host
        fraction and its profile CI, and the provenance of the background
        rates.
    """
    markers = donor_hom_markers(informative_markers, artifact_thresholds)
    fallback_e = error_rate / N_OTHER_BASES

    if artifact_filter:
        rows = [m for m in markers if not m.artifact]
        n_artifact_filtered = len(markers) - len(rows)
    else:
        rows = markers
        n_artifact_filtered = 0

    if not rows:
        return HostPresenceResult(
            n_markers=0,
            n_donor_absent_reads=0,
            expected_background=0.0,
            poisson_pval=1.0,
            lrt_pval=1.0,
            f_host_mle=0.0,
            f_host_ci=(0.0, 0.0),
            used_per_site_error=False,
            error_rate_source="none",
            n_artifact_filtered=n_artifact_filtered,
        )

    e_list, n_per_site, n_fallback = _resolve_e_per_marker(
        rows, marker_errors, fallback_e, error_floor,
    )

    if n_per_site and n_fallback:
        source: ErrorRateSource = "mixed"
    elif n_per_site:
        source = "per-site"
    else:
        source = "global-fallback"

    ys = np.asarray([r.y for r in rows], dtype=float)
    ns = np.asarray([r.n for r in rows], dtype=float)
    hs = np.asarray([r.h for r in rows], dtype=float)
    coef = hs / 2.0
    e = np.asarray(e_list, dtype=float)

    # Pooled Poisson under H0: y_i ~ Binomial(n_i, e_i), summed -> Poisson(Lam)
    # is accurate because e_i is tiny and n_i is large.
    Y = int(ys.sum())
    Lam = float((ns * e).sum())
    if Lam <= 0.0:
        # Should not happen with a positive error_floor, but keep the path
        # honest.
        p_pois = 1.0
    else:
        p_pois = 1.0 if Y == 0 else float(poisson.sf(Y - 1, Lam))

    # Bounded MLE for f_h. The likelihood is concave in f_h for fixed
    # (positive) backgrounds, so a bounded scalar minimiser is enough.
    ll0 = _loglik(ys, ns, coef, e, 0.0)
    res = minimize_scalar(
        lambda fh: -_loglik(ys, ns, coef, e, fh),
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-9},
    )
    f_hat = float(max(0.0, min(1.0, res.x)))
    ll_hat = -float(res.fun)

    # Chi-bar-square boundary correction: f_h >= 0, so under H0 the LRT is a
    # 50:50 mixture of chi2_0 (point mass at 0) and chi2_1. If the MLE is at
    # the boundary, the statistic is 0 and p = 1.
    if ll_hat <= ll0 + 1e-9 or f_hat <= 0.0:
        return HostPresenceResult(
            n_markers=len(rows),
            n_donor_absent_reads=Y,
            expected_background=Lam,
            poisson_pval=p_pois,
            lrt_pval=1.0,
            f_host_mle=0.0,
            f_host_ci=(0.0, 0.0),
            used_per_site_error=n_per_site > 0,
            error_rate_source=source,
            n_artifact_filtered=n_artifact_filtered,
        )

    D = 2.0 * (ll_hat - ll0)
    p_lrt = 0.5 * float(chi2.sf(D, 1)) if D > 0 else 1.0
    f_lo, f_hi = _profile_ci_for_f(ys, ns, coef, e, f_hat, ll_hat)

    return HostPresenceResult(
        n_markers=len(rows),
        n_donor_absent_reads=Y,
        expected_background=Lam,
        poisson_pval=p_pois,
        lrt_pval=p_lrt,
        f_host_mle=f_hat,
        f_host_ci=(f_lo, f_hi),
        used_per_site_error=n_per_site > 0,
        error_rate_source=source,
        n_artifact_filtered=n_artifact_filtered,
    )


__all__ = [
    "ArtifactThresholds",
    "Direction",
    "DonorHomMarker",
    "ErrorRateSource",
    "HostPresenceResult",
    "MarkerKey",
    "donor_hom_markers",
    "host_presence_test",
    "select_donor_hom_markers",
]
