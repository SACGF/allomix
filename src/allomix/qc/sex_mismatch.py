"""Independent sex-chromosome cross-check of the donor fraction (#48).

For a sex-mismatched host/donor pair the sex chromosomes carry a second,
independent readout of the mixing fraction, which this module turns into an
estimate with its own n and confidence interval and compares with the
autosomal MLE. The headline ``donor_pct`` is never blended with it; a
discordant pair of estimates is a QC signal (``qc.assess_quality`` promotes
a QC warning without changing the status), not a correction.

Active only when ``SexResult.pair`` is ``PairStatus.MISMATCHED`` and there is
exactly one donor. With two donors the two sex-chromosome readouts below have
two unknown fractions to split between the parties and are not identifiable
from a single chrY depth or a single chrX copy-number contrast, so the check
returns None and the report shows NA.

Two readouts, in priority order:

1. **chrY depth ratio** (experimental; no data in hand has non-PAR chrY
   content, so this path is validated on synthetic depth tables only). The lab
   statistic ``R = max(DP at non-PAR chrY sites) / median(DP at autosomal
   sites)`` in the admixture sample, normalised by the same ratio in the pair's
   male reference sample (same panel, same chemistry), gives the male fraction;
   the CI is a bootstrap over the chrY sites. Depth comes from a forced pileup
   at the panel's interval midpoints (``region_depths_from_vcf``), since GATK
   emits no record at an invariant chrY site. Skipped when no depth input is
   given or fewer than ``SEXCHROM_MIN_CHRY_SITES`` chrY sites are present.
2. **chrX copy-number-weighted allele fraction.** At non-PAR chrX markers where
   host and donor genotypes differ, a female contributes two X copies and a
   male one, so the expected REF fraction at donor fraction ``f`` is

       w(f) = [(1 - f) cn_h ref_h + f cn_d ref_d] / [(1 - f) cn_h + f cn_d]

   with ``ref_x`` the party's REF allele fraction (1, 0.5 or 0 for a diploid
   party; 1 or 0 for a hemizygous one). ``f`` is fitted by maximising the
   beta-binomial log-likelihood over the admixture allele counts with the same
   4-state (or per-marker asymmetric) error model, logit-space bias correction
   and profiled overdispersion the main estimator uses (grid, then Brent, then
   a profile-likelihood CI). A het call in the male party is a genotyping error
   on a hemizygous chromosome and that marker is dropped and counted. The
   signal is asymmetric: a male minor at a female-homozygous site gives an ALT
   fraction of ``f / (2 - f)`` (half the autosomal contrast), a female hom-alt
   minor at a male hom-ref site gives ``2f / (1 + f)`` (double).

``basis`` records which readout the headline ``frac_donor`` came from:
``"chrY-depth"`` when it ran with enough sites, else ``"chrX-cn"`` when at
least ``SEXCHROM_MIN_MARKERS`` usable chrX markers were found, else None. Both
readouts are kept on the result when both ran. Concordance with the autosomal
MLE is CI overlap, or an absolute difference within
``SEXCHROM_CONCORDANCE_PP`` percentage points.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from cyvcf2 import VCF
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit, gammaln, logit
from scipy.stats import chi2

from allomix.constants import (
    CI_LEVEL,
    DEFAULT_ERROR_RATE,
    N_OTHER_BASES,
    PLOIDY,
    SEXCHROM_CONCORDANCE_PP,
    SEXCHROM_MIN_CHRY_SITES,
    SEXCHROM_MIN_MARKERS,
)
from allomix.contigs import ContigClass, classify_contig
from allomix.estimate.likelihood import W_EPS, PanelCalibration
from allomix.genotype import ContigPolicy, InformativeMarker, MarkerData, classify_markers
from allomix.qc.sex import RegionDepth, SexResult, chry_depth_ratio, split_region_depths
from allomix.sex_types import PairStatus, Sex

#: chrX copy number per sex (non-PAR).
CHRX_COPY_NUMBER: dict[Sex, int] = {Sex.FEMALE: 2, Sex.MALE: 1}

# Grid resolution for the chrX fit before Brent refinement, and the rho range
# profiled at each f (same bounds as the main estimator's profile).
_GRID_STEPS = 201
_RHO_MIN = 1.0
_RHO_MAX = 50000.0
_P_EPS = 1e-6
# Bootstrap replicates for the chrY depth-ratio CI (deterministic seed).
_CHRY_BOOTSTRAP = 2000
_CHRY_SEED = 20260048

CHRY_BASIS = "chrY-depth"
CHRX_BASIS = "chrX-cn"


@dataclass
class SexMismatchResult:
    """Sex-chromosome cross-check of the donor fraction for a mismatched pair.

    Attributes:
        male_party: ``"host"`` or ``"donor"``.
        basis: Readout behind ``frac_donor``: ``"chrY-depth"``, ``"chrX-cn"``,
            or None when neither readout had enough data.
        frac_donor: Donor fraction (0-1) implied by the sex chromosomes, directly
            comparable with the MLE ``donor_fraction``; None when ``basis`` is None.
        ci_low: Lower 95% bound of ``frac_donor``, or None.
        ci_high: Upper 95% bound, or None.
        n: Markers (chrX) or chrY sites behind the headline estimate.
        concordant: True when the headline estimate agrees with the autosomal
            MLE (CI overlap, or within ``SEXCHROM_CONCORDANCE_PP`` points); None
            when there is no estimate.
        mle_frac_donor: The autosomal MLE it was compared with.
        mle_ci: The MLE's 95% CI.
        chrx_n: Usable non-PAR chrX markers (host and donor differ, male party
            homozygous) that went into the chrX fit.
        chrx_n_male_het_dropped: Informative chrX markers dropped because the male
            party has a het call there.
        chrx_frac_donor: chrX copy-number-weighted estimate, or None when fewer
            than ``SEXCHROM_MIN_MARKERS`` markers were usable.
        chrx_ci: Its profile-likelihood 95% CI, or None.
        chrx_rho: Beta-binomial concentration fitted at the chrX estimate, or None.
        chry_n_sites: Non-PAR chrY depth sites seen in the admixture depth input
            (0 when no depth input was given).
        chry_frac_donor: chrY depth-ratio estimate, or None when not run.
        chry_ci: Its bootstrap 95% CI, or None. The whole [0, 1] when fewer than
            ``SEXCHROM_MIN_CHRY_SITES`` sites were available.
        chry_ratio_admix: ``max(chrY DP) / median(autosomal DP)`` in the admixture.
        chry_ratio_ref: The same ratio in the male reference sample.
        chry_ci_unreliable: True when the chrY CI is the [0, 1] placeholder.
    """

    male_party: str
    basis: str | None
    frac_donor: float | None
    ci_low: float | None
    ci_high: float | None
    n: int
    concordant: bool | None
    mle_frac_donor: float
    mle_ci: tuple[float, float]
    chrx_n: int = 0
    chrx_n_male_het_dropped: int = 0
    chrx_frac_donor: float | None = None
    chrx_ci: tuple[float, float] | None = None
    chrx_rho: float | None = None
    chry_n_sites: int = 0
    chry_frac_donor: float | None = None
    chry_ci: tuple[float, float] | None = None
    chry_ratio_admix: float | None = None
    chry_ratio_ref: float | None = None
    chry_ci_unreliable: bool = False


# ---------------------------------------------------------------------------
# chrX copy-number-weighted allele-fraction readout
# ---------------------------------------------------------------------------


def ref_fraction(gt: tuple[int, int]) -> float:
    """REF allele fraction of a party's own genome at a bi-allelic marker.

    ``1 - alt_dose / PLOIDY``: 1, 0.5 or 0 for a diploid party. A hemizygous
    party is written diploid-homozygous by GATK, so its 1 or 0 comes out the
    same way; a het call on a hemizygous chromosome must be dropped first.
    """
    return 1.0 - (gt[0] + gt[1]) / PLOIDY


def expected_ref_fraction(
    f_donor: float | np.ndarray,
    ref_host: float | np.ndarray,
    ref_donor: float | np.ndarray,
    cn_host: int,
    cn_donor: int,
) -> float | np.ndarray:
    """Copy-number-weighted expected REF fraction in the mixture.

    ``w(f) = [(1 - f) cn_h ref_h + f cn_d ref_d] / [(1 - f) cn_h + f cn_d]``.
    With ``cn_host == cn_donor`` the copy numbers cancel and this is the diploid
    ``expected_weight`` of the main estimator.
    """
    f_host = 1.0 - f_donor
    num = f_host * cn_host * ref_host + f_donor * cn_donor * ref_donor
    den = f_host * cn_host + f_donor * cn_donor
    return num / den


@dataclass(frozen=True)
class _ChrXArrays:
    """Per-marker arrays for the chrX fit (everything that does not depend on f, rho)."""

    ref_host: np.ndarray
    ref_donor: np.ndarray
    n: np.ndarray
    k: np.ndarray  # ALT reads
    bias: np.ndarray
    bias_mask: np.ndarray
    e_refalt: np.ndarray
    e_altref: np.ndarray
    error_mask: np.ndarray
    cn_host: int
    cn_donor: int


def _chrx_arrays(
    markers: list[InformativeMarker],
    cn_host: int,
    cn_donor: int,
    calibration: PanelCalibration,
) -> _ChrXArrays:
    n_m = len(markers)
    bias = np.array([calibration.bias_for(m) for m in markers], dtype=float)
    e_refalt = np.full(n_m, np.nan)
    e_altref = np.full(n_m, np.nan)
    for i, m in enumerate(markers):
        entry = calibration.error_for(m)
        if entry is not None and entry.e_refalt is not None and entry.e_altref is not None:
            e_refalt[i] = entry.e_refalt
            e_altref[i] = entry.e_altref
    return _ChrXArrays(
        ref_host=np.array([ref_fraction(m.host_gt) for m in markers], dtype=float),
        ref_donor=np.array([ref_fraction(m.donor_gts[0]) for m in markers], dtype=float),
        n=np.array([m.admix_ad_ref + m.admix_ad_alt for m in markers], dtype=float),
        k=np.array([m.admix_ad_alt for m in markers], dtype=float),
        bias=bias,
        bias_mask=bias != 0.0,
        e_refalt=e_refalt,
        e_altref=e_altref,
        error_mask=~np.isnan(e_refalt),
        cn_host=cn_host,
        cn_donor=cn_donor,
    )


def _p_alt(arr: _ChrXArrays, f_donor: float, error_rate: float) -> np.ndarray:
    """Per-marker P(observe ALT) at donor fraction ``f_donor``.

    Mirrors ``estimate.likelihood._p_alt_for_f`` with the copy-number-weighted
    expectation in place of the diploid one: logit-space bias correction at
    the biased markers, then the 4-state symmetric error model, replaced by the
    asymmetric REF/ALT-only model where both per-direction rates are known.
    """
    w = expected_ref_fraction(f_donor, arr.ref_host, arr.ref_donor, arr.cn_host, arr.cn_donor)
    w = np.asarray(w, dtype=float)
    if arr.bias_mask.any():
        wm = np.clip(w[arr.bias_mask], W_EPS, 1.0 - W_EPS)
        pb = np.clip(0.5 + arr.bias[arr.bias_mask], W_EPS, 1.0 - W_EPS)
        w[arr.bias_mask] = np.clip(expit(logit(wm) - logit(pb)), W_EPS, 1.0 - W_EPS)
    e = error_rate
    e_specific = e / N_OTHER_BASES
    p_alt_raw = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref_raw = w * (1.0 - e) + (1.0 - w) * e_specific
    p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)
    if arr.error_mask.any():
        p_alt_asym = w * arr.e_refalt + (1.0 - w) * (1.0 - arr.e_altref)
        p_alt = np.where(arr.error_mask, p_alt_asym, p_alt)
    return np.clip(p_alt, _P_EPS, 1.0 - _P_EPS)


def _bb_loglik(n: np.ndarray, k: np.ndarray, p_alt: np.ndarray, rho: float) -> float:
    """Beta-binomial log-likelihood (constant ``log C(n, k)`` dropped)."""
    a = np.maximum(p_alt * rho, 1e-10)
    b = np.maximum((1.0 - p_alt) * rho, 1e-10)
    ll = gammaln(k + a) + gammaln(n - k + b) - gammaln(n + rho) - gammaln(a) - gammaln(b)
    ll += gammaln(rho)
    return float(ll.sum())


def _profile_ll(arr: _ChrXArrays, f_donor: float, error_rate: float) -> tuple[float, float]:
    """Max log-likelihood over rho at fixed f. Returns ``(ll, rho)``."""
    p_alt = _p_alt(arr, f_donor, error_rate)
    opt = minimize_scalar(
        lambda log_r: -_bb_loglik(arr.n, arr.k, p_alt, math.exp(log_r)),
        bounds=(math.log(_RHO_MIN), math.log(_RHO_MAX)),
        method="bounded",
    )
    return -float(opt.fun), math.exp(float(opt.x))


def fit_chrx_fraction(
    markers: list[InformativeMarker],
    cn_host: int,
    cn_donor: int,
    *,
    error_rate: float = DEFAULT_ERROR_RATE,
    calibration: PanelCalibration | None = None,
) -> tuple[float, tuple[float, float], float]:
    """Fit the donor fraction from chrX allele counts under the copy-number model.

    Grid search over ``f`` in [0, 1] with rho profiled out at each point, Brent
    refinement around the best grid point, then a profile-likelihood CI at
    ``CI_LEVEL``.

    Args:
        markers: Non-PAR chrX markers where host and donor differ, the male party
            homozygous at each (see ``usable_chrx_markers``).
        cn_host: Host chrX copy number (1 or 2).
        cn_donor: Donor chrX copy number (1 or 2).
        error_rate: Symmetric 4-state sequencing error rate (fallback where a
            marker has no per-direction rates in ``calibration``).
        calibration: Optional per-marker bias and error tables, applied as in
            the main estimator.

    Returns:
        ``(f_mle, (ci_low, ci_high), rho_mle)``.

    Raises:
        ValueError: If ``markers`` is empty.
    """
    if not markers:
        raise ValueError("fit_chrx_fraction needs at least one marker")
    arr = _chrx_arrays(markers, cn_host, cn_donor, calibration or PanelCalibration())

    grid = np.linspace(0.0, 1.0, _GRID_STEPS)
    lls = [_profile_ll(arr, float(f), error_rate)[0] for f in grid]
    i_best = int(np.argmax(lls))
    step = 1.0 / (_GRID_STEPS - 1)
    lo_b = max(0.0, grid[i_best] - step)
    hi_b = min(1.0, grid[i_best] + step)
    opt = minimize_scalar(
        lambda f: -_profile_ll(arr, f, error_rate)[0],
        bounds=(lo_b, hi_b),
        method="bounded",
        options={"xatol": 1e-7},
    )
    f_mle = float(min(max(opt.x, 0.0), 1.0))
    ll_max, rho_mle = _profile_ll(arr, f_mle, error_rate)

    half_threshold = chi2.ppf(CI_LEVEL, df=1) / 2.0

    def gap(f: float) -> float:
        return ll_max - _profile_ll(arr, f, error_rate)[0] - half_threshold

    if f_mle <= 0.0 or gap(0.0) <= 0.0:
        f_lo = 0.0
    else:
        f_lo = float(brentq(gap, 0.0, f_mle, xtol=1e-6))
    if f_mle >= 1.0 or gap(1.0) <= 0.0:
        f_hi = 1.0
    else:
        f_hi = float(brentq(gap, f_mle, 1.0, xtol=1e-6))
    return f_mle, (f_lo, f_hi), rho_mle


def _is_het(gt: tuple[int, int]) -> bool:
    return gt[0] != gt[1]


def usable_chrx_markers(
    host: list[MarkerData],
    donor: list[MarkerData],
    admix: list[MarkerData],
    male_party: str,
    *,
    min_dp: int,
    min_gq: int,
) -> tuple[list[InformativeMarker], int]:
    """Select the non-PAR chrX markers for the copy-number fit.

    Runs ``classify_markers`` under the diagnostic ``ALL_PRIMARY`` policy (so
    the same PASS, GQ and depth filters apply as in the main analysis) and
    keeps the informative non-PAR chrX markers at which the male party is
    homozygous.

    Returns:
        ``(markers, n_male_het_dropped)``.
    """
    g = classify_markers(
        host,
        [donor],
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        contig_policy=ContigPolicy.ALL_PRIMARY,
    )
    kept: list[InformativeMarker] = []
    n_dropped = 0
    for m in g.informative:
        if classify_contig(m.chrom, m.pos) is not ContigClass.X_NONPAR:
            continue
        male_gt = m.host_gt if male_party == "host" else m.donor_gts[0]
        if _is_het(male_gt):
            n_dropped += 1
            continue
        kept.append(m)
    return kept, n_dropped


# ---------------------------------------------------------------------------
# chrY depth-ratio readout (experimental)
# ---------------------------------------------------------------------------


def region_depths_from_vcf(path: str, sample: str) -> list[RegionDepth]:
    """Read per-site depth for one sample from a forced-pileup VCF.

    The pipeline's midpoint pileups (``refs/midpoints.vcf.gz`` for the reference
    samples, ``<patient>.admix.midpoints.vcf.gz`` for the admixture samples)
    carry one record per panel interval with ``FORMAT/DP``. Sites with a missing
    DP are skipped.

    Raises:
        ValueError: If ``sample`` is not in the VCF.
    """
    vcf = VCF(str(path))
    if sample not in vcf.samples:
        raise ValueError(f"sample {sample!r} not found in {path} (samples: {vcf.samples})")
    idx = vcf.samples.index(sample)
    out: list[RegionDepth] = []
    for rec in vcf:
        dp_arr = rec.format("DP")
        if dp_arr is None:
            continue
        dp = int(dp_arr[idx][0])
        if dp < 0:  # cyvcf2 missing sentinel
            continue
        out.append((rec.CHROM, rec.POS, dp))
    return out


@dataclass(frozen=True)
class ChrYEstimate:
    """chrY depth-ratio estimate of the male fraction in the admixture."""

    frac_male: float
    ci: tuple[float, float]
    n_sites: int
    ratio_admix: float
    ratio_ref: float
    ci_unreliable: bool


def chry_male_fraction(
    admix_depths: list[RegionDepth],
    male_ref_depths: list[RegionDepth],
) -> ChrYEstimate | None:
    """Male fraction of the admixture from the chrY depth ratio (experimental).

    ``frac_male = R_admix / R_male_ref`` clipped to [0, 1], where ``R`` is
    ``chry_depth_ratio``. The CI is a percentile bootstrap over the chrY sites
    (the same resampled sites in both samples, since they are the same
    regions). With fewer than ``SEXCHROM_MIN_CHRY_SITES`` chrY sites the CI is
    the whole [0, 1] and ``ci_unreliable`` is set.

    Returns None when either sample has no chrY or no autosomal depth, or the
    reference ratio is zero (a male reference with no chrY coverage cannot
    normalise anything).
    """
    r_admix, n_admix = chry_depth_ratio(admix_depths)
    r_ref, _ = chry_depth_ratio(male_ref_depths)
    if r_admix is None or r_ref is None or r_ref <= 0.0:
        return None
    frac = min(max(r_admix / r_ref, 0.0), 1.0)

    y_a, auto_a = split_region_depths(admix_depths)
    y_r, auto_r = split_region_depths(male_ref_depths)
    n_sites = int(min(y_a.size, y_r.size))
    if n_sites < SEXCHROM_MIN_CHRY_SITES:
        return ChrYEstimate(frac, (0.0, 1.0), n_sites, r_admix, r_ref, True)

    # Site-matched bootstrap: resample chrY site indices with replacement and
    # recompute both ratios on the same draw. The site lists are matched by
    # position so a draw refers to the same regions in both samples.
    pos_a = {
        (c, p): dp for c, p, dp in admix_depths if classify_contig(c, p) is ContigClass.Y_NONPAR
    }
    pos_r = {
        (c, p): dp for c, p, dp in male_ref_depths if classify_contig(c, p) is ContigClass.Y_NONPAR
    }
    shared = sorted(set(pos_a) & set(pos_r))
    if len(shared) < SEXCHROM_MIN_CHRY_SITES:
        return ChrYEstimate(frac, (0.0, 1.0), len(shared), r_admix, r_ref, True)
    ya = np.array([pos_a[s] for s in shared], dtype=float)
    yr = np.array([pos_r[s] for s in shared], dtype=float)
    med_a = float(np.median(auto_a))
    med_r = float(np.median(auto_r))
    rng = np.random.default_rng(_CHRY_SEED)
    draws = rng.integers(0, len(shared), size=(_CHRY_BOOTSTRAP, len(shared)))
    r_a = ya[draws].max(axis=1) / med_a
    r_r = yr[draws].max(axis=1) / med_r
    with np.errstate(divide="ignore", invalid="ignore"):
        boots = np.clip(np.where(r_r > 0, r_a / r_r, np.nan), 0.0, 1.0)
    boots = boots[~np.isnan(boots)]
    if boots.size == 0:
        return ChrYEstimate(frac, (0.0, 1.0), len(shared), r_admix, r_ref, True)
    alpha = (1.0 - CI_LEVEL) / 2.0
    lo = float(np.quantile(boots, alpha))
    hi = float(np.quantile(boots, 1.0 - alpha))
    return ChrYEstimate(frac, (min(lo, frac), max(hi, frac)), len(shared), r_admix, r_ref, False)


# ---------------------------------------------------------------------------
# Concordance and the top-level check
# ---------------------------------------------------------------------------


def is_concordant(
    frac: float,
    ci: tuple[float, float],
    mle_frac: float,
    mle_ci: tuple[float, float],
    tolerance_pp: float = SEXCHROM_CONCORDANCE_PP,
) -> bool:
    """True when the two estimates' CIs overlap or they differ by at most ``tolerance_pp``."""
    overlap = ci[0] <= mle_ci[1] and mle_ci[0] <= ci[1]
    return overlap or abs(frac - mle_frac) * 100.0 <= tolerance_pp


def _to_donor(frac_male: float, ci: tuple[float, float], male_party: str):
    """Convert a male-fraction estimate to a donor-fraction one."""
    if male_party == "donor":
        return frac_male, ci
    return 1.0 - frac_male, (1.0 - ci[1], 1.0 - ci[0])


def sex_mismatch_check(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    sex: SexResult | None,
    *,
    mle_frac_donor: float,
    mle_ci: tuple[float, float],
    min_dp: int,
    min_gq: int,
    error_rate: float = DEFAULT_ERROR_RATE,
    calibration: PanelCalibration | None = None,
    admix_region_depth: list[RegionDepth] | None = None,
    male_ref_region_depth: list[RegionDepth] | None = None,
) -> SexMismatchResult | None:
    """Run the sex-chromosome cross-check for one admixture sample.

    Args:
        host: Parsed host markers.
        donors: Parsed donor markers (one list per donor).
        admix: Parsed admixture markers (``min_dp`` is applied here).
        sex: Sex inference for the pair; the check runs only for
            ``PairStatus.MISMATCHED``.
        mle_frac_donor: The autosomal MLE donor fraction to compare with.
        mle_ci: Its 95% CI.
        min_dp: Minimum admixture depth at a chrX marker.
        min_gq: Minimum reference-sample GQ at a chrX marker.
        error_rate: Symmetric sequencing error rate (fallback).
        calibration: Optional per-marker bias and error tables.
        admix_region_depth: Forced-pileup depths for the admixture sample, for the
            chrY readout; None skips it.
        male_ref_region_depth: The same for the pair's male reference sample.

    Returns:
        A ``SexMismatchResult``, or None when the pair is not sex-mismatched,
        when sex was not assessed, or when there is more than one donor.
    """
    if sex is None or sex.pair is not PairStatus.MISMATCHED or len(donors) != 1:
        return None
    host_sex = sex.host.effective
    donor_sex = sex.donors[0].effective
    male_party = "host" if host_sex is Sex.MALE else "donor"
    result = SexMismatchResult(
        male_party=male_party,
        basis=None,
        frac_donor=None,
        ci_low=None,
        ci_high=None,
        n=0,
        concordant=None,
        mle_frac_donor=mle_frac_donor,
        mle_ci=mle_ci,
    )

    # Readout 2: chrX copy-number-weighted allele fraction.
    markers, n_het_dropped = usable_chrx_markers(
        host, donors[0], admix, male_party, min_dp=min_dp, min_gq=min_gq
    )
    result.chrx_n = len(markers)
    result.chrx_n_male_het_dropped = n_het_dropped
    if len(markers) >= SEXCHROM_MIN_MARKERS:
        f, ci, rho = fit_chrx_fraction(
            markers,
            CHRX_COPY_NUMBER[host_sex],
            CHRX_COPY_NUMBER[donor_sex],
            error_rate=error_rate,
            calibration=calibration,
        )
        result.chrx_frac_donor = f
        result.chrx_ci = ci
        result.chrx_rho = rho

    # Readout 1 (experimental): chrY depth ratio, when depth input is present.
    if admix_region_depth is not None and male_ref_region_depth is not None:
        result.chry_n_sites = chry_depth_ratio(admix_region_depth)[1]
        est = chry_male_fraction(admix_region_depth, male_ref_region_depth)
        if est is not None:
            fd, cid = _to_donor(est.frac_male, est.ci, male_party)
            result.chry_frac_donor = fd
            result.chry_ci = cid
            result.chry_ratio_admix = est.ratio_admix
            result.chry_ratio_ref = est.ratio_ref
            result.chry_ci_unreliable = est.ci_unreliable
            result.chry_n_sites = est.n_sites

    # Headline: chrY when it ran with enough sites, else chrX.
    if result.chry_frac_donor is not None and not result.chry_ci_unreliable:
        result.basis = CHRY_BASIS
        result.frac_donor = result.chry_frac_donor
        result.ci_low, result.ci_high = result.chry_ci
        result.n = result.chry_n_sites
    elif result.chrx_frac_donor is not None:
        result.basis = CHRX_BASIS
        result.frac_donor = result.chrx_frac_donor
        result.ci_low, result.ci_high = result.chrx_ci
        result.n = result.chrx_n
    if result.frac_donor is not None:
        result.concordant = is_concordant(
            result.frac_donor, (result.ci_low, result.ci_high), mle_frac_donor, mle_ci
        )
    return result


__all__ = [
    "CHRX_BASIS",
    "CHRX_COPY_NUMBER",
    "CHRY_BASIS",
    "ChrYEstimate",
    "RegionDepth",
    "SexMismatchResult",
    "chry_depth_ratio",
    "chry_male_fraction",
    "expected_ref_fraction",
    "fit_chrx_fraction",
    "is_concordant",
    "ref_fraction",
    "region_depths_from_vcf",
    "sex_mismatch_check",
    "usable_chrx_markers",
]
