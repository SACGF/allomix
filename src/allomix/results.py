"""Output data types for chimerism estimation.

The result objects produced by the estimators in ``allomix.chimerism`` and
consumed by ``report.py``, ``qc.py``, and ``analysis.py``. None of those
consumers call an estimator; they only read these result objects, so the data
types live in their own module. They aggregate result types from the identity
and contamination QC modules (``detect``, ``contamination``, ``relatedness``,
``runmeta``), none of which import this module, so there is no import cycle.
"""

from dataclasses import dataclass

from allomix.contamination import ContaminationResult
from allomix.detect import HostPresenceResult
from allomix.relatedness import AdmixConsistencyResult, RelatednessResult
from allomix.runmeta import RunUnitInfo


@dataclass
class MarkerResult:
    """Per-marker contribution to the chimerism estimate."""

    chrom: str
    pos: int
    marker_type: int
    expected_vaf: float
    observed_vaf: float
    residual: float
    ad_ref: int
    ad_alt: int
    dp: int
    included: bool  # False if outlier-excluded


@dataclass
class ChimerismResult:
    """Result of single-donor chimerism estimation."""

    donor_fraction: float  # MLE point estimate (0.0-1.0)
    donor_fraction_ci: tuple[float, float]  # 95% CI
    host_fraction: float  # 1 - donor_fraction
    log_likelihood: float  # at MLE
    n_informative: int
    n_markers_used: int  # after outlier exclusion
    per_marker: list[MarkerResult]
    error_rate: float
    rho: float = float("inf")  # beta-binomial concentration; inf = no overdispersion
    # Per-marker-type overdispersion (issue #33). None unless that mode ran with
    # both classes above MIN_CLASS_MARKERS. rho_het governs the donor-het class
    # (background VAF ~0.5, the over-dispersed one); rho_hom the donor-hom class.
    # In two-rho mode the headline ``rho`` above is set to rho_het.
    rho_hom: float | None = None
    rho_het: float | None = None
    # Set to a human-readable reason when --marker-type-overdispersion was
    # requested but a class was too sparse to identify its rho, so the estimator
    # fell back to shared rho. None otherwise. Surfaced as a QC warning because
    # stderr is lost when the CLI runs as a subprocess.
    marker_type_overdispersion_fallback: str | None = None
    # Per-sample analytical detection limits (donor fractions, 0.0-1.0), computed
    # from the Fisher information of this sample's own markers. inf = nothing detectable.
    lob_fraction: float = float("inf")  # limit of blank
    lod_fraction: float = float("inf")  # limit of detection
    # Host-presence detector output (see ``allomix.detect``). None when the
    # caller disabled the detector or when there were no donor-homozygous
    # markers to run it on.
    host_presence: HostPresenceResult | None = None
    # Robust-refit accounting (see ``estimate_single_donor_bb`` robust mode).
    # n_robust_excluded is the count of markers dropped as residual outliers and
    # excluded from the final fit; robust_drop_fraction is that count over
    # n_informative. Both are 0 when robust mode is off or nothing was dropped.
    n_robust_excluded: int = 0
    robust_drop_fraction: float = 0.0
    # Identity QC (see ``allomix.relatedness``), attached by ``analyse_sample``.
    # relatedness holds one entry per reference-sample pair (host vs each donor);
    # admix_consistency is the consensus-homozygote swap check. None when the
    # caller did not compute them.
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None
    # In-data contamination estimate (see ``allomix.contamination``), attached by
    # ``analyse_sample``. None when the caller did not compute it.
    contamination: ContaminationResult | None = None
    # Sequencing run-unit metadata for this admix sample, read from the admix VCF
    # header (see ``allomix.runmeta``). None when the VCF carried no run metadata.
    run_unit: RunUnitInfo | None = None


@dataclass
class MultiDonorResult:
    """Result of multi-donor chimerism estimation."""

    donor_fractions: list[float]  # [f_donor1, f_donor2, ...]
    donor_fraction_cis: list[tuple[float, float]]  # [(lo, hi), ...] per donor
    host_fraction: float  # 1 - sum(donor_fractions)
    log_likelihood: float
    n_informative: int
    n_markers_used: int
    per_marker: list[MarkerResult]
    error_rate: float
    per_donor_n_informative: list[int] | None = None  # informative markers per donor
    rho: float = float("inf")  # beta-binomial concentration; inf = no overdispersion
    host_presence: HostPresenceResult | None = None
    n_robust_excluded: int = 0
    robust_drop_fraction: float = 0.0
    # Identity QC; see ChimerismResult above. For multi-donor, relatedness also
    # includes donor-vs-donor pairs.
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None
    contamination: ContaminationResult | None = None
    run_unit: RunUnitInfo | None = None
