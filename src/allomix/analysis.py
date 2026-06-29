"""Single-sample analysis pipeline shared by the CLI and diagnostic scripts.

Given pre-parsed host, donor and admixture markers, ``analyse_sample`` runs the
full per-sample path once: classify markers, estimate the chimerism fraction
(single- or multi-donor), run the host-presence detector, assess QC, and select
the donor-homozygous markers with their artifact flags. Both ``allomix.cli`` and
the ``scripts/`` diagnostics call this so the genotype -> classify ->
select-donor-hom -> presence path lives in exactly one place and every knob
(sex chromosomes, artifact filter, bias/error tables) is handled identically.

This is library code: it does not read VCFs (callers parse first, since the
input conventions differ) and it does not print. Callers own I/O and messaging.
"""

from dataclasses import dataclass

from allomix.chimerism import estimate_multi_donor, estimate_single_donor_bb
from allomix.constants import ROBUST_K_DEFAULT
from allomix.contamination import estimate_contamination
from allomix.detect import DonorHomMarker, donor_hom_markers, host_presence_test
from allomix.genotype import MarkerData, MarkerGenotypes, classify_markers
from allomix.likelihood import PanelCalibration
from allomix.qc import QCReport, assess_quality
from allomix.relatedness import (
    RelatednessResult,
    admix_consistency,
    relatedness_coefficient,
)
from allomix.results import ChimerismResult, MultiDonorResult
from allomix.runmeta import RunUnitInfo


@dataclass
class SampleAnalysis:
    """Everything one admixture sample produces.

    Attributes:
        genotypes: Classified markers (informative subset, sex-chrom tally).
        result: The chimerism estimate, single- or multi-donor. The
            host-presence result is attached as ``result.host_presence`` when
            ``run_host_presence`` was True.
        qc: Quality-control verdict for ``result``.
        donor_hom_markers: Donor-homozygous host-presence markers with their
            artifact flags, the per-marker detail the diagnostics plot. Empty
            when ``run_host_presence`` is False.
    """

    genotypes: MarkerGenotypes
    result: ChimerismResult | MultiDonorResult
    qc: QCReport
    donor_hom_markers: list[DonorHomMarker]


def _floor_detection_limits(
    result: ChimerismResult | MultiDonorResult, contamination_fraction: float
) -> None:
    """Raise the reported LoB/LoD to the in-data contamination floor in place.

    The analytical ``lob_fraction`` / ``lod_fraction`` come from sequencing error
    and Fisher information alone (see ``allomix.chimerism.detection_limit``). They
    do not know about the co-pooled contamination floor allomix estimates per
    sample, which is a second noise term competing with sub-1% host detection: at
    a 0.2% floor a 0.5% host call is within noise, and the reported limit of
    detection should say so. This floors both limits at the contamination
    fraction so a sample whose contamination floor exceeds its analytical LoD
    reports the floor. See ``claude/further_improvements.md``, Obs 2.

    A no-op when ``contamination_fraction`` is 0 (a clean sample) or when the
    result carries no LoB/LoD fields (multi-donor results do not). Monotonicity
    ``lod >= lob`` is preserved because both are floored at the same value.

    Args:
        result: Chimerism result to update in place.
        contamination_fraction: In-data contamination floor (donor-fraction
            scale, 0.0-1.0) from ``allomix.contamination.estimate_contamination``.
    """
    if contamination_fraction <= 0.0 or not hasattr(result, "lob_fraction"):
        return
    result.lob_fraction = max(result.lob_fraction, contamination_fraction)
    result.lod_fraction = max(result.lod_fraction, contamination_fraction)


def analyse_sample(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    *,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    calibration: PanelCalibration | None = None,
    run_host_presence: bool = True,
    use_sex_chroms: bool = False,
    artifact_filter: bool = True,
    sample_name: str | None = None,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
    marker_type_overdispersion: bool = True,
    expected_relatedness: list[str] | None = None,
    relatedness_tolerance: int = 1,
    run_unit: RunUnitInfo | None = None,
) -> SampleAnalysis:
    """Run the chimerism pipeline for one pre-parsed admixture sample.

    Single-donor estimation is used when ``donors`` has one entry, multi-donor
    otherwise. The host-presence detector (on by default) is cheap and
    complementary to the MLE; see ``allomix.detect``.

    Args:
        host: Parsed host markers.
        donors: One parsed marker list per donor.
        admix: Parsed admixture markers (parse with ``min_dp=0``; filtering is
            applied here via ``min_dp``).
        min_dp: Minimum admix depth for a marker to be used.
        min_gq: Minimum host/donor genotype quality.
        error_rate: Global symmetric sequencing error rate.
        calibration: Optional per-marker bias and error tables (see
            ``allomix.likelihood.PanelCalibration``).
        run_host_presence: Run the host-presence detector and select the
            donor-homozygous markers. When False, ``result.host_presence`` is
            left unset and ``donor_hom_markers`` is empty.
        use_sex_chroms: Keep sex-chromosome markers (default drops them; their
            allele dosage is wrong in sex-mismatched transplants).
        artifact_filter: Drop alignment-artifact markers from the
            host-presence test (the returned ``donor_hom_markers`` still lists
            them, flagged).
        sample_name: Optional name stamped onto ``genotypes.sample_name``.
        robust: Robust-refit mode passed to the estimator ("off"/"auto"/"force";
            see ``estimate_single_donor_bb``). Drops host copy-number/LoH-
            inconsistent markers and refits; "auto" is the recommended policy.
        robust_k: Robust residual cut (robust SDs) for the refit.
        marker_type_overdispersion: Fit a separate beta-binomial rho per marker
            class (donor-hom vs donor-het) in single-donor estimation (issue #33).
            On by default; set False for the legacy shared-rho path. Ignored for
            multi-donor (later phase).
        expected_relatedness: Optional declared relationship per donor (one entry
            per ``donors`` list, value in ``allomix.relatedness.VALID_DECLARATIONS``
            or "NA"/None for no expectation). Compared against the estimated
            host-vs-donor relatedness in QC.
        relatedness_tolerance: Allowed degree distance before a declared-vs-detected
            relatedness mismatch is flagged (default 1; see ``evaluate_expected``).

    Returns:
        A ``SampleAnalysis`` bundling the genotypes, estimate, QC and the
        flagged donor-homozygous markers.
    """
    cal = calibration or PanelCalibration()
    genotypes = classify_markers(
        host, donors, admix, min_dp=min_dp, min_gq=min_gq, use_sex_chroms=use_sex_chroms
    )
    if sample_name is not None:
        genotypes.sample_name = sample_name

    if len(donors) == 1:
        result: ChimerismResult | MultiDonorResult = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            calibration=cal,
            robust=robust,
            robust_k=robust_k,
            marker_type_overdispersion=marker_type_overdispersion,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            calibration=cal,
            robust=robust,
            robust_k=robust_k,
        )

    # In-data contamination estimate at consensus-homozygous markers. Independent
    # of the MLE and of any sequencing-run metadata (issue #12). Computed before
    # the host-presence test and the LoD flooring below so its floor can feed
    # both: a co-pooled genome is a second noise term competing with sub-1% host
    # detection (see ``claude/further_improvements.md``, Obs 2). QC reads it off
    # the result like the other identity checks.
    result.contamination = estimate_contamination(
        host,
        donors,
        admix,
        marker_errors=cal.errors,
        error_rate=error_rate,
        min_dp=min_dp,
    )
    contamination_floor = result.contamination.contamination_fraction
    _floor_detection_limits(result, contamination_floor)

    dh_markers: list[DonorHomMarker] = []
    if run_host_presence:
        # Attached to the result before QC so the QC step can read it. The
        # contamination floor is added to the per-marker H0 background: a
        # co-pooled genome carrying the donor-absent allele inflates exactly the
        # counts this test reads, so raising the background guards against
        # calling that contamination as host signal.
        result.host_presence = host_presence_test(
            genotypes.informative,
            marker_errors=cal.errors,
            error_rate=error_rate,
            contamination_floor=contamination_floor,
            artifact_filter=artifact_filter,
        )
        dh_markers = donor_hom_markers(genotypes.informative)

    # Identity QC over the raw reference/admix markers (not the informative set,
    # which excludes the shared and consensus-homozygous sites these checks need).
    # Ordering invariant: host-vs-donor pairs first, in donor order, so QC can
    # align them with ``expected_relatedness``; donor-vs-donor pairs follow.
    donor_labels = ["donor"] if len(donors) == 1 else [f"donor{i + 1}" for i in range(len(donors))]
    relatedness: list[RelatednessResult] = [
        relatedness_coefficient(host, donors[i], "host", donor_labels[i])
        for i in range(len(donors))
    ]
    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            relatedness.append(
                relatedness_coefficient(donors[i], donors[j], donor_labels[i], donor_labels[j])
            )
    result.relatedness = relatedness
    result.admix_consistency = admix_consistency(
        host, donors, admix, error_rate=error_rate, min_dp=min_dp
    )
    # Run-unit metadata (index-hopping provenance) read from the admix VCF header
    # by the caller; attached before QC so the shared-run flag can be reported.
    result.run_unit = run_unit

    qc = assess_quality(
        result,
        genotypes,
        expected_relatedness=expected_relatedness,
        relatedness_tolerance=relatedness_tolerance,
    )

    return SampleAnalysis(
        genotypes=genotypes,
        result=result,
        qc=qc,
        donor_hom_markers=dh_markers,
    )


__all__ = ["SampleAnalysis", "analyse_sample"]
