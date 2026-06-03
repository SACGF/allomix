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

from allomix.chimerism import (
    ChimerismResult,
    MultiDonorResult,
    PanelCalibration,
    estimate_multi_donor,
    estimate_single_donor_bb,
)
from allomix.constants import ROBUST_K_DEFAULT
from allomix.detect import DonorHomMarker, donor_hom_markers, host_presence_test
from allomix.genotype import MarkerData, MarkerGenotypes, classify_markers
from allomix.qc import QCReport, assess_quality
from allomix.relatedness import (
    RelatednessResult,
    admix_consistency,
    relatedness_coefficient,
)


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
    expected_relatedness: list[str] | None = None,
    relatedness_tolerance: int = 1,
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
            ``allomix.chimerism.PanelCalibration``).
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

    dh_markers: list[DonorHomMarker] = []
    if run_host_presence:
        # Attached to the result before QC so the QC step can read it.
        result.host_presence = host_presence_test(
            genotypes.informative,
            marker_errors=cal.errors,
            error_rate=error_rate,
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
