"""Single-sample analysis pipeline shared by the CLI and diagnostic scripts.

``analyse_sample`` runs the full per-sample path once (classify markers, estimate
chimerism, run the host-presence detector, assess QC, select the donor-homozygous
markers) so that path and every knob live in one place for both ``allomix.cli``
and the ``scripts/`` diagnostics.

Library code: it does not read VCFs (callers parse first, conventions differ) and
it does not print. Callers own I/O and messaging.
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
    Relatedness,
    RelatednessResult,
    admix_consistency,
    relatedness_coefficient,
)
from allomix.results import ChimerismResult, MultiDonorResult
from allomix.runmeta import RunUnitInfo


@dataclass
class AdmixtureSampleAnalysis:
    genotypes: MarkerGenotypes
    result: ChimerismResult | MultiDonorResult  # host_presence attached when run_host_presence
    qc: QCReport
    donor_hom_markers: list[DonorHomMarker]  # empty when run_host_presence is False


def _floor_detection_limits(
    result: ChimerismResult | MultiDonorResult, contamination_fraction: float
) -> None:
    """Floor LoB/LoD at the in-data contamination level (further_improvements.md, Obs 2).

    Analytical limits (see ``allomix.chimerism.detection_limit``) come from
    sequencing error and Fisher information alone; they ignore the co-pooled
    contamination floor, a second noise term competing with sub-1% host detection.
    No-op when the fraction is 0 or the result has no LoB/LoD fields (multi-donor).
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
    expected_relatedness: list[Relatedness | None] | None = None,
    relatedness_tolerance: int = 1,
    run_unit: RunUnitInfo | None = None,
) -> AdmixtureSampleAnalysis:
    """Run the chimerism pipeline for one pre-parsed admixture sample.

    Single-donor estimation when ``donors`` has one entry, multi-donor otherwise.
    The host-presence detector (on by default) is cheap and complementary to the
    MLE; see ``allomix.detect``.

    Args:
        admix: Parsed admixture markers (parse with ``min_dp=0``; filtering is
            applied here via ``min_dp``).
        run_host_presence: When False, ``result.host_presence`` is left unset and
            ``donor_hom_markers`` is empty.
        artifact_filter: Drop alignment-artifact markers from the host-presence
            test (returned ``donor_hom_markers`` still lists them, flagged).
        robust: Robust-refit mode ("off"/"auto"/"force"; see
            ``estimate_single_donor_bb``). Drops host copy-number/LoH-inconsistent
            markers and refits; "auto" is the recommended policy.
        marker_type_overdispersion: Fit a separate beta-binomial rho per marker
            class (donor-hom vs donor-het) in single-donor estimation (issue #33).
            Ignored for multi-donor.
        expected_relatedness: Declared relationship per donor as a ``Relatedness``
            member (one entry per ``donors``; None for no expectation). Compared
            against estimated host-vs-donor relatedness in QC.
        relatedness_tolerance: Allowed degree distance before a declared-vs-detected
            mismatch is flagged (see ``evaluate_expected``).
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

    # In-data contamination estimate at consensus-homozygous markers, independent
    # of the MLE and run metadata (issue #12). Computed before the host-presence
    # test and LoD flooring so its floor feeds both (further_improvements.md, Obs 2).
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
        # Attached before QC so QC can read it. The contamination floor raises the
        # per-marker H0 background, guarding against calling a co-pooled genome's
        # donor-absent allele as host signal.
        result.host_presence = host_presence_test(
            genotypes.informative,
            marker_errors=cal.errors,
            error_rate=error_rate,
            contamination_floor=contamination_floor,
            artifact_filter=artifact_filter,
        )
        dh_markers = donor_hom_markers(genotypes.informative)

    # Identity QC over the raw reference/admix markers, not the informative set
    # (which excludes the shared and consensus-hom sites these checks need).
    # Ordering invariant: host-vs-donor pairs first in donor order (so QC aligns
    # them with ``expected_relatedness``), then donor-vs-donor pairs.
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
    # Run-unit metadata (index-hopping provenance); attached before QC so the
    # shared-run flag can be reported.
    result.run_unit = run_unit

    qc = assess_quality(
        result,
        genotypes,
        expected_relatedness=expected_relatedness,
        relatedness_tolerance=relatedness_tolerance,
    )

    return AdmixtureSampleAnalysis(
        genotypes=genotypes,
        result=result,
        qc=qc,
        donor_hom_markers=dh_markers,
    )


__all__ = ["AdmixtureSampleAnalysis", "analyse_sample"]
