"""Single-sample analysis pipeline shared by the CLI and diagnostic scripts.

``analyse_sample`` runs the full per-sample path once (classify markers, estimate
chimerism, run the host-presence detector, assess QC, select the donor-homozygous
markers) so that path and every knob live in one place for both ``allomix.cli``
and the ``scripts/`` diagnostics.

Library code: it does not read VCFs (callers parse first, conventions differ) and
it does not print. Callers own I/O and messaging.
"""

from dataclasses import dataclass

from allomix.constants import ROBUST_K_DEFAULT, SEX_MIN_REF_DP
from allomix.estimate.chimerism import estimate_multi_donor, estimate_single_donor_bb
from allomix.estimate.likelihood import PanelCalibration
from allomix.genotype import ContigPolicy, MarkerData, MarkerGenotypes, classify_markers
from allomix.qc.host_presence import DonorHomMarker, donor_hom_markers, host_presence_test
from allomix.qc.qc import QCReport, assess_quality
from allomix.qc.relatedness import (
    Relatedness,
    RelatednessResult,
    admix_consistency,
    relatedness_coefficient,
    shared_het_balance,
)
from allomix.qc.runmeta import RunUnitInfo
from allomix.qc.sample_contamination import estimate_contamination
from allomix.qc.sex import Sex, assess_sex
from allomix.results import ChimerismResult, MultiDonorResult


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

    Analytical limits (see ``allomix.estimate.chimerism.detection_limit``) come from
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
    contig_policy: ContigPolicy = ContigPolicy.SEX_AWARE,
    declared_host_sex: Sex | None = None,
    declared_donor_sexes: list[Sex | None] | None = None,
    artifact_filter: bool = True,
    sample_name: str | None = None,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
    marker_type_overdispersion: bool = True,
    expected_relatedness: list[Relatedness | None] | None = None,
    relatedness_tolerance: int = 1,
    run_unit: RunUnitInfo | None = None,
    clinical_gating: bool = True,
) -> AdmixtureSampleAnalysis:
    """Run the chimerism pipeline for one pre-parsed admixture sample.

    Single-donor estimation when ``donors`` has one entry, multi-donor otherwise.
    The host-presence detector (on by default) is cheap and complementary to the
    MLE; see ``allomix.qc.host_presence``.

    Args:
        admix: Parsed admixture markers (parse with ``min_dp=0``; filtering is
            applied here via ``min_dp``).
        run_host_presence: When False, ``result.host_presence`` is left unset and
            ``donor_hom_markers`` is empty.
        contig_policy: Which contig classes enter the informative set (see
            ``ContigPolicy``). The default ``SEX_AWARE`` admits autosomes and
            routes non-PAR chrX on the inferred/declared host-donor sex pair
            (#46); ``AUTOSOMES_ONLY`` never uses chrX; ``ALL_PRIMARY`` is for
            diagnostic genomic views.
        declared_host_sex: Declared host sex (``Sex.FEMALE`` / ``Sex.MALE``) or
            None. Compared with the chrX-inferred sex in QC; see ``allomix.qc.sex``.
        declared_donor_sexes: Declared sex per donor, aligned with ``donors``
            (None entries for no declaration); None for nothing declared.
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
        clinical_gating: When True (default), gate the reliability REVIEW flags on
            effect size, not bare significance: promote a goodness-of-fit misfit only
            when the reduced chi-sq is large, the pre-trim full-set guard only near
            the detection floor, and the consensus-hom swap test only when the
            discordant fraction is high. When False, use the legacy p-value-only rules.
    """
    cal = calibration or PanelCalibration()
    # Sex of each reference sample from its own non-PAR chrX heterozygosity,
    # reconciled with any declared sex, plus the host/donor pair status (#50).
    # Computed before classification because the ``SEX_AWARE`` contig policy
    # routes non-PAR chrX markers on the pair status (#46), and attached to the
    # result before QC, which fails a confident inference that contradicts the
    # declaration. The depth floor is the reference-sample one, not the
    # admixture ``min_dp`` (see ``SEX_MIN_REF_DP``).
    sex = assess_sex(
        host,
        donors,
        min_dp=SEX_MIN_REF_DP,
        min_gq=min_gq,
        declared_host=declared_host_sex,
        declared_donors=declared_donor_sexes,
    )
    genotypes = classify_markers(
        host,
        donors,
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        contig_policy=contig_policy,
        pair_status=sex.pair,
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
    # Consensus-het allele balance: orthogonal to the consensus-hom swap and
    # contamination checks above (issue #38). Flags contamination, CNV/allelic
    # imbalance, or a sample mix-up via VAF skew at sites het in all parties.
    result.shared_het_balance = shared_het_balance(host, donors, admix, min_dp=min_dp)
    result.sex = sex
    # Run-unit metadata (index-hopping provenance); attached before QC so the
    # shared-run flag can be reported.
    result.run_unit = run_unit

    qc = assess_quality(
        result,
        genotypes,
        expected_relatedness=expected_relatedness,
        relatedness_tolerance=relatedness_tolerance,
        clinical_gating=clinical_gating,
    )

    return AdmixtureSampleAnalysis(
        genotypes=genotypes,
        result=result,
        qc=qc,
        donor_hom_markers=dh_markers,
    )


__all__ = ["AdmixtureSampleAnalysis", "analyse_sample"]
