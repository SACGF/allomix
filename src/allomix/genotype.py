"""VCF parsing and marker genotype classification for chimerism analysis.

Reads per-sample VCFs using cyvcf2, extracts genotype and allele depth
information at each marker, and classifies markers as informative or
non-informative based on host vs donor genotype comparison.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path

from cyvcf2 import VCF

from allomix.constants import (
    DEFAULT_MIN_DP,
    DEFAULT_MIN_GQ,
    HOM_ALT_MIN_VAF,
    HOM_REF_MAX_VAF,
)

# ``is_sex_chrom`` moved to ``allomix.contigs``; re-exported here so existing
# ``from allomix.genotype import is_sex_chrom`` callers keep working.
from allomix.contigs import ContigClass, classify_contig, is_sex_chrom  # noqa: F401
from allomix.sex_types import PairStatus


class ContigPolicy(Enum):
    """Which contig classes ``classify_markers`` admits to the informative set.

    PAR X/Y sites and ``OTHER`` (non-primary) contigs are excluded under every
    policy; the policy only decides what happens to non-PAR chrX, non-PAR chrY
    and MT.

    - ``SEX_AWARE`` (default): autosomes always; non-PAR chrX routed on the
      host/donor sex pair status from ``allomix.qc.sex`` (#46). For a
      ``MATCHED_FEMALE`` pair chrX markers go through the normal diploid path.
      For a ``MATCHED_MALE`` pair they are admitted only where the host and
      every donor are homozygous (a het call on a hemizygous chromosome is a
      genotyping error; those markers are dropped and counted in
      ``MarkerGenotypes.n_chrx_male_het_dropped``). For a ``MISMATCHED`` or
      ``UNKNOWN`` pair, or when no pair status is given, chrX is excluded.
      Non-PAR chrY and MT are never used in the estimate. Informative markers
      excluded this way are counted in
      ``MarkerGenotypes.n_informative_sex_chrom_excluded``.
    - ``AUTOSOMES_ONLY``: autosomes only, whatever the pair status. Informative
      non-PAR X/Y/MT markers are dropped and counted as above.
    - ``ALL_PRIMARY``: diagnostic only. Also admits non-PAR chrX, non-PAR chrY
      and MT so genomic views can show them. The diploid dosage model is wrong
      on those contigs for a sex-mismatched pair, so this is never the clinical
      path.
    """

    SEX_AWARE = "sex_aware"
    AUTOSOMES_ONLY = "autosomes_only"
    ALL_PRIMARY = "all_primary"


class MarkerType(IntEnum):
    """Vynck host-vs-donor informative marker class.

    The single definition of the Vynck et al. classification used across allomix.
    The member value is the type code stored on results and serialized output, so
    a ``MarkerType`` compares and serialises as its code (it is an ``int``). Each
    member also carries the host and donor ALT dose it represents
    (``host_dose``/``donor_dose``); those drive both the readable ``label`` and the
    ``classify`` reverse lookup, so the host/donor genotype pairing is defined once
    on the members rather than restated in a parallel table.

    Reference: Vynck et al., classification of bi-allelic markers by host/donor
    genotype contrast for chimerism estimation.
    """

    # name = (code, host ALT dose, donor ALT dose)
    HOST_HOMREF_DONOR_HOMALT = (0, 0, 2)
    HOST_HOMALT_DONOR_HOMREF = (1, 2, 0)
    HOST_HET_DONOR_HOMREF = (10, 1, 0)
    HOST_HET_DONOR_HOMALT = (11, 1, 2)
    HOST_HOMREF_DONOR_HET = (20, 0, 1)
    HOST_HOMALT_DONOR_HET = (21, 2, 1)

    def __new__(cls, code: int, host_dose: int, donor_dose: int) -> MarkerType:
        member = int.__new__(cls, code)
        member._value_ = code
        member.host_dose = host_dose
        member.donor_dose = donor_dose
        return member

    @property
    def label(self) -> str:
        """Readable ``host G/G, donor G/G`` genotype label for this type."""
        return f"host {_DOSE_GT[self.host_dose]}, donor {_DOSE_GT[self.donor_dose]}"

    @classmethod
    def classify(cls, host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> MarkerType | None:
        """Classify marker informativeness from the host and donor genotypes.

        Returns the matching member, or None when the marker is non-informative
        (host and donor carry the same ALT dose). Allele frequencies are not
        needed, only the ALT dose contrast, which suits an arbitrary panel.
        """
        key = (host_gt[0] + host_gt[1], donor_gt[0] + donor_gt[1])
        return _MARKER_TYPE_BY_DOSE.get(key)

    @classmethod
    def label_for(cls, code: int) -> str:
        """Readable host/donor genotype label for a marker-type code.

        ``code`` may be a bare integer (as found in serialized output). Returns the
        raw integer as a string for an unrecognised code, so output never breaks on
        an unexpected value.
        """
        try:
            return cls(code).label
        except ValueError:
            return str(code)


# Diploid ALT dose -> displayed genotype, for MarkerType.label.
_DOSE_GT = {0: "0/0", 1: "0/1", 2: "1/1"}
# Reverse lookup (host_dose, donor_dose) -> MarkerType, derived from the members so
# the genotype pairing has a single source of truth. Non-informative pairs (equal
# dose) are simply absent, so classify() returns None for them.
_MARKER_TYPE_BY_DOSE: dict[tuple[int, int], MarkerType] = {
    (m.host_dose, m.donor_dose): m for m in MarkerType
}


@dataclass
class MarkerData:
    """Genotype and depth data at a single marker for one sample."""

    chrom: str
    pos: int
    ref: str
    alt: str
    gt: tuple[int, int]  # allele indices, e.g. (0,0), (0,1), (1,1)
    ad_ref: int
    ad_alt: int
    dp: int
    gq: int | None = None
    filter: str = "PASS"
    # Read-level bias annotations from bcftools mpileup INFO. Present on admix
    # samples; None for GATK-called panel samples (which do not emit them) and
    # at sites without contrasting reads. Consumed by the host-presence
    # artifact filter in ``allomix.qc.host_presence``.
    dp4: tuple[int, int, int, int] | None = None  # ref_fwd, ref_rev, alt_fwd, alt_rev
    rpbz: float | None = None  # read-position bias z-score (|.| large => artifact)
    scbz: float | None = None  # soft-clip-length bias z-score
    bqbz: float | None = None  # base-quality bias z-score


@dataclass
class InformativeMarker:
    """A marker where host and at least one donor have different genotypes."""

    chrom: str
    pos: int
    ref: str
    alt: str
    host_gt: tuple[int, int]
    donor_gts: list[tuple[int, int]]
    marker_type: int  # Vynck classification for first donor: 0,1,10,11,20,21
    admix_ad_ref: int
    admix_ad_alt: int
    admix_dp: int
    marker_types: list[int | None] | None = None  # Vynck type per donor (None = non-informative)
    informative_for: list[bool] | None = None  # True per donor if informative
    # Admix-side read-level bias (copied from the admix MarkerData), passed to
    # the host-presence artifact filter in ``allomix.qc.host_presence``.
    admix_dp4: tuple[int, int, int, int] | None = None
    admix_rpbz: float | None = None
    admix_scbz: float | None = None
    admix_bqbz: float | None = None


@dataclass
class MarkerCounts:
    """Per-input marker counts explaining how the informative set was reached.

    Diagnostic data only: produced as a byproduct of classification and
    interpreted by ``allomix.qc.qc``; the estimator never reads it.
    """

    n_host: int = 0
    n_donor_markers: list[int] = field(default_factory=list)  # per donor
    n_admix: int = 0  # the universe
    n_admix_in_host: int = 0
    n_admix_in_donor: list[int] = field(default_factory=list)  # per donor
    n_dropped_failed_filter: int = 0
    n_dropped_low_gq_host: int = 0
    n_dropped_low_gq_donor: int = 0
    n_dropped_low_admix_dp: int = 0


@dataclass
class MarkerGenotypes:
    """Result of parsing and classifying markers across host, donor(s), and admixture."""

    informative: list[InformativeMarker]
    non_informative: list[MarkerData]
    n_total: int
    n_shared: int
    n_filtered: int
    sample_name: str = ""
    marker_counts: MarkerCounts | None = None  # per-input diagnostic counts
    # Informative non-PAR X/Y/MT markers dropped by the contig policy: all of
    # them under ``AUTOSOMES_ONLY``; under ``SEX_AWARE`` the non-PAR chrY and MT
    # ones plus chrX for a mismatched or unknown pair.
    n_informative_sex_chrom_excluded: int = 0
    # Shared markers excluded by contig class under every ``ContigPolicy``:
    # pseudoautosomal X/Y sites, and sites on non-primary contigs (alt, decoy,
    # unplaced, random). See ``allomix.contigs``.
    n_par_excluded: int = 0
    n_other_contig_excluded: int = 0
    # Non-PAR chrX markers admitted to the informative set: a sex-matched pair
    # under ``SEX_AWARE`` (#46), or anything under ``ALL_PRIMARY``. Always 0
    # under ``AUTOSOMES_ONLY``.
    n_chrx_used: int = 0
    # Non-PAR chrX markers dropped from a male/male pair because the host or a
    # donor carries a het call there (a genotyping error on a hemizygous
    # chromosome). Counted informative or not; only ``SEX_AWARE`` produces them.
    n_chrx_male_het_dropped: int = 0


# Reference-sample GT/AD consistency thresholds (see parse_vcf, gt_ad_consistency).
# A called genotype is dropped when its AD-derived VAF contradicts the call.
# Bounds are deliberately loose to tolerate genuine capture bias. The het bounds
# are symmetric about 0.5; the hom cutoffs (shared with the simulator's caller)
# live in allomix.constants.
VAF_CHECK_MIN_DEPTH = 20  # need this many total reads before judging the VAF
HET_MIN_VAF = 0.35
HET_MAX_VAF = 1.0 - HET_MIN_VAF  # 0.65


def parse_vcf(
    path: Path | str,
    sample: str | int = 0,
    min_dp: int = 0,
    min_gq: int = 0,
    gt_ad_consistency: bool = False,
    include_sites: set[tuple[str, int]] | None = None,
    exclude_sites: set[tuple[str, int]] | None = None,
) -> list[MarkerData]:
    """Read a VCF and extract MarkerData for a specific sample.

    Args:
        sample: Sample name (str) or 0-based column index (int). Default 0.
        gt_ad_consistency: If True, drop markers where the called GT contradicts
            the AD-derived VAF (het outside [0.35, 0.65], hom-ref VAF > 0.05,
            hom-alt VAF < 0.95). Use for reference samples (host/donor) whose GT
            must match their reads. Do NOT use for admix: a mixture's VAF is not
            expected to land at 0/0.5/1.
        include_sites: If given, keep only markers whose ``(chrom, pos)`` is in
            this set (1-based positions). Every other site is dropped.
        exclude_sites: If given, drop markers whose ``(chrom, pos)`` is in this
            set. Applied after ``include_sites``. Pass at most one of the two
            (the CLI enforces mutual exclusivity via a panel-qc sites BED).

    Raises:
        ValueError: If a string sample name is not found in the VCF.
    """
    markers: list[MarkerData] = []
    vcf = VCF(str(path))

    if isinstance(sample, str):
        if sample not in vcf.samples:
            raise ValueError(f"Sample '{sample}' not found in VCF. Available: {list(vcf.samples)}")
        sample_idx = list(vcf.samples).index(sample)
    else:
        sample_idx = sample

    for variant in vcf:
        if include_sites is not None and (variant.CHROM, variant.POS) not in include_sites:
            continue
        if exclude_sites is not None and (variant.CHROM, variant.POS) in exclude_sites:
            continue

        if len(variant.ALT) > 1:
            continue

        alt = variant.ALT[0] if variant.ALT else "."

        # Skip indels: admix-side AD comes from straight pileup, which cannot
        # count indel reads the way local-reassembly callers (GATK
        # HaplotypeCaller) do, giving systematic admix=0-ALT at sites where the
        # panel sample is genuinely het/hom-alt.
        if alt != "." and (len(variant.REF) != 1 or len(alt) != 1):
            continue

        gt_arr = variant.genotypes[sample_idx]  # [allele1, allele2, phased]
        a1, a2 = gt_arr[0], gt_arr[1]

        if a1 < 0 or a2 < 0:
            continue  # skip no-calls

        gt = (min(a1, a2), max(a1, a2))

        ad = variant.format("AD")
        if ad is not None:
            ad_vals = ad[sample_idx]
            ad_ref = int(ad_vals[0]) if ad_vals[0] >= 0 else 0
            ad_alt = int(ad_vals[1]) if len(ad_vals) > 1 and ad_vals[1] >= 0 else 0
        else:
            continue  # AD is required

        dp_arr = variant.format("DP")
        if dp_arr is not None:
            dp = int(dp_arr[sample_idx][0])
        else:
            dp = ad_ref + ad_alt

        # Extract GQ. cyvcf2 raises KeyError when GQ is not declared in the
        # VCF header (e.g. bcftools call -C alleles output), as opposed to
        # returning None for declared-but-missing — treat both as "no GQ".
        try:
            gq_arr = variant.format("GQ")
        except KeyError:
            gq_arr = None
        gq = int(gq_arr[sample_idx][0]) if gq_arr is not None else None

        if dp < min_dp:
            continue
        if gq is not None and min_gq > 0 and gq < min_gq:
            continue

        # Reference-sample GT/AD consistency check (host/donor only; an admix
        # mixture's VAF need not track its GT). A het call with VAF outside
        # 35-65% is almost certainly a GATK miscall rescued from marginal
        # evidence in a 2-sample joint call, and would feed a systematic bias
        # into the estimator. Thresholds are loose to tolerate genuine capture
        # bias (median |bias| ~0.5%, 95th pct ~4%).
        if gt_ad_consistency and (ad_ref + ad_alt) >= VAF_CHECK_MIN_DEPTH and alt != ".":
            vaf = ad_alt / (ad_ref + ad_alt)
            if gt == (0, 1):
                if vaf < HET_MIN_VAF or vaf > HET_MAX_VAF:
                    continue
            elif gt == (0, 0):
                if vaf > HOM_REF_MAX_VAF:
                    continue
            elif gt == (1, 1):
                if vaf < HOM_ALT_MIN_VAF:
                    continue

        filt = variant.FILTER
        if filt is None:
            filt = "PASS"

        # Read-level bias annotations (bcftools mpileup). Absent in GATK VCFs
        # and at sites without contrasting reads; missing -> None.
        info = variant.INFO
        dp4_raw = info.get("DP4")
        dp4: tuple[int, int, int, int] | None = None
        if dp4_raw is not None:
            dp4_vals = tuple(int(x) for x in dp4_raw)
            if len(dp4_vals) == 4:
                dp4 = dp4_vals  # type: ignore[assignment]
        rpbz = info.get("RPBZ")
        scbz = info.get("SCBZ")
        bqbz = info.get("BQBZ")

        markers.append(
            MarkerData(
                chrom=variant.CHROM,
                pos=variant.POS,
                ref=variant.REF,
                alt=alt,
                gt=gt,
                ad_ref=ad_ref,
                ad_alt=ad_alt,
                dp=dp,
                gq=gq,
                filter=filt,
                dp4=dp4,
                rpbz=float(rpbz) if rpbz is not None else None,
                scbz=float(scbz) if scbz is not None else None,
                bqbz=float(bqbz) if bqbz is not None else None,
            )
        )

    vcf.close()
    return markers


#: Marker key shape ``(chrom, pos, ref, alt)``. Canonical definition for the
#: whole package: the bias table (``allomix.calibration.bias``), the error table
#: (``allomix.calibration.error_rates``), and the host-presence detector
#: (``allomix.qc.host_presence``) all key markers by this tuple and import it from here.
MarkerKey = tuple[str, int, str, str]


def marker_key(m: MarkerData) -> MarkerKey:
    """Key for joining the same marker across samples and tables."""
    return (m.chrom, m.pos, m.ref, m.alt)


def _is_hom(gt: tuple[int, int]) -> bool:
    return gt[0] == gt[1]


def admit_contig(
    contig_class: ContigClass,
    policy: ContigPolicy,
    pair_status: PairStatus | None,
    gts: list[tuple[int, int]],
) -> tuple[bool, bool]:
    """Decide whether a non-PAR, primary-contig marker may enter the informative set.

    The routing table of the sex-marker plan (#46), for the contig classes that
    survive the unconditional PAR and non-primary exclusions. Pure, so it can
    be tested row by row.

    Args:
        contig_class: ``AUTOSOME``, ``X_NONPAR``, ``Y_NONPAR`` or ``MT``.
        policy: The ``ContigPolicy`` in force.
        pair_status: Host/donor sex pair status from ``allomix.qc.sex``, or None
            when sex was not assessed (treated like ``UNKNOWN``).
        gts: Host genotype followed by every donor genotype at the marker.

    Returns:
        ``(admit, male_het_drop)``. ``admit`` is True when the marker may be
        used. ``male_het_drop`` is True only for a non-PAR chrX marker of a
        ``MATCHED_MALE`` pair that is rejected because some party is
        heterozygous there.
    """
    if contig_class is ContigClass.AUTOSOME or policy is ContigPolicy.ALL_PRIMARY:
        return True, False
    if policy is ContigPolicy.AUTOSOMES_ONLY or contig_class is not ContigClass.X_NONPAR:
        return False, False
    # SEX_AWARE, non-PAR chrX.
    if pair_status is PairStatus.MATCHED_FEMALE:
        return True, False
    if pair_status is PairStatus.MATCHED_MALE:
        if all(_is_hom(gt) for gt in gts):
            return True, False
        return False, True
    return False, False


def classify_markers(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admixture: list[MarkerData],
    min_dp: int = DEFAULT_MIN_DP,
    min_gq: int = DEFAULT_MIN_GQ,
    pass_only: bool = True,
    sample_name: str = "",
    contig_policy: ContigPolicy = ContigPolicy.SEX_AWARE,
    pair_status: PairStatus | None = None,
) -> MarkerGenotypes:
    """Classify shared markers as informative or non-informative.

    Joins markers across host, donor(s), and admixture by (chrom, pos, ref, alt),
    applies depth and quality filters, and assigns Vynck marker types for the
    first donor (multi-donor types are stored in the donor_gts list).

    ``contig_policy`` decides which contig classes may enter the informative
    set (see ``ContigPolicy``); under the default ``SEX_AWARE`` policy
    ``pair_status`` (from ``allomix.qc.sex.assess_sex``) drives the non-PAR
    chrX routing, and None is treated as an unknown pair (chrX excluded). PAR
    and non-primary contigs are excluded under every policy.
    """
    n_total = len(admixture)

    host_idx = {marker_key(m): m for m in host}
    donor_idxs = [{marker_key(m): m for m in d} for d in donors]
    admix_idx = {marker_key(m): m for m in admixture}

    # Per-input coverage of the admixture marker set (the universe), to show
    # which input genotyping is sparse.
    admix_keys = set(admix_idx.keys())
    n_admix_in_host = len(admix_keys & set(host_idx.keys()))
    n_admix_in_donor = [len(admix_keys & set(di.keys())) for di in donor_idxs]

    # Shared keys: present in host, all donors, and admixture.
    shared_keys = set(host_idx.keys()) & admix_keys
    for di in donor_idxs:
        shared_keys &= set(di.keys())

    n_shared = len(shared_keys)

    informative: list[InformativeMarker] = []
    non_informative: list[MarkerData] = []
    n_dropped_failed_filter = 0
    n_dropped_low_gq_host = 0
    n_dropped_low_gq_donor = 0
    n_dropped_low_admix_dp = 0
    n_informative_sex_chrom_excluded = 0
    n_par_excluded = 0
    n_other_contig_excluded = 0
    n_chrx_used = 0
    n_chrx_male_het_dropped = 0

    for key in sorted(shared_keys):
        h = host_idx[key]
        ds = [di[key] for di in donor_idxs]
        a = admix_idx[key]

        # Contig-class exclusions that no option can override: non-primary
        # contigs are mapping hazards, and pseudoautosomal X/Y sites have
        # ambiguous copy number. Both are counted so the loss is visible.
        contig_class = classify_contig(key[0], key[1])
        if contig_class is ContigClass.OTHER:
            n_other_contig_excluded += 1
            continue
        if contig_class in (ContigClass.X_PAR, ContigClass.Y_PAR):
            n_par_excluded += 1
            continue

        if pass_only and (
            h.filter != "PASS" or a.filter != "PASS" or any(d.filter != "PASS" for d in ds)
        ):
            n_dropped_failed_filter += 1
            continue

        if min_gq > 0:
            if h.gq is not None and h.gq < min_gq:
                n_dropped_low_gq_host += 1
                continue
            if any(d.gq is not None and d.gq < min_gq for d in ds):
                n_dropped_low_gq_donor += 1
                continue

        if a.dp < min_dp:
            n_dropped_low_admix_dp += 1
            continue

        # Informative if the host differs from ANY donor.
        donor_gts = [d.gt for d in ds]
        mtypes = [MarkerType.classify(h.gt, d.gt) for d in ds]
        any_informative = any(mt is not None for mt in mtypes)

        # Non-PAR sex / mitochondrial contigs: route on the policy and the
        # host/donor sex pair (see ``admit_contig``), counting what is lost so
        # the cost is visible.
        admit, male_het_drop = admit_contig(
            contig_class, contig_policy, pair_status, [h.gt, *donor_gts]
        )
        if not admit:
            if male_het_drop:
                n_chrx_male_het_dropped += 1
            elif any_informative:
                n_informative_sex_chrom_excluded += 1
            continue

        if any_informative:
            if contig_class is ContigClass.X_NONPAR:
                n_chrx_used += 1
            # Use first donor's type for backward compat; fall back to first non-None
            mtype_first = mtypes[0]
            if mtype_first is None:
                mtype_first = next(mt for mt in mtypes if mt is not None)
            informative.append(
                InformativeMarker(
                    chrom=key[0],
                    pos=key[1],
                    ref=key[2],
                    alt=key[3],
                    host_gt=h.gt,
                    donor_gts=donor_gts,
                    marker_type=mtype_first,
                    admix_ad_ref=a.ad_ref,
                    admix_ad_alt=a.ad_alt,
                    admix_dp=a.dp,
                    marker_types=mtypes,
                    informative_for=[mt is not None for mt in mtypes],
                    admix_dp4=a.dp4,
                    admix_rpbz=a.rpbz,
                    admix_scbz=a.scbz,
                    admix_bqbz=a.bqbz,
                )
            )
        else:
            non_informative.append(a)

    counts = MarkerCounts(
        n_host=len(host),
        n_donor_markers=[len(d) for d in donors],
        n_admix=len(admixture),
        n_admix_in_host=n_admix_in_host,
        n_admix_in_donor=n_admix_in_donor,
        n_dropped_failed_filter=n_dropped_failed_filter,
        n_dropped_low_gq_host=n_dropped_low_gq_host,
        n_dropped_low_gq_donor=n_dropped_low_gq_donor,
        n_dropped_low_admix_dp=n_dropped_low_admix_dp,
    )

    return MarkerGenotypes(
        informative=informative,
        non_informative=non_informative,
        n_total=n_total,
        n_shared=n_shared,
        n_filtered=(
            n_dropped_failed_filter
            + n_dropped_low_gq_host
            + n_dropped_low_gq_donor
            + n_dropped_low_admix_dp
        ),
        sample_name=sample_name,
        marker_counts=counts,
        n_informative_sex_chrom_excluded=n_informative_sex_chrom_excluded,
        n_par_excluded=n_par_excluded,
        n_other_contig_excluded=n_other_contig_excluded,
        n_chrx_used=n_chrx_used,
        n_chrx_male_het_dropped=n_chrx_male_het_dropped,
    )
