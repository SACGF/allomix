"""VCF parsing and marker genotype classification for chimerism analysis.

Reads per-sample VCFs using cyvcf2, extracts genotype and allele depth
information at each marker, and classifies markers as informative or
non-informative based on host vs donor genotype comparison.
"""

from dataclasses import dataclass, field
from pathlib import Path

from cyvcf2 import VCF


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
    # artifact filter in ``allomix.detect``.
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
    # the host-presence artifact filter in ``allomix.detect``.
    admix_dp4: tuple[int, int, int, int] | None = None
    admix_rpbz: float | None = None
    admix_scbz: float | None = None
    admix_bqbz: float | None = None


@dataclass
class MarkerCounts:
    """Per-input marker counts explaining how the informative set was reached.

    Diagnostic data only: produced as a byproduct of classification and
    interpreted by ``allomix.qc``; the estimator never reads it.
    """

    n_host: int = 0  # markers genotyped in host
    n_donor_markers: list[int] = field(default_factory=list)  # per donor
    n_admix: int = 0  # markers in the admixture sample (the universe)
    n_admix_in_host: int = 0  # admix markers also genotyped in host
    n_admix_in_donor: list[int] = field(default_factory=list)  # admix markers also in each donor
    # Reasons shared markers were dropped before classification.
    n_drop_pass: int = 0  # failed PASS filter
    n_drop_gq_host: int = 0  # host GQ below threshold
    n_drop_gq_donor: int = 0  # a donor GQ below threshold
    n_drop_admix_dp: int = 0  # admixture depth below threshold


@dataclass
class MarkerGenotypes:
    """Result of parsing and classifying markers across host, donor(s), and admixture."""

    informative: list[InformativeMarker]
    non_informative: list[MarkerData]
    n_total: int
    n_shared: int
    n_filtered: int
    sample_name: str = ""
    marker_counts: MarkerCounts | None = None  # per-input diagnostic counts (see MarkerCounts)
    n_sex_chrom_excluded: int = 0  # informative sex-chrom markers dropped (use_sex_chroms=False)


_SEX_CHROM_NAMES = {"X", "Y", "M", "MT"}

# Reference-sample GT/AD consistency thresholds (see parse_vcf, gt_ad_consistency).
# A called genotype is dropped when its AD-derived VAF contradicts the call.
# Bounds are deliberately loose to tolerate genuine capture bias. The het and
# hom bounds are symmetric about 0.5, so each pair is defined from one number.
VAF_CHECK_MIN_DEPTH = 20  # need this many total reads before judging the VAF
HET_MIN_VAF = 0.35
HET_MAX_VAF = 1.0 - HET_MIN_VAF  # 0.65
HOM_REF_MAX_VAF = 0.05
HOM_ALT_MIN_VAF = 1.0 - HOM_REF_MAX_VAF  # 0.95


def is_sex_chrom(chrom: str) -> bool:
    """True if ``chrom`` is a sex or mitochondrial contig (chr-prefix optional)."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return c.upper() in _SEX_CHROM_NAMES


def parse_vcf(
    path: Path | str,
    sample: str | int = 0,
    min_dp: int = 0,
    min_gq: int = 0,
    gt_ad_consistency: bool = False,
) -> list[MarkerData]:
    """Read a VCF and extract MarkerData for a specific sample.

    Args:
        path: Path to VCF or VCF.gz file.
        sample: Sample name (str) or column index (int, 0-based). Default 0.
        min_dp: Minimum depth filter. Records below this are excluded.
        min_gq: Minimum genotype quality. Records below this are excluded.
        gt_ad_consistency: If True, drop markers where the called GT
            contradicts the AD-derived VAF (het outside [0.35, 0.65],
            hom-ref VAF > 0.05, hom-alt VAF < 0.95). Use for reference
            samples (host/donor) whose GT must match their reads. Do
            NOT use for admix samples — admix is a mixture, so its VAF
            is not expected to land at 0/0.5/1 by definition.

    Returns:
        List of MarkerData, one per passing record.

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
        # Skip multiallelic sites
        if len(variant.ALT) > 1:
            continue

        alt = variant.ALT[0] if variant.ALT else "."

        # Skip indels — admix-side AD comes from straight pileup which
        # cannot count indel reads the way local-reassembly callers
        # (GATK HaplotypeCaller) do, producing systematic admix=0-ALT
        # at sites where the panel sample is genuinely het/hom-alt.
        if alt != "." and (len(variant.REF) != 1 or len(alt) != 1):
            continue

        # Extract genotype for the selected sample
        gt_arr = variant.genotypes[sample_idx]  # [allele1, allele2, phased]
        a1, a2 = gt_arr[0], gt_arr[1]

        # Skip no-calls
        if a1 < 0 or a2 < 0:
            continue

        gt = (min(a1, a2), max(a1, a2))

        # Extract AD
        ad = variant.format("AD")
        if ad is not None:
            ad_vals = ad[sample_idx]
            ad_ref = int(ad_vals[0]) if ad_vals[0] >= 0 else 0
            ad_alt = int(ad_vals[1]) if len(ad_vals) > 1 and ad_vals[1] >= 0 else 0
        else:
            continue  # AD is required

        # Extract DP
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

        # Apply filters
        if dp < min_dp:
            continue
        if gq is not None and min_gq > 0 and gq < min_gq:
            continue

        # Reference-sample GT/AD consistency check. Caller-only filter:
        # admix is a mixture so its VAF is not expected to track its GT.
        # For host/donor a het call with VAF below 35% or above 65% is
        # almost certainly a miscall (GATK rescued it from marginal
        # evidence in a 2-sample joint call), and using it would feed a
        # systematic bias into the chimerism estimator at the recovered
        # marker. The thresholds are loose to tolerate genuine capture
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

        # Filter status
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
#: whole package: the bias table (``allomix.bias``), the error table
#: (``allomix.error_rates``), and the host-presence detector (``allomix.detect``)
#: all key markers by this tuple and import it from here.
MarkerKey = tuple[str, int, str, str]


def marker_key(m: MarkerData) -> MarkerKey:
    """Key for joining the same marker across samples and tables."""
    return (m.chrom, m.pos, m.ref, m.alt)


def _alt_dose(gt: tuple[int, int]) -> int:
    """Count of ALT alleles in a diploid genotype (0, 1, or 2)."""
    return gt[0] + gt[1]


def marker_type(host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> int | None:
    """Classify marker informativeness using Vynck et al. types.

    Types:
        0:  host hom-ref (0/0), donor hom-alt (1/1) — fully informative
        1:  host hom-alt (1/1), donor hom-ref (0/0) — fully informative
        10: host het (0/1), donor hom-ref (0/0) — partially informative
        11: host het (0/1), donor hom-alt (1/1) — partially informative
        20: host hom-ref (0/0), donor het (0/1) — partially informative
        21: host hom-alt (1/1), donor het (0/1) — partially informative

    Returns:
        Integer type code, or None if the marker is non-informative
        (same genotype in host and donor).
    """
    h_dose = _alt_dose(host_gt)
    d_dose = _alt_dose(donor_gt)

    if h_dose == d_dose:
        return None  # non-informative

    if h_dose == 0 and d_dose == 2:
        return 0
    if h_dose == 2 and d_dose == 0:
        return 1
    if h_dose == 1 and d_dose == 0:
        return 10
    if h_dose == 1 and d_dose == 2:
        return 11
    if h_dose == 0 and d_dose == 1:
        return 20
    if h_dose == 2 and d_dose == 1:
        return 21

    return None  # shouldn't reach here for biallelic diploid


def classify_markers(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admixture: list[MarkerData],
    min_dp: int = 100,
    min_gq: int = 20,
    pass_only: bool = True,
    sample_name: str = "",
    use_sex_chroms: bool = False,
) -> MarkerGenotypes:
    """Classify shared markers as informative or non-informative.

    Joins markers across host, donor(s), and admixture by (chrom, pos, ref, alt).
    Applies depth and quality filters. Assigns Vynck marker types for the first
    donor (multi-donor types are stored in donor_gts list).

    Args:
        host: Markers from host genotyping VCF.
        donors: List of marker lists, one per donor.
        admixture: Markers from post-HSCT admixture VCF.
        min_dp: Minimum depth for admixture sample at a marker.
        min_gq: Minimum GQ for host/donor genotyping samples.
        pass_only: Only use PASS-filtered markers.
        use_sex_chroms: If False (default), markers on the sex and mitochondrial
            contigs (X/Y/M) are excluded. They are unreliable for chimerism in
            sex-mismatched transplants, where the host/donor allele dosage on
            chrX/chrY is wrong. Enable per run only once host and donor sex are
            known to match. The count of informative sex-chrom markers dropped is
            reported in ``MarkerGenotypes.n_sex_chrom_excluded``.

    Returns:
        MarkerGenotypes with classified markers.
    """
    n_total = len(admixture)

    # Index all inputs by key
    host_idx = {marker_key(m): m for m in host}
    donor_idxs = [{marker_key(m): m for m in d} for d in donors]
    admix_idx = {marker_key(m): m for m in admixture}

    # Per-input coverage of the admixture marker set (the universe), to show
    # which input genotyping is sparse.
    admix_keys = set(admix_idx.keys())
    n_admix_in_host = len(admix_keys & set(host_idx.keys()))
    n_admix_in_donor = [len(admix_keys & set(di.keys())) for di in donor_idxs]

    # Find shared keys (present in host, all donors, and admixture)
    shared_keys = set(host_idx.keys()) & admix_keys
    for di in donor_idxs:
        shared_keys &= set(di.keys())

    n_shared = len(shared_keys)

    informative: list[InformativeMarker] = []
    non_informative: list[MarkerData] = []
    n_drop_pass = 0
    n_drop_gq_host = 0
    n_drop_gq_donor = 0
    n_drop_admix_dp = 0
    n_sex_chrom_excluded = 0

    for key in sorted(shared_keys):
        h = host_idx[key]
        ds = [di[key] for di in donor_idxs]
        a = admix_idx[key]

        # Filter: PASS only
        if pass_only and (
            h.filter != "PASS" or a.filter != "PASS" or any(d.filter != "PASS" for d in ds)
        ):
            n_drop_pass += 1
            continue

        # Filter: host/donor GQ
        if min_gq > 0:
            if h.gq is not None and h.gq < min_gq:
                n_drop_gq_host += 1
                continue
            if any(d.gq is not None and d.gq < min_gq for d in ds):
                n_drop_gq_donor += 1
                continue

        # Filter: admixture depth
        if a.dp < min_dp:
            n_drop_admix_dp += 1
            continue

        # Classify: informative if host differs from ANY donor
        donor_gts = [d.gt for d in ds]
        mtypes = [marker_type(h.gt, d.gt) for d in ds]
        any_informative = any(mt is not None for mt in mtypes)

        # Drop sex / mitochondrial contigs unless explicitly enabled, counting
        # the informative ones lost so the cost is visible.
        if not use_sex_chroms and is_sex_chrom(key[0]):
            if any_informative:
                n_sex_chrom_excluded += 1
            continue

        if any_informative:
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
        n_drop_pass=n_drop_pass,
        n_drop_gq_host=n_drop_gq_host,
        n_drop_gq_donor=n_drop_gq_donor,
        n_drop_admix_dp=n_drop_admix_dp,
    )

    return MarkerGenotypes(
        informative=informative,
        non_informative=non_informative,
        n_total=n_total,
        n_shared=n_shared,
        n_filtered=n_drop_pass + n_drop_gq_host + n_drop_gq_donor + n_drop_admix_dp,
        sample_name=sample_name,
        marker_counts=counts,
        n_sex_chrom_excluded=n_sex_chrom_excluded,
    )
