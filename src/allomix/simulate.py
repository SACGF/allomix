"""Synthetic chimeric VCF generator for testing allomix.

Blends two genotype VCFs at a specified mixture fraction to produce a synthetic
chimeric VCF with realistic allele counts drawn from a binomial distribution.

Uses plain-text VCF parsing only (no cyvcf2 dependency) so this module can be
used in test environments without compiled libraries.
"""

import math
import random
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from allomix.constants import (
    DEFAULT_ERROR_RATE,
    HOM_ALT_MIN_VAF,
    HOM_REF_MAX_VAF,
    N_OTHER_BASES,
)
from allomix.estimate.likelihood import inject_bias
from allomix.genotype import InformativeMarker, MarkerData


def _chrom_sort_key(chrom: str) -> tuple[int, int]:
    """Sort key for chromosome names (chr1, chr2, ... chr22, chrX, chrY)."""
    name = chrom.replace("chr", "")
    try:
        return (0, int(name))
    except ValueError:
        return (1, ord(name[0]) if name else 0)


@dataclass
class VcfRecord:
    """A single VCF data line, lightly parsed."""

    chrom: str
    pos: int
    id_: str
    ref: str
    alt: str
    qual: str
    filter_: str
    info: str
    format_: str
    sample: str

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.pos}"


def parse_text_vcf(path: str | Path) -> tuple[list[str], list[VcfRecord]]:
    """Read a plain-text VCF and return (header_lines, records).

    Kept separate from ``allomix.genotype.parse_vcf`` (cyvcf2-backed) so the
    simulator stays dependency-light and round-trips raw VCF text.
    """
    header: list[str] = []
    records: list[VcfRecord] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                continue
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            records.append(
                VcfRecord(
                    chrom=fields[0],
                    pos=int(fields[1]),
                    id_=fields[2],
                    ref=fields[3],
                    alt=fields[4],
                    qual=fields[5],
                    filter_=fields[6],
                    info=fields[7],
                    format_=fields[8],
                    sample=fields[9],
                )
            )
    return header, records


def extract_gt(record: VcfRecord) -> tuple[int, int] | None:
    """Extract the diploid genotype allele indices, or None if missing/nocall."""
    fmt_keys = record.format_.split(":")
    fmt_vals = record.sample.split(":")
    gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else None
    if gt_idx is None:
        return None
    gt_str = fmt_vals[gt_idx]
    if gt_str in ("./.", ".|.", "."):
        return None
    sep = "|" if "|" in gt_str else "/"
    parts = gt_str.split(sep)
    if len(parts) != 2:
        return None
    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        return None


def extract_depth(record: VcfRecord) -> int | None:
    """Total depth (DP), falling back to summing AD, or None if not determinable."""
    fmt_keys = record.format_.split(":")
    fmt_vals = record.sample.split(":")
    lookup = dict(zip(fmt_keys, fmt_vals))

    if "DP" in lookup and lookup["DP"] != ".":
        try:
            return int(lookup["DP"])
        except ValueError:
            pass

    if "AD" in lookup and lookup["AD"] != ".":
        try:
            return sum(int(x) for x in lookup["AD"].split(","))
        except ValueError:
            pass

    return None


def alt_dose(gt: tuple[int, int]) -> int:
    """Count of ALT alleles (non-zero) in a diploid genotype (0, 1, or 2)."""
    return (1 if gt[0] != 0 else 0) + (1 if gt[1] != 0 else 0)


def expected_vaf(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    donor_fraction: float,
) -> float:
    """Expected ALT VAF in a diploid chimeric mixture.

        expected_vaf = ((1 - f) * host_alt_dose + f * donor_alt_dose) / 2

    where f (``donor_fraction``, 0.0-1.0) is the donor cell fraction.
    """
    h = alt_dose(host_gt)
    d = alt_dose(donor_gt)
    return ((1.0 - donor_fraction) * h + donor_fraction * d) / 2.0


def expected_vaf_multi(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
) -> float:
    """Expected ALT VAF in a multi-donor chimeric mixture.

    VAF = ((1 - f1 - f2 - ...) * host_dose + f1 * d1_dose + f2 * d2_dose + ...) / 2.
    ``donor_fractions`` must sum to <= 1.0.
    """
    f_host = 1.0 - sum(donor_fractions)
    vaf = f_host * alt_dose(host_gt)
    for dgt, f in zip(donor_gts, donor_fractions):
        vaf += f * alt_dose(dgt)
    return vaf / 2.0


@dataclass
class HostAberration:
    """A somatic copy-number aberration carried by the host (recipient) clone.

    At one marker the host is a mixture: a fraction ``clonal_fraction`` (0.0-1.0)
    of cells carry the aberration (copy number ``cn`` with ``alt_copies`` ALT
    alleles), the rest are normal diploid germline.

    Copy-neutral LoH (CN-LoH, acquired uniparental disomy) is ``cn=2`` with
    ``alt_copies`` forced to 0 or 2 (clone retains two copies of one germline
    homolog). The same dataclass expresses deletions (``cn=1``) and gains
    (``cn=3``).
    """

    cn: int
    alt_copies: int  # ALT alleles among the cn copies (0..cn)
    clonal_fraction: float


def cn_weighted_vaf(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
    host_aberration: HostAberration | None = None,
) -> float:
    """Expected ALT VAF in a chimeric mixture, weighted by copy number.

    ``expected_vaf_multi`` divides by 2 assuming every genome is diploid. A host
    copy-number aberration changes the clone's allele balance and (for
    non-copy-neutral changes) how much host DNA the locus contributes, so the
    local mixing fraction differs from the genome-wide fractions:

        VAF = sum_i (frac_i * cn_i * alt_frac_i) / sum_i (frac_i * cn_i)

    over contributors i (normal-diploid host, host clone, each donor). With
    ``host_aberration`` None this reduces exactly to ``expected_vaf_multi``.
    """
    f_host = 1.0 - sum(donor_fractions)
    if host_aberration is None:
        return expected_vaf_multi(host_gt, donor_gts, donor_fractions)

    c = host_aberration.clonal_fraction
    # Host splits into a normal diploid sub-population and the aberrant clone.
    num = f_host * (1.0 - c) * alt_dose(host_gt) + f_host * c * host_aberration.alt_copies
    den = f_host * (1.0 - c) * 2.0 + f_host * c * host_aberration.cn
    for dgt, fr in zip(donor_gts, donor_fractions):
        num += fr * alt_dose(dgt)
        den += fr * 2.0
    return num / den if den > 0 else 0.0


CNV_KINDS = ("cnloh", "deletion", "gain")


def _clone_state(
    host_gt: tuple[int, int],
    kind: str,
    clonal_fraction: float,
    rng: random.Random,
) -> HostAberration | None:
    """Derive the clone's copy-number state at one marker from the germline GT.

    Mutates one randomly chosen germline homolog:

      - ``cnloh``: copy-neutral LoH; retained homolog duplicated, clone homozygous
        (cn=2). Invisible at a homozygous germline genotype, returns None there.
      - ``deletion``: one homolog lost (cn=1).
      - ``gain``: one homolog duplicated (cn=3).

    Deletions and gains change the locus DNA mass, so they shift the local mixing
    fraction even at homozygous markers and apply to any germline genotype.
    """
    a = list(host_gt)
    if kind == "cnloh":
        if a[0] == a[1]:
            return None
        retained = a[rng.randint(0, 1)]
        return HostAberration(cn=2, alt_copies=2 * retained, clonal_fraction=clonal_fraction)
    if kind == "deletion":
        kept = a[rng.randint(0, 1)]
        return HostAberration(cn=1, alt_copies=kept, clonal_fraction=clonal_fraction)
    if kind == "gain":
        dup = a[rng.randint(0, 1)]
        return HostAberration(
            cn=3, alt_copies=alt_dose(host_gt) + dup, clonal_fraction=clonal_fraction
        )
    raise ValueError(f"kind must be one of {CNV_KINDS}, got {kind!r}")


def assign_cnv_aberrations(
    markers: list[dict],
    fraction_affected: float,
    clonal_fraction: float,
    rng: random.Random,
    kind: str = "cnloh",
    host_key: str = "host_gt",
) -> list[HostAberration | None]:
    """Assign host copy-number aberrations of one kind to a fraction of markers.

    Each marker is independently eligible with probability ``fraction_affected``.
    For ``cnloh`` only heterozygous germline markers show an effect (see
    ``_clone_state``); for ``deletion`` and ``gain`` every eligible marker is
    affected. ``markers`` must be in the order the blender iterates them.

    Returns a list aligned to ``markers``, each entry a HostAberration or None.
    """
    if kind not in CNV_KINDS:
        raise ValueError(f"kind must be one of {CNV_KINDS}, got {kind!r}")
    out: list[HostAberration | None] = []
    for m in markers:
        if rng.random() < fraction_affected:
            out.append(_clone_state(m[host_key], kind, clonal_fraction, rng))
        else:
            out.append(None)
    return out


def assign_cnloh_aberrations(
    markers: list[dict],
    fraction_affected: float,
    clonal_fraction: float,
    rng: random.Random,
    host_key: str = "host_gt",
) -> list[HostAberration | None]:
    """Assign copy-neutral LoH aberrations to a fraction of host het markers.

    Thin wrapper over ``assign_cnv_aberrations`` with ``kind="cnloh"``. CN-LoH is
    only observable at host heterozygous markers; an affected het marker retains
    one germline homolog at random, so the clone becomes homozygous REF
    (``alt_copies=0``) or homozygous ALT (``alt_copies=2``).
    """
    return assign_cnv_aberrations(
        markers, fraction_affected, clonal_fraction, rng, kind="cnloh", host_key=host_key
    )


def is_informative(host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> bool:
    """Informative when host and donor differ in ALT allele dose.

    Only then does the mixed sample show a VAF shift relative to the host.
    """
    return alt_dose(host_gt) != alt_dose(donor_gt)


def sample_allele_counts(
    vaf: float,
    depth: int,
    rng: random.Random | None = None,
    error_rate: float = 0.0,
    rho: float = float("inf"),
    rho_marker_type: str = "all",
) -> tuple[int, int]:
    """Sample (ref_count, alt_count) from a (beta-)binomial with sequencing errors.

    Matches the error model in ``chimerism.log_likelihood_marker_bb()`` for a
    consistent generative model.

    When ``error_rate`` > 0, reads are mis-called under the 4-state
    (trinucleotide) model: a correct read survives with probability
    ``(1 - error_rate)``, an error is uniform over the 3 other bases, so the
    effective per-read REF->ALT (or ALT->REF) rate at a biallelic site is
    ``error_rate / 3``.

    When ``rho`` is finite, reads are beta-binomial: a per-marker ALT probability
    is drawn from Beta(p*rho, (1-p)*rho), then reads binomial given it. This
    injects the overdispersion of real sequencing (PCR/capture jitter); variance
    is inflated by ``(n+rho)/(rho+1)``, ``rho -> inf`` recovers the pure binomial,
    and ``rho`` is the dominant control on the achievable LoD at high depth.

    ``rho_marker_type`` chooses where overdispersion applies:

      - ``"all"`` (default): every marker is beta-binomial when ``rho`` is finite.
      - ``"het_only"``: overdispersion only at intermediate VAF (``0.05 < vaf <
        0.95``, evaluated on the *unbiased* input before the error transform pulls
        ``p`` off the boundaries); boundary VAF stays binomial. Overdispersion is a
        het/intermediate amplification effect; at a donor-absent allele at the
        sequencing-error background the residual variance is a binomial error
        floor. This is the regime presence-detection at donor-homozygous markers
        targets (see ``claude/20_host_presence_detection_plan.md``).
    """
    if rng is None:
        rng = random.Random()
    if depth <= 0:
        return (0, 0)
    p = max(0.0, min(1.0, vaf))

    # Decide overdispersion from the unbiased VAF, before the error transform
    # shifts p off the boundary.
    if rho_marker_type == "all":
        apply_rho = True
    elif rho_marker_type == "het_only":
        eps = 0.05
        apply_rho = eps < p < 1.0 - eps
    else:
        raise ValueError(f"rho_marker_type must be 'all' or 'het_only', got {rho_marker_type!r}")

    # 4-state error model (matching chimerism.log_likelihood_marker). Estimator
    # allocates p_alt = p*(1-e) + (1-p)*e/3 and p_ref similarly, with
    # p_alt + p_ref = 1 - 2e/3 (rest goes to the 2 non-REF/ALT bases). A binomial
    # simulator classifies every read as REF or ALT, so normalise to the
    # conditional p_alt / (1 - 2e/3): symmetric e/3 error floors and an unbiased
    # MLE under the estimator's likelihood.
    if error_rate > 0:
        e = error_rate
        p_alt = p * (1.0 - e) + (1.0 - p) * e / N_OTHER_BASES
        p = p_alt / (1.0 - 2.0 * e / N_OTHER_BASES)

    # Beta-binomial: draw per-marker ALT probability from Beta(mean p,
    # concentration rho), then sample binomially. At p in {0, 1} the Beta is
    # degenerate, so keep p.
    if apply_rho and math.isfinite(rho) and 0.0 < p < 1.0:
        p = rng.betavariate(p * rho, (1.0 - p) * rho)

    if hasattr(rng, "binomialvariate"):
        alt_count = rng.binomialvariate(depth, p)
    else:
        seed = rng.getrandbits(32)
        alt_count = int(np.random.default_rng(seed).binomial(depth, p))
    ref_count = depth - alt_count
    return (ref_count, alt_count)


def thin_informative_markers(
    markers: list[InformativeMarker],
    rate: float,
    rng: np.random.Generator,
) -> list[InformativeMarker]:
    """Binomially down-sample admix AD by one global keep-rate (samtools -s).

    One ``rate`` (0 < rate <= 1) applied to every marker preserves the
    locus-to-locus depth CV: each marker lands near ``rate`` times its own depth,
    so deep markers stay deeper and shallow ones drop out first. This is the
    analog of ``samtools view -s`` / ``seqtk sample`` on the reads, where each
    read survives independently with probability ``rate``, making the surviving
    ref/alt counts binomial draws.

    Host/donor genotypes, marker type, and bias annotations are preserved; only
    admix counts are resampled. ``rate == 1.0`` is a no-op (cannot upsample).
    Input markers are never mutated (fresh copies returned).

    Raises:
        ValueError: If ``rate`` is not in (0, 1].
    """
    if not 0.0 < rate <= 1.0:
        raise ValueError(f"rate must be in (0, 1], got {rate!r}")
    if rate == 1.0:
        return list(markers)
    out: list[InformativeMarker] = []
    for m in markers:
        new_ref = int(rng.binomial(m.admix_ad_ref, rate))
        new_alt = int(rng.binomial(m.admix_ad_alt, rate))
        out.append(
            replace(
                m,
                admix_ad_ref=new_ref,
                admix_ad_alt=new_alt,
                admix_dp=new_ref + new_alt,
            )
        )
    return out


def gt_from_counts(ref_count: int, alt_count: int) -> str:
    """Call a genotype string ('./.', '0/0', '0/1', '1/1') from allele counts."""
    total = ref_count + alt_count
    if total == 0:
        return "./."
    af = alt_count / total
    if af < HOM_REF_MAX_VAF:
        return "0/0"
    if af > HOM_ALT_MIN_VAF:
        return "1/1"
    return "0/1"


@dataclass
class BlendResult:
    """Result of blending two VCFs at a given mixture fraction."""

    header: list[str]
    records: list[str]
    num_markers: int
    num_informative: int
    marker_biases: list[tuple[str, int, str, str, float]] | None = None
    # (chrom, pos, ref, alt, bias) per shared marker, or None if no bias
    markers: list[MarkerData] | None = None
    # Populated only when blend_vcfs(return_markers=True). Identical to
    # parse_vcf(write_vcf(result)) for the SNP panels the simulator emits, but
    # built in memory to skip the disk write/parse round-trip.


def sample_marker_depths(
    n_markers: int,
    mean_depth: int,
    depth_cv: float,
    rng: random.Random,
) -> list[int]:
    """Draw per-marker depths from a log-normal with E[X]=mean_depth, CV[X]=depth_cv.

    Real panels vary depth substantially across markers (empirically CV=0.43 on
    the 76-SNP rhAmpSeq panel). ``depth_cv=0`` gives uniform depth. Depths are
    integers, minimum 1.
    """
    if depth_cv <= 0:
        return [mean_depth] * n_markers
    sigma2 = math.log(1 + depth_cv**2)
    mu = math.log(mean_depth) - sigma2 / 2
    sigma = math.sqrt(sigma2)
    return [max(1, round(math.exp(rng.gauss(mu, sigma)))) for _ in range(n_markers)]


def generate_marker_biases(
    n_markers: int,
    rng: random.Random,
    bias_sd: float = 0.02,
) -> list[float]:
    """Generate per-marker capture/amplification biases, one fixed N(0, bias_sd) each.

    Models the systematic REF/ALT capture efficiency difference of real
    hybridisation-capture and amplicon data (Vynck et al.). Bias is in het-site
    VAF units: +0.02 shifts a true het's observed ALT VAF to 0.52. Injected
    multiplicatively in logit space (``allomix.estimate.likelihood.inject_bias``), so away
    from VAF 0.5 the shift is proportional, not additive, matching how the
    estimator corrects for it (issue #20).

    ``bias_sd`` typical values: 0.0 = none, 0.02 = realistic (empirically 0.019 on
    the 76-SNP rhAmpSeq panel across 210 joint-called VCFs), 0.05 = poor panel.
    """
    if bias_sd <= 0:
        return [0.0] * n_markers
    return [rng.gauss(0.0, bias_sd) for _ in range(n_markers)]


def generate_marker_biases_realistic(
    n_markers: int,
    rng: random.Random,
    sd: float = 0.012,
    outlier_frac: float = 0.05,
    outlier_sd: float = 0.08,
) -> list[float]:
    """Generate biases with a heavy-tailed mixture distribution.

    The empirical bias distribution is heavy-tailed (median |bias| 0.005, 95th
    pct 0.041, max 0.10), which a simple Gaussian underestimates. Mixture:
    ``1 - outlier_frac`` of markers from N(0, sd), the rest from N(0, outlier_sd).
    Defaults calibrated from 71 markers across 210 joint-called VCFs on the 76-SNP
    rhAmpSeq panel, giving overall SD ~0.018 matching the empirical measurement.
    """
    biases = []
    for _ in range(n_markers):
        if rng.random() < outlier_frac:
            biases.append(rng.gauss(0, outlier_sd))
        else:
            biases.append(rng.gauss(0, sd))
    return biases


# IBD sharing probabilities: (P(IBD=0), P(IBD=1), P(IBD=2))
RELATEDNESS_IBD = {
    "unrelated": (1.0, 0.0, 0.0),
    "cousin": (0.75, 0.25, 0.0),  # first cousins: 1/8 kinship
    "half-sibling": (0.5, 0.5, 0.0),  # half-siblings: 1/4 kinship
    "parent-child": (0.0, 1.0, 0.0),  # parent-child: always share 1 allele
    "sibling": (0.25, 0.5, 0.25),  # full siblings: 1/4 kinship
}


def _draw_genotype(p_alt: float, rng: random.Random) -> tuple[int, int]:
    """Draw a diploid genotype under Hardy-Weinberg at ALT frequency ``p_alt``."""
    a1 = 1 if rng.random() < p_alt else 0
    a2 = 1 if rng.random() < p_alt else 0
    return (a1, a2)


def _draw_related_genotype(
    host_gt: tuple[int, int],
    p_alt: float,
    ibd_probs: tuple[float, float, float],
    rng: random.Random,
) -> tuple[int, int]:
    """Draw a donor genotype conditional on host genotype and IBD sharing.

    ``ibd_probs`` is (P(IBD=0), P(IBD=1), P(IBD=2)).
    """
    r = rng.random()
    if r < ibd_probs[0]:
        # IBD=0: independent draw
        return _draw_genotype(p_alt, rng)
    if r < ibd_probs[0] + ibd_probs[1]:
        # IBD=1: share one allele, draw the other independently
        shared = host_gt[rng.randint(0, 1)]
        other = 1 if rng.random() < p_alt else 0
        return (shared, other) if rng.random() < 0.5 else (other, shared)
    # IBD=2: identical genotype
    return host_gt


def generate_related_genotypes(
    n_markers: int,
    relatedness: str,
    rng: random.Random,
    maf_range: tuple[float, float] = (0.2, 0.5),
) -> list[dict]:
    """Generate synthetic host-donor genotype pairs with specified relatedness.

    Marker allele frequencies are drawn uniformly from ``maf_range``, then host
    and donor genotypes are generated with the matching IBD sharing.
    ``relatedness`` is a key of ``RELATEDNESS_IBD``.

    Returns a list of dicts keyed: chrom, pos, ref, alt, host_gt, donor_gt,
    p_alt, informative.
    """
    if relatedness not in RELATEDNESS_IBD:
        raise ValueError(
            f"Unknown relatedness '{relatedness}'. Choose from: {list(RELATEDNESS_IBD.keys())}"
        )
    ibd_probs = RELATEDNESS_IBD[relatedness]

    markers = []
    for i in range(n_markers):
        p_alt = rng.uniform(*maf_range)
        host_gt = _draw_genotype(p_alt, rng)
        donor_gt = _draw_related_genotype(host_gt, p_alt, ibd_probs, rng)

        markers.append(
            {
                "chrom": "chr1",
                "pos": 10000 + i * 1000,
                "ref": "A",
                "alt": "G",
                "host_gt": host_gt,
                "donor_gt": donor_gt,
                "p_alt": p_alt,
                "informative": alt_dose(host_gt) != alt_dose(donor_gt),
            }
        )

    return markers


def generate_paired_related_genotypes(
    n_markers: int,
    relatedness_levels: list[str],
    rng: random.Random,
    maf_range: tuple[float, float] = (0.2, 0.5),
) -> dict[str, list[dict]]:
    """Generate host-donor panels for several relatedness levels sharing one host.

    One host panel (allele frequencies and host genotypes) is generated, and each
    level's donor is derived from it using the same per-marker random draws mapped
    through the level's IBD probabilities. Only the IBD sharing varies across
    levels, so the informative-marker count is monotone non-increasing as
    relatedness rises (more sharing leaves fewer differing markers). This removes
    the pair-to-pair noise that, with independent draws, can leave a more related
    level (e.g. first cousin) showing more informative markers than a less related
    one (e.g. unrelated) by chance.

    Each level in ``relatedness_levels`` is a key of ``RELATEDNESS_IBD``. Returns
    a dict mapping each level to its marker-dict list (same keys as
    :func:`generate_related_genotypes`); host genotype and ``p_alt`` at marker i
    are identical across levels.
    """
    unknown = [r for r in relatedness_levels if r not in RELATEDNESS_IBD]
    if unknown:
        raise ValueError(
            f"Unknown relatedness {unknown}. Choose from: {list(RELATEDNESS_IBD.keys())}"
        )

    panels: dict[str, list[dict]] = {rel: [] for rel in relatedness_levels}
    for i in range(n_markers):
        p_alt = rng.uniform(*maf_range)
        host_gt = _draw_genotype(p_alt, rng)

        # Pre-draw per-marker randomness once, shared across all levels.
        r_ibd = rng.random()
        share_idx = rng.randint(0, 1)
        other_allele = 1 if rng.random() < p_alt else 0
        swap = rng.random() < 0.5
        indep_gt = _draw_genotype(p_alt, rng)

        chrom = "chr1"
        pos = 10000 + i * 1000
        for rel in relatedness_levels:
            ibd_probs = RELATEDNESS_IBD[rel]
            if r_ibd < ibd_probs[0]:
                donor_gt = indep_gt  # IBD=0: independent draw
            elif r_ibd < ibd_probs[0] + ibd_probs[1]:
                shared = host_gt[share_idx]  # IBD=1: share one allele
                donor_gt = (shared, other_allele) if swap else (other_allele, shared)
            else:
                donor_gt = host_gt  # IBD=2: identical genotype

            panels[rel].append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": "A",
                    "alt": "G",
                    "host_gt": host_gt,
                    "donor_gt": donor_gt,
                    "p_alt": p_alt,
                    "informative": alt_dose(host_gt) != alt_dose(donor_gt),
                }
            )

    return panels


def _mendelian_child(
    parent1: tuple[int, int],
    parent2: tuple[int, int],
    rng: random.Random,
) -> tuple[int, int]:
    """Draw a child genotype by Mendelian segregation (one allele per parent).

    Returned sorted, smaller allele first.
    """
    a1 = parent1[rng.randint(0, 1)]
    a2 = parent2[rng.randint(0, 1)]
    return (min(a1, a2), max(a1, a2))


def generate_sibling_trio_genotypes(
    n_markers: int,
    rng: random.Random,
    maf_range: tuple[float, float] = (0.2, 0.5),
) -> list[dict]:
    """Generate genotypes for 3 siblings (host + 2 donors) from shared parents.

    Per marker: draw an ALT frequency, draw two parent genotypes from
    Hardy-Weinberg, then derive each sibling independently by Mendelian
    segregation. This preserves the 3-way sibling correlation: each pair has IBD
    distribution (0.25, 0.5, 0.25) and all three are correlated through the shared
    parents.

    Returns a list of dicts keyed: chrom, pos, ref, alt, host_gt, donor1_gt,
    donor2_gt, p_alt, informative_d1, informative_d2, informative_any,
    donors_distinguishable.
    """
    markers = []
    for i in range(n_markers):
        p_alt = rng.uniform(*maf_range)

        parent1 = _draw_genotype(p_alt, rng)
        parent2 = _draw_genotype(p_alt, rng)

        host_gt = _mendelian_child(parent1, parent2, rng)
        donor1_gt = _mendelian_child(parent1, parent2, rng)
        donor2_gt = _mendelian_child(parent1, parent2, rng)

        markers.append(
            {
                "chrom": f"chr{(i % 22) + 1}",
                "pos": 1_000_000 + i * 100_000,
                "ref": "A",
                "alt": "G",
                "host_gt": host_gt,
                "donor1_gt": donor1_gt,
                "donor2_gt": donor2_gt,
                "p_alt": p_alt,
                "informative_d1": alt_dose(host_gt) != alt_dose(donor1_gt),
                "informative_d2": alt_dose(host_gt) != alt_dose(donor2_gt),
                "informative_any": (
                    alt_dose(host_gt) != alt_dose(donor1_gt)
                    or alt_dose(host_gt) != alt_dose(donor2_gt)
                ),
                "donors_distinguishable": alt_dose(donor1_gt) != alt_dose(donor2_gt),
            }
        )

    return markers


def write_genotype_vcf(
    markers: list[dict],
    path: str | Path,
    sample_name: str,
    key: str = "host_gt",
    depth: int = 100,
) -> None:
    """Write a synthetic genotype VCF from generated marker data.

    ``key`` selects which genotype to write ('host_gt' or 'donor_gt').
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    chroms = sorted(set(m["chrom"] for m in markers), key=_chrom_sort_key)

    info_format_header = """\
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
"""
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name]

    with open(path, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        for chrom in chroms:
            f.write(f"##contig=<ID={chrom}>\n")
        f.write(info_format_header)
        f.write("\t".join(columns) + "\n")

        for m in markers:
            gt = m[key]
            gt_str = f"{gt[0]}/{gt[1]}"
            n_alt = alt_dose(gt)
            ad_alt = round(depth * n_alt / 2)
            ad_ref = depth - ad_alt
            sample_field = f"{gt_str}:{ad_ref},{ad_alt}:{depth}:99"
            row = [
                m["chrom"],
                str(m["pos"]),
                ".",
                m["ref"],
                m["alt"],
                ".",
                "PASS",
                ".",
                "GT:AD:DP:GQ",
                sample_field,
            ]
            f.write("\t".join(row) + "\n")


def blend_vcfs(
    host_path: str | Path,
    donor_path: str | Path,
    donor_fraction: float,
    target_depth: int | None = None,
    sample_name: str | None = None,
    seed: int | None = None,
    marker_bias_sd: float = 0.0,
    fixed_biases: list[float] | None = None,
    error_rate: float = DEFAULT_ERROR_RATE,
    allele_dropout_rate: float = 0.0,
    locus_dropout_rate: float = 0.0,
    depth_cv: float = 0.0,
    realistic_biases: bool = False,
    rho: float = float("inf"),
    rho_marker_type: str = "all",
    host_aberrations: list[HostAberration | None] | None = None,
    return_markers: bool = False,
) -> BlendResult:
    """Blend two genotype VCFs to create a synthetic chimeric VCF.

    Args:
        donor_fraction: Fraction of donor DNA (0.0 to 1.0).
        target_depth: Fixed depth for all markers; None uses the host VCF depth.
        marker_bias_sd: Per-marker capture bias SD (0.0 ideal, 0.02 realistic).
            Ignored when ``fixed_biases`` is given or ``realistic_biases`` is set.
        fixed_biases: Pre-generated per-marker biases used directly. Length must
            match the number of shared markers.
        allele_dropout_rate: Per-marker probability (0.0-1.0) that one allele's
            reads are entirely lost at a het site, making a het look like a hom.
        locus_dropout_rate: Per-marker probability (0.0-1.0) of zero reads.
        depth_cv: Per-marker depth CV (0.0 uniform, 0.43 empirical rhAmpSeq).
            Applied only when > 0 and ``target_depth`` is set (log-normal draw).
        realistic_biases: Use the heavy-tailed mixture bias distribution instead
            of a Gaussian. Ignored when ``fixed_biases`` is given.
        rho: Beta-binomial overdispersion (inf = pure binomial); smaller raises
            the achievable LoD. See ``sample_allele_counts``.
        rho_marker_type: ``"all"`` (default) or ``"het_only"``.
        host_aberrations: Optional per-shared-marker host copy-number aberrations
            (see ``HostAberration``). Affected markers use the copy-number-weighted
            mixture instead of the diploid model. Aligns with the shared markers in
            iteration order like ``fixed_biases``; ``None`` entries stay diploid.
    """
    if not 0.0 <= donor_fraction <= 1.0:
        raise ValueError(f"donor_fraction must be 0.0-1.0, got {donor_fraction}")

    rng = random.Random(seed)
    host_header, host_records = parse_text_vcf(host_path)
    _, donor_records = parse_text_vcf(donor_path)

    donor_by_locus: dict[str, VcfRecord] = {}
    for rec in donor_records:
        donor_by_locus[rec.locus] = rec

    if sample_name is None:
        sample_name = "simulated"

    # Replace the sample name in the #CHROM line, copy other header lines as-is.
    out_header = []
    for line in host_header:
        if line.startswith("#CHROM"):
            parts = line.split("\t")
            parts[-1] = sample_name
            out_header.append("\t".join(parts))
        else:
            out_header.append(line)

    # Per-marker biases, one per shared locus in iteration order.
    n_shared = sum(1 for hr in host_records if hr.locus in donor_by_locus)
    if fixed_biases is not None:
        if len(fixed_biases) != n_shared:
            raise ValueError(
                f"fixed_biases length ({len(fixed_biases)}) != shared markers ({n_shared})"
            )
        marker_biases = fixed_biases
    elif realistic_biases:
        marker_biases = generate_marker_biases_realistic(n_shared, rng)
    else:
        marker_biases = generate_marker_biases(n_shared, rng, marker_bias_sd)

    if host_aberrations is not None and len(host_aberrations) != n_shared:
        raise ValueError(
            f"host_aberrations length ({len(host_aberrations)}) != shared markers ({n_shared})"
        )

    if depth_cv > 0 and target_depth is not None:
        marker_depths = sample_marker_depths(n_shared, target_depth, depth_cv, rng)
    else:
        marker_depths = None  # use flat target_depth or host depth

    out_records: list[str] = []
    bias_info: list[tuple[str, int, str, str, float]] = []
    num_markers = 0
    num_informative = 0
    bias_idx = 0
    out_markers: list[MarkerData] = []

    for host_rec in host_records:
        donor_rec = donor_by_locus.get(host_rec.locus)
        if donor_rec is None:
            continue

        host_gt = extract_gt(host_rec)
        donor_gt = extract_gt(donor_rec)
        if host_gt is None or donor_gt is None:
            continue

        # Must share the same REF allele.
        if host_rec.ref != donor_rec.ref:
            continue

        if marker_depths is not None:
            depth = marker_depths[bias_idx]
        elif target_depth is not None:
            depth = target_depth
        else:
            depth = extract_depth(host_rec) or 1000

        # Locus dropout: marker produces zero reads.
        if locus_dropout_rate > 0 and rng.random() < locus_dropout_rate:
            bias_idx += 1
            continue

        num_markers += 1
        if is_informative(host_gt, donor_gt):
            num_informative += 1

        # A host copy-number aberration replaces the diploid model at this marker.
        aberr = host_aberrations[bias_idx] if host_aberrations is not None else None
        if aberr is not None:
            vaf = cn_weighted_vaf(host_gt, [donor_gt], [donor_fraction], aberr)
        else:
            vaf = expected_vaf(host_gt, donor_gt, donor_fraction)
        this_bias = marker_biases[bias_idx]
        vaf_biased = float(inject_bias(vaf, this_bias)) if this_bias != 0.0 else vaf
        alt_allele_bias = host_rec.alt if host_rec.alt != "." else donor_rec.alt
        bias_info.append((host_rec.chrom, host_rec.pos, host_rec.ref, alt_allele_bias, this_bias))
        bias_idx += 1

        # Allele dropout at a het-like site: one allele lost, VAF pushed to 0 or 1.
        if allele_dropout_rate > 0 and 0.05 < vaf_biased < 0.95:
            if rng.random() < allele_dropout_rate:
                vaf_biased = 0.0 if rng.random() < 0.5 else 1.0

        ref_count, alt_count = sample_allele_counts(
            vaf_biased, depth, rng, error_rate, rho, rho_marker_type
        )

        gt = gt_from_counts(ref_count, alt_count)
        total = ref_count + alt_count
        af_val = f"{alt_count / total:.4f}" if total > 0 else "0"
        gq = 99
        # Simplified PL: 0 for the called genotype.
        if gt == "0/0":
            pl = f"0,{gq},{gq * 10}"
        elif gt == "1/1":
            pl = f"{gq * 10},{gq},0"
        else:
            pl = f"{gq * 5},0,{gq * 5}"

        sample_field = f"{gt}:{ref_count},{alt_count}:{total}:{gq}:{pl}:{af_val}"
        format_field = "GT:AD:DP:GQ:PL:AF"

        # Use the host's ALT allele; if host was ref-only, use the donor's.
        alt_allele = host_rec.alt if host_rec.alt != "." else donor_rec.alt
        if alt_allele == ".":
            # Both hom-ref, no ALT listed: still emit the site, but no AF field.
            sample_field = f"{gt}:{ref_count}:{total}:{gq}:{pl}"
            format_field = "GT:AD:DP:GQ:PL"

        info_parts = [f"DP={total}"]
        if alt_allele != ".":
            ac = alt_count
            an = 2
            info_parts.extend([f"AC={ac}", f"AN={an}"])

        line = "\t".join(
            [
                host_rec.chrom,
                str(host_rec.pos),
                host_rec.id_,
                host_rec.ref,
                alt_allele,
                str(host_rec.qual),
                "PASS",
                ";".join(info_parts),
                format_field,
                sample_field,
            ]
        )
        out_records.append(line)

        # In-memory MarkerData matching parse_vcf(write_vcf(...)). Mirror
        # parse_vcf's skips: drop multiallelic (never emitted here) and indels;
        # keep alt="." sites (ad_alt reads 0, as the round-trip would).
        if return_markers and (
            alt_allele == "." or (len(host_rec.ref) == 1 and len(alt_allele) == 1)
        ):
            a1, a2 = (int(x) for x in gt.split("/"))
            out_markers.append(
                MarkerData(
                    chrom=host_rec.chrom,
                    pos=host_rec.pos,
                    ref=host_rec.ref,
                    alt=alt_allele,
                    gt=(min(a1, a2), max(a1, a2)),
                    ad_ref=ref_count,
                    ad_alt=alt_count if alt_allele != "." else 0,
                    dp=total,
                    gq=gq,
                    filter="PASS",
                )
            )

    return BlendResult(
        header=out_header,
        records=out_records,
        num_markers=num_markers,
        num_informative=num_informative,
        marker_biases=bias_info
        if (marker_bias_sd > 0 or fixed_biases is not None or realistic_biases)
        else None,
        markers=out_markers if return_markers else None,
    )


def write_vcf(result: BlendResult, path: str | Path) -> None:
    """Write a BlendResult to a VCF file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for line in result.header:
            fh.write(line + "\n")
        for line in result.records:
            fh.write(line + "\n")


def blend_from_genotype_dicts(
    markers: list[dict],
    donor_fractions: list[float],
    target_depth: int = 1000,
    seed: int | None = None,
    error_rate: float = DEFAULT_ERROR_RATE,
    depth_cv: float = 0.0,
    sample_name: str = "simulated",
    rho: float = float("inf"),
    rho_marker_type: str = "all",
) -> BlendResult:
    """Create a synthetic chimeric VCF directly from genotype dicts.

    Designed for ``generate_sibling_trio_genotypes()`` output. Supports 1 or 2
    donors via the length of ``donor_fractions`` ([f1] or [f1, f2]).
    """
    if sum(donor_fractions) > 1.0 + 1e-9:
        raise ValueError(f"donor_fractions sum to {sum(donor_fractions):.4f}, must be <= 1.0")

    rng = random.Random(seed)
    n = len(markers)

    if depth_cv > 0:
        depths = sample_marker_depths(n, target_depth, depth_cv, rng)
    else:
        depths = [target_depth] * n

    chroms = sorted(set(m["chrom"] for m in markers), key=_chrom_sort_key)
    header = [
        "##fileformat=VCFv4.2",
        *[f"##contig=<ID={c}>" for c in chroms],
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
    ]

    out_records = []
    n_informative = 0

    n_donors = len(donor_fractions)
    donor_keys = [f"donor{i + 1}_gt" for i in range(n_donors)]

    for i, m in enumerate(markers):
        host_gt = m["host_gt"]
        donor_gts = [m[k] for k in donor_keys]

        vaf = expected_vaf_multi(host_gt, donor_gts, donor_fractions)
        ref_count, alt_count = sample_allele_counts(
            vaf, depths[i], rng, error_rate, rho, rho_marker_type
        )

        if m.get("informative_any", False):
            n_informative += 1

        total = ref_count + alt_count
        gt = gt_from_counts(ref_count, alt_count)
        af_val = f"{alt_count / total:.4f}" if total > 0 else "0"
        sample_field = f"{gt}:{ref_count},{alt_count}:{total}:99:{af_val}"

        line = (
            f"{m['chrom']}\t{m['pos']}\t.\t{m['ref']}\t{m['alt']}\t"
            f".\tPASS\tDP={total}\tGT:AD:DP:GQ:AF\t{sample_field}"
        )
        out_records.append(line)

    return BlendResult(
        header=header,
        records=out_records,
        num_markers=n,
        num_informative=n_informative,
    )


@dataclass
class JointVcfResult:
    """Result of building a multi-sample joint VCF."""

    header: list[str]
    records: list[str]
    num_markers: int
    num_informative: int
    sample_names: list[str]


def _simulate_genotype_sample(
    gt: tuple[int, int],
    depth: int,
    rng: random.Random,
) -> str:
    """Simulate a FORMAT sample field for a genotyping sample.

    Draws allele counts binomially from the true genotype, as sequencing of a
    pure (non-admixed) sample would produce.
    """
    vaf = alt_dose(gt) / 2.0
    ref_count, alt_count = sample_allele_counts(vaf, depth, rng, error_rate=0.01)
    total = ref_count + alt_count
    gt_str = f"{gt[0]}/{gt[1]}"
    af_val = f"{alt_count / total:.4f}" if total > 0 else "0"
    return f"{gt_str}:{ref_count},{alt_count}:{total}:99:{af_val}"


def build_joint_vcf(
    host_path: str | Path,
    donor_paths: list[str | Path],
    admix_fractions: list[float],
    admix_sample_names: list[str],
    host_sample_name: str = "HOST",
    donor_sample_names: list[str] | None = None,
    target_depth: int | None = None,
    seed: int | None = None,
    error_rate: float = DEFAULT_ERROR_RATE,
    depth_cv: float = 0.0,
    marker_bias_sd: float = 0.0,
) -> JointVcfResult:
    """Build a multi-sample joint VCF of host, donor(s), and admixture samples.

    Simulates GATK joint calling: all samples in one VCF, with ALT alleles
    discovered anywhere propagated to every sample's AD field.

    Args:
        admix_fractions: Donor fraction per admixture sample. Single-donor only
            (one float each); multi-donor per-sample lists are not yet supported.
        donor_sample_names: Defaults to DONOR, or DONOR1/DONOR2/... for multiple.
        target_depth: Fixed depth for all markers/samples; None uses host VCF depth.
    """
    if len(admix_fractions) != len(admix_sample_names):
        raise ValueError(
            f"admix_fractions length ({len(admix_fractions)}) != "
            f"admix_sample_names length ({len(admix_sample_names)})"
        )

    rng = random.Random(seed)

    host_header, host_records = parse_text_vcf(host_path)
    donor_record_lists = [parse_text_vcf(dp)[1] for dp in donor_paths]

    if donor_sample_names is None:
        if len(donor_paths) == 1:
            donor_sample_names = ["DONOR"]
        else:
            donor_sample_names = [f"DONOR{i + 1}" for i in range(len(donor_paths))]

    all_sample_names = [host_sample_name] + donor_sample_names + admix_sample_names

    donor_by_locus = [{rec.locus: rec for rec in donor_recs} for donor_recs in donor_record_lists]

    # Loci present in host and all donors.
    shared_loci = []
    for host_rec in host_records:
        if all(host_rec.locus in dbl for dbl in donor_by_locus):
            shared_loci.append(host_rec)

    n_shared = len(shared_loci)

    marker_biases = generate_marker_biases(n_shared, rng, marker_bias_sd)
    if depth_cv > 0 and target_depth is not None:
        marker_depths = sample_marker_depths(n_shared, target_depth, depth_cv, rng)
    else:
        marker_depths = None

    out_header = []
    for line in host_header:
        if line.startswith("#CHROM"):
            parts = line.split("\t")[:9]
            parts.extend(all_sample_names)
            out_header.append("\t".join(parts))
        else:
            out_header.append(line)

    out_records = []
    num_informative = 0
    format_field = "GT:AD:DP:GQ:AF"

    for idx, host_rec in enumerate(shared_loci):
        host_gt = extract_gt(host_rec)
        donor_recs = [dbl[host_rec.locus] for dbl in donor_by_locus]
        donor_gts = [extract_gt(dr) for dr in donor_recs]

        if host_gt is None or any(dg is None for dg in donor_gts):
            continue

        if any(dr.ref != host_rec.ref for dr in donor_recs):
            continue

        # ALT allele: first non-"." ALT from any sample.
        alt_allele = host_rec.alt
        if alt_allele == ".":
            for dr in donor_recs:
                if dr.alt != ".":
                    alt_allele = dr.alt
                    break
        if alt_allele == ".":
            alt_allele = "."

        if marker_depths is not None:
            depth = marker_depths[idx]
        elif target_depth is not None:
            depth = target_depth
        else:
            depth = extract_depth(host_rec) or 1000

        bias = marker_biases[idx]

        if any(is_informative(host_gt, dg) for dg in donor_gts):
            num_informative += 1

        sample_fields = []
        sample_fields.append(_simulate_genotype_sample(host_gt, depth, rng))
        for dg in donor_gts:
            sample_fields.append(_simulate_genotype_sample(dg, depth, rng))

        for frac in admix_fractions:
            if len(donor_gts) == 1:
                vaf = expected_vaf(host_gt, donor_gts[0], frac)
            else:
                # Multi-donor with a single fraction: split it equally.
                per_donor = [frac / len(donor_gts)] * len(donor_gts)
                vaf = expected_vaf_multi(host_gt, donor_gts, per_donor)

            vaf_biased = float(inject_bias(vaf, bias)) if bias != 0.0 else vaf
            ref_count, alt_count = sample_allele_counts(vaf_biased, depth, rng, error_rate)
            total = ref_count + alt_count
            gt = gt_from_counts(ref_count, alt_count)
            af_val = f"{alt_count / total:.4f}" if total > 0 else "0"
            sample_fields.append(f"{gt}:{ref_count},{alt_count}:{total}:99:{af_val}")

        info_parts = [f"DP={depth}"]

        line = "\t".join(
            [
                host_rec.chrom,
                str(host_rec.pos),
                host_rec.id_,
                host_rec.ref,
                alt_allele,
                str(host_rec.qual),
                "PASS",
                ";".join(info_parts),
                format_field,
                *sample_fields,
            ]
        )
        out_records.append(line)

    return JointVcfResult(
        header=out_header,
        records=out_records,
        num_markers=len(out_records),
        num_informative=num_informative,
        sample_names=all_sample_names,
    )


def write_joint_vcf(result: JointVcfResult, path: str | Path) -> None:
    """Write a JointVcfResult to a VCF file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for line in result.header:
            fh.write(line + "\n")
        for line in result.records:
            fh.write(line + "\n")
