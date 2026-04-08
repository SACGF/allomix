"""Synthetic chimeric VCF generator for testing allomix.

Blends two genotype VCFs at a specified mixture fraction to produce a synthetic
chimeric VCF with realistic allele counts drawn from a binomial distribution.

Uses plain-text VCF parsing only (no cyvcf2 dependency) so this module can be
used in test environments without compiled libraries.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------


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
        """Return a unique key for the genomic position."""
        return f"{self.chrom}:{self.pos}"


# ---------------------------------------------------------------------------
# VCF I/O helpers
# ---------------------------------------------------------------------------


def parse_vcf(path: str | Path) -> tuple[list[str], list[VcfRecord]]:
    """Read a VCF file and return (header_lines, records).

    Args:
        path: Path to a plain-text VCF file (not gzipped).

    Returns:
        A tuple of (header_lines, records) where header_lines includes all
        lines starting with '#' and records is a list of VcfRecord.
    """
    header: list[str] = []
    records: list[VcfRecord] = []
    with open(path) as fh:
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
    """Extract the diploid genotype allele indices from a VcfRecord.

    Args:
        record: A VcfRecord with GT as the first FORMAT field.

    Returns:
        Tuple of (allele1, allele2) as ints, or None if GT is missing/nocall.
    """
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
    """Extract the total depth (DP) from a VcfRecord.

    Falls back to summing AD if DP is not available.

    Args:
        record: A VcfRecord.

    Returns:
        Total read depth as int, or None if not determinable.
    """
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
    """Count the number of ALT alleles (non-zero) in a diploid genotype.

    Args:
        gt: Tuple of two allele indices.

    Returns:
        Number of ALT alleles (0, 1, or 2).
    """
    return (1 if gt[0] != 0 else 0) + (1 if gt[1] != 0 else 0)


# ---------------------------------------------------------------------------
# Core simulation logic
# ---------------------------------------------------------------------------


def expected_vaf(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    donor_fraction: float,
) -> float:
    """Calculate the expected ALT VAF in a chimeric mixture.

    The mixture model assumes diploid genomes:
        expected_vaf = ((1 - f) * host_alt_dose + f * donor_alt_dose) / 2

    where f is the fraction of cells from the donor.

    Args:
        host_gt: Host diploid genotype as (allele1, allele2).
        donor_gt: Donor diploid genotype as (allele1, allele2).
        donor_fraction: Fraction of donor DNA in mixture (0.0 to 1.0).

    Returns:
        Expected ALT allele frequency (0.0 to 1.0).
    """
    h = alt_dose(host_gt)
    d = alt_dose(donor_gt)
    return ((1.0 - donor_fraction) * h + donor_fraction * d) / 2.0


def expected_vaf_multi(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
) -> float:
    """Calculate the expected ALT VAF in a multi-donor chimeric mixture.

    VAF = ((1 - f1 - f2 - ...) * host_dose + f1 * d1_dose + f2 * d2_dose + ...) / 2

    Args:
        host_gt: Host diploid genotype.
        donor_gts: List of donor diploid genotypes.
        donor_fractions: List of donor fractions (must sum to <= 1.0).

    Returns:
        Expected ALT allele frequency (0.0 to 1.0).
    """
    f_host = 1.0 - sum(donor_fractions)
    vaf = f_host * alt_dose(host_gt)
    for dgt, f in zip(donor_gts, donor_fractions):
        vaf += f * alt_dose(dgt)
    return vaf / 2.0


def is_informative(host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> bool:
    """Determine whether a marker is informative for chimerism detection.

    A marker is informative when host and donor have different ALT allele doses,
    which means the mixed sample will show a VAF shift relative to the host.

    Args:
        host_gt: Host diploid genotype.
        donor_gt: Donor diploid genotype.

    Returns:
        True if the marker is informative.
    """
    return alt_dose(host_gt) != alt_dose(donor_gt)


def sample_allele_counts(
    vaf: float,
    depth: int,
    rng: random.Random | None = None,
    error_rate: float = 0.0,
) -> tuple[int, int]:
    """Sample allele counts from a binomial distribution with sequencing errors.

    When ``error_rate`` > 0, each simulated read has a chance of being
    mis-called: a true REF read becomes ALT (or vice-versa) with
    probability ``error_rate``.  This mirrors the flat error model used
    in the MLE estimator.

    Args:
        vaf: Expected variant allele frequency (before sequencing error).
        depth: Total read depth to simulate.
        rng: Optional Random instance for reproducibility.
        error_rate: Per-read sequencing error probability (0.0–1.0).
            Each read of the wrong allele is flipped with this probability.
            Default 0.0 (no sequencing error).

    Returns:
        Tuple of (ref_count, alt_count).
    """
    if rng is None:
        rng = random.Random()
    if depth <= 0:
        return (0, 0)
    # Clamp VAF to valid probability range
    p = max(0.0, min(1.0, vaf))

    # Apply sequencing error: a read from the true allele is mis-called
    # with probability error_rate.  The effective observed ALT probability
    # becomes: p_obs = p*(1-e) + (1-p)*e  (symmetric error model).
    if error_rate > 0:
        p = p * (1.0 - error_rate) + (1.0 - p) * error_rate

    if hasattr(rng, "binomialvariate"):
        alt_count = rng.binomialvariate(depth, p)
    else:
        alt_count = _binomial(rng, depth, p)
    ref_count = depth - alt_count
    return (ref_count, alt_count)


def _binomial(rng: random.Random, n: int, p: float) -> int:
    """Fallback binomial sampling using inverse-transform for small n.

    For large n, uses the normal approximation. This avoids depending on numpy
    just for test infrastructure.

    Args:
        rng: Random instance.
        n: Number of trials.
        p: Probability of success.

    Returns:
        Number of successes.
    """
    if p <= 0.0:
        return 0
    if p >= 1.0:
        return n
    if n <= 100:
        # Direct simulation
        count = 0
        for _ in range(n):
            if rng.random() < p:
                count += 1
        return count
    # Normal approximation for large n
    mu = n * p
    sigma = math.sqrt(n * p * (1 - p))
    result = round(rng.gauss(mu, sigma))
    return max(0, min(n, result))


def gt_from_counts(ref_count: int, alt_count: int) -> str:
    """Assign a genotype string from allele counts using simple thresholds.

    Args:
        ref_count: Reference allele read count.
        alt_count: Alternative allele read count.

    Returns:
        Genotype string: '0/0', '0/1', or '1/1'.
    """
    total = ref_count + alt_count
    if total == 0:
        return "./."
    af = alt_count / total
    if af < 0.05:
        return "0/0"
    elif af > 0.95:
        return "1/1"
    else:
        return "0/1"


# ---------------------------------------------------------------------------
# High-level blending
# ---------------------------------------------------------------------------


@dataclass
class BlendResult:
    """Result of blending two VCFs at a given mixture fraction."""

    header: list[str]
    records: list[str]
    num_markers: int
    num_informative: int
    marker_biases: list[tuple[str, int, str, str, float]] | None = None
    # List of (chrom, pos, ref, alt, bias) for each shared marker, or None if no bias


def sample_marker_depths(
    n_markers: int,
    mean_depth: int,
    depth_cv: float,
    rng: random.Random,
) -> list[int]:
    """Draw per-marker depths from a log-normal matching empirical CV.

    In real sequencing panels, depth varies substantially across markers
    (empirically CV=0.43 on 76-SNP rhAmpSeq panel). The log-normal is
    parameterised so that E[X] = mean_depth and CV[X] = depth_cv.

    Args:
        n_markers: Number of markers.
        mean_depth: Target mean depth across markers.
        depth_cv: Coefficient of variation for depth across markers.
            0.0 = uniform depth (all markers get mean_depth).
        rng: Random instance for reproducibility.

    Returns:
        List of per-marker depths (integers, minimum 1).
    """
    if depth_cv <= 0:
        return [mean_depth] * n_markers
    # Log-normal parameters from desired mean and CV
    sigma2 = math.log(1 + depth_cv**2)
    mu = math.log(mean_depth) - sigma2 / 2
    sigma = math.sqrt(sigma2)
    return [max(1, round(math.exp(rng.gauss(mu, sigma)))) for _ in range(n_markers)]


def generate_marker_biases(
    n_markers: int,
    rng: random.Random,
    bias_sd: float = 0.02,
) -> list[float]:
    """Generate per-marker capture/amplification biases.

    Each marker gets a fixed bias drawn from N(0, bias_sd). This models
    the systematic reference/alt allele capture efficiency difference seen
    in real hybridisation capture and amplicon data (Vynck et al.).

    A bias of +0.02 means the observed ALT VAF is shifted +0.02 relative
    to the true VAF (i.e. ALT allele is preferentially captured).

    Args:
        n_markers: Number of markers.
        rng: Random instance for reproducibility.
        bias_sd: Standard deviation of the bias distribution. Typical values:
            0.0 = no bias (ideal), 0.02 = realistic (empirically measured as
            0.019 on 76-SNP rhAmpSeq panel across 210 joint-called VCFs),
            0.05 = high (poor panel design).

    Returns:
        List of per-marker bias values.
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
    """Generate biases with a heavy-tailed distribution.

    The empirical bias distribution is heavy-tailed: median |bias| is 0.005
    but 95th percentile is 0.041 and max is 0.10. A simple Gaussian
    underestimates the tails. This uses a mixture model:

        95% of markers: N(0, sd)         — typical markers
        5% of markers:  N(0, outlier_sd) — outlier markers with extreme bias

    The default parameters are calibrated from 71 markers across 210
    joint-called VCFs on a 76-SNP rhAmpSeq panel, yielding an overall
    SD of ~0.018 matching the empirical measurement.

    Args:
        n_markers: Number of markers.
        rng: Random instance for reproducibility.
        sd: Standard deviation for the bulk of markers (default 0.012).
        outlier_frac: Fraction of outlier markers (default 0.05).
        outlier_sd: Standard deviation for outlier markers (default 0.08).

    Returns:
        List of per-marker bias values.
    """
    biases = []
    for _ in range(n_markers):
        if rng.random() < outlier_frac:
            biases.append(rng.gauss(0, outlier_sd))
        else:
            biases.append(rng.gauss(0, sd))
    return biases


# ---------------------------------------------------------------------------
# Relatedness-based genotype generation
# ---------------------------------------------------------------------------

# IBD sharing probabilities: (P(IBD=0), P(IBD=1), P(IBD=2))
RELATEDNESS_IBD = {
    "unrelated": (1.0, 0.0, 0.0),
    "cousin": (0.75, 0.25, 0.0),  # first cousins: 1/8 kinship
    "half-sibling": (0.5, 0.5, 0.0),  # half-siblings: 1/4 kinship
    "parent-child": (0.0, 1.0, 0.0),  # parent-child: always share 1 allele
    "sibling": (0.25, 0.5, 0.25),  # full siblings: 1/4 kinship
}


def _draw_genotype(p_alt: float, rng: random.Random) -> tuple[int, int]:
    """Draw a diploid genotype from Hardy-Weinberg equilibrium.

    Args:
        p_alt: Population ALT allele frequency.
        rng: Random instance.

    Returns:
        Diploid genotype as (allele1, allele2), each 0 or 1.
    """
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

    Args:
        host_gt: Host diploid genotype.
        p_alt: Population ALT allele frequency.
        ibd_probs: (P(IBD=0), P(IBD=1), P(IBD=2)).
        rng: Random instance.

    Returns:
        Donor diploid genotype.
    """
    r = rng.random()
    if r < ibd_probs[0]:
        # IBD=0: independent draw
        return _draw_genotype(p_alt, rng)
    elif r < ibd_probs[0] + ibd_probs[1]:
        # IBD=1: share one allele, draw the other independently
        shared = host_gt[rng.randint(0, 1)]
        other = 1 if rng.random() < p_alt else 0
        return (shared, other) if rng.random() < 0.5 else (other, shared)
    else:
        # IBD=2: identical genotype
        return host_gt


def generate_related_genotypes(
    n_markers: int,
    relatedness: str,
    rng: random.Random,
    maf_range: tuple[float, float] = (0.2, 0.5),
) -> list[dict]:
    """Generate synthetic host-donor genotype pairs with specified relatedness.

    Marker allele frequencies are drawn uniformly from ``maf_range``, then
    host and donor genotypes are generated with appropriate IBD sharing.

    Args:
        n_markers: Number of markers to generate.
        relatedness: One of 'unrelated', 'cousin', 'half-sibling',
            'parent-child', 'sibling'.
        rng: Random instance for reproducibility.
        maf_range: (min, max) minor allele frequency range for markers.

    Returns:
        List of dicts with keys: chrom, pos, ref, alt, host_gt, donor_gt,
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


def _mendelian_child(
    parent1: tuple[int, int],
    parent2: tuple[int, int],
    rng: random.Random,
) -> tuple[int, int]:
    """Draw a child genotype by Mendelian segregation from two parents.

    Each parent transmits one allele (chosen uniformly at random).

    Args:
        parent1: First parent diploid genotype.
        parent2: Second parent diploid genotype.
        rng: Random instance.

    Returns:
        Child diploid genotype (sorted so smaller allele first).
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

    For each marker:
    1. Draw population ALT allele frequency
    2. Draw two parent genotypes from Hardy-Weinberg
    3. Derive each sibling independently by Mendelian segregation

    This preserves the correct 3-way sibling correlation structure:
    each pair has IBD distribution (0.25, 0.5, 0.25) and the three
    genotypes are correlated through shared parents.

    Args:
        n_markers: Number of biallelic markers to generate.
        rng: Random instance for reproducibility.
        maf_range: (min, max) minor allele frequency range for markers.

    Returns:
        List of dicts with keys: chrom, pos, ref, alt, host_gt, donor1_gt,
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

    Args:
        markers: List of marker dicts from generate_related_genotypes().
        path: Output VCF path.
        sample_name: Sample name for VCF header.
        key: Which genotype to write ('host_gt' or 'donor_gt').
        depth: Simulated read depth for FORMAT fields.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=248956422>\n")
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
        f.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
        f.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">\n')
        f.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        for m in markers:
            gt = m[key]
            gt_str = f"{gt[0]}/{gt[1]}"
            n_alt = alt_dose(gt)
            ad_alt = round(depth * n_alt / 2)
            ad_ref = depth - ad_alt
            sample_field = f"{gt_str}:{ad_ref},{ad_alt}:{depth}:99"
            f.write(
                f"{m['chrom']}\t{m['pos']}\t.\t{m['ref']}\t{m['alt']}\t"
                f".\tPASS\t.\tGT:AD:DP:GQ\t{sample_field}\n"
            )


def blend_vcfs(
    host_path: str | Path,
    donor_path: str | Path,
    donor_fraction: float,
    target_depth: int | None = None,
    sample_name: str | None = None,
    seed: int | None = None,
    marker_bias_sd: float = 0.0,
    fixed_biases: list[float] | None = None,
    error_rate: float = 0.01,
    allele_dropout_rate: float = 0.0,
    locus_dropout_rate: float = 0.0,
    depth_cv: float = 0.0,
    realistic_biases: bool = False,
) -> BlendResult:
    """Blend two genotype VCFs to create a synthetic chimeric VCF.

    Args:
        host_path: Path to host genotype VCF (plain text).
        donor_path: Path to donor genotype VCF (plain text).
        donor_fraction: Fraction of donor DNA (0.0 to 1.0).
        target_depth: Fixed depth for all markers. If None, uses actual
            depth from the host VCF.
        sample_name: Sample name for the output VCF column header.
            Defaults to 'simulated'.
        seed: Random seed for reproducibility.
        marker_bias_sd: Standard deviation of per-marker capture bias.
            0.0 = no bias (ideal simulation), 0.02 = realistic.
            Ignored when ``fixed_biases`` is provided. When
            ``realistic_biases`` is True, this is also ignored (the
            realistic mixture model is used instead).
        fixed_biases: Pre-generated per-marker bias values. When provided,
            these biases are used directly instead of generating random ones.
            Length must match the number of shared markers.
        error_rate: Per-read sequencing error rate (default 0.01 = 1%).
            Each read has this probability of being mis-called.
        allele_dropout_rate: Per-marker probability of allele dropout
            (0.0–1.0). When dropout occurs at a heterozygous site, one
            allele's reads are entirely lost, making a het look like a hom.
        locus_dropout_rate: Per-marker probability of complete locus
            dropout (0.0–1.0). Affected markers produce zero reads.
        depth_cv: Coefficient of variation for per-marker depth
            (0.0 = uniform depth, 0.43 = empirical value from rhAmpSeq).
            When > 0 and target_depth is set, per-marker depths are drawn
            from a log-normal distribution.
        realistic_biases: If True, use the heavy-tailed mixture bias
            distribution (generate_marker_biases_realistic) instead of
            a simple Gaussian. Ignored when ``fixed_biases`` is provided.

    Returns:
        BlendResult containing the header, VCF record lines, and statistics.
    """
    if not 0.0 <= donor_fraction <= 1.0:
        raise ValueError(f"donor_fraction must be 0.0-1.0, got {donor_fraction}")

    rng = random.Random(seed)
    host_header, host_records = parse_vcf(host_path)
    _, donor_records = parse_vcf(donor_path)

    # Index donor records by locus
    donor_by_locus: dict[str, VcfRecord] = {}
    for rec in donor_records:
        donor_by_locus[rec.locus] = rec

    if sample_name is None:
        sample_name = "simulated"

    # Rewrite the header: replace sample name in #CHROM line
    out_header = []
    for line in host_header:
        if line.startswith("#CHROM"):
            parts = line.split("\t")
            parts[-1] = sample_name
            out_header.append("\t".join(parts))
        else:
            out_header.append(line)

    # Pre-generate per-marker biases (one per shared locus, in order)
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

    # Pre-generate per-marker depths
    if depth_cv > 0 and target_depth is not None:
        marker_depths = sample_marker_depths(n_shared, target_depth, depth_cv, rng)
    else:
        marker_depths = None  # use flat target_depth or host depth

    out_records: list[str] = []
    bias_info: list[tuple[str, int, str, str, float]] = []
    num_markers = 0
    num_informative = 0
    bias_idx = 0

    for host_rec in host_records:
        donor_rec = donor_by_locus.get(host_rec.locus)
        if donor_rec is None:
            continue

        host_gt = extract_gt(host_rec)
        donor_gt = extract_gt(donor_rec)
        if host_gt is None or donor_gt is None:
            continue

        # Must share the same REF allele
        if host_rec.ref != donor_rec.ref:
            continue

        num_markers += 1
        if is_informative(host_gt, donor_gt):
            num_informative += 1

        # Determine depth
        if marker_depths is not None:
            depth = marker_depths[bias_idx]
        elif target_depth is not None:
            depth = target_depth
        else:
            depth = extract_depth(host_rec) or 1000

        # Locus dropout: marker produces zero reads
        if locus_dropout_rate > 0 and rng.random() < locus_dropout_rate:
            bias_idx += 1
            continue

        # Calculate expected VAF, apply per-marker capture bias, then sample
        vaf = expected_vaf(host_gt, donor_gt, donor_fraction)
        this_bias = marker_biases[bias_idx]
        vaf_biased = max(0.0, min(1.0, vaf + this_bias))
        alt_allele_bias = host_rec.alt if host_rec.alt != "." else donor_rec.alt
        bias_info.append((host_rec.chrom, host_rec.pos, host_rec.ref, alt_allele_bias, this_bias))
        bias_idx += 1

        # Allele dropout: at a het-like site, one allele is entirely lost.
        # This pushes the observed VAF to 0.0 or 1.0.
        if allele_dropout_rate > 0 and 0.05 < vaf_biased < 0.95:
            if rng.random() < allele_dropout_rate:
                vaf_biased = 0.0 if rng.random() < 0.5 else 1.0

        ref_count, alt_count = sample_allele_counts(vaf_biased, depth, rng, error_rate)

        # Build output record
        gt = gt_from_counts(ref_count, alt_count)
        total = ref_count + alt_count
        af_val = f"{alt_count / total:.4f}" if total > 0 else "0"
        gq = 99
        # Simplified PL: just put 0 for the called genotype
        if gt == "0/0":
            pl = f"0,{gq},{gq * 10}"
        elif gt == "1/1":
            pl = f"{gq * 10},{gq},0"
        else:
            pl = f"{gq * 5},0,{gq * 5}"

        sample_field = f"{gt}:{ref_count},{alt_count}:{total}:{gq}:{pl}:{af_val}"
        format_field = "GT:AD:DP:GQ:PL:AF"

        # Use the host's ALT allele; if host was ref-only, use donor's ALT
        alt_allele = host_rec.alt if host_rec.alt != "." else donor_rec.alt
        if alt_allele == ".":
            # Both are hom-ref with no ALT listed -- skip or emit as-is
            # For chimerism we still want to emit the site for completeness
            alt_allele = "."
            # If no ALT allele, strip AF from format
            sample_field = f"{gt}:{ref_count}:{total}:{gq}:{pl}"
            format_field = "GT:AD:DP:GQ:PL"

        # Build minimal INFO
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

    return BlendResult(
        header=out_header,
        records=out_records,
        num_markers=num_markers,
        num_informative=num_informative,
        marker_biases=bias_info
        if (marker_bias_sd > 0 or fixed_biases is not None or realistic_biases)
        else None,
    )


def write_vcf(result: BlendResult, path: str | Path) -> None:
    """Write a BlendResult to a VCF file.

    Args:
        result: BlendResult from blend_vcfs().
        path: Output file path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for line in result.header:
            fh.write(line + "\n")
        for line in result.records:
            fh.write(line + "\n")


def blend_from_genotype_dicts(
    markers: list[dict],
    donor_fractions: list[float],
    target_depth: int = 1000,
    seed: int | None = None,
    error_rate: float = 0.01,
    depth_cv: float = 0.0,
    sample_name: str = "simulated",
) -> BlendResult:
    """Create a synthetic chimeric VCF directly from genotype dicts.

    Designed for use with generate_sibling_trio_genotypes() output.
    Supports 1 or 2 donors via donor_fractions length.

    Args:
        markers: List of marker dicts with host_gt, donor1_gt, donor2_gt.
        donor_fractions: [f_donor1] or [f_donor1, f_donor2].
        target_depth: Mean sequencing depth.
        seed: Random seed.
        error_rate: Sequencing error rate.
        depth_cv: Depth CV across markers.
        sample_name: Sample name for VCF header.

    Returns:
        BlendResult with synthetic chimeric VCF data.
    """
    if sum(donor_fractions) > 1.0 + 1e-9:
        raise ValueError(f"donor_fractions sum to {sum(donor_fractions):.4f}, must be <= 1.0")

    rng = random.Random(seed)
    n = len(markers)

    if depth_cv > 0:
        depths = sample_marker_depths(n, target_depth, depth_cv, rng)
    else:
        depths = [target_depth] * n

    # Build header
    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=248956422>",
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
        ref_count, alt_count = sample_allele_counts(vaf, depths[i], rng, error_rate)

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
