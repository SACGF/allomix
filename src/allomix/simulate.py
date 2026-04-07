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
            records.append(VcfRecord(
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
            ))
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
) -> tuple[int, int]:
    """Sample allele counts from a binomial distribution.

    Args:
        vaf: Expected variant allele frequency.
        depth: Total read depth to simulate.
        rng: Optional Random instance for reproducibility.

    Returns:
        Tuple of (ref_count, alt_count).
    """
    if rng is None:
        rng = random.Random()
    if depth <= 0:
        return (0, 0)
    # Clamp VAF to valid probability range
    p = max(0.0, min(1.0, vaf))
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
            0.0 = no bias (ideal), 0.02 = moderate (typical for capture panels),
            0.05 = high (poor panel design).

    Returns:
        List of per-marker bias values.
    """
    if bias_sd <= 0:
        return [0.0] * n_markers
    return [rng.gauss(0.0, bias_sd) for _ in range(n_markers)]


def blend_vcfs(
    host_path: str | Path,
    donor_path: str | Path,
    donor_fraction: float,
    target_depth: int | None = None,
    sample_name: str | None = None,
    seed: int | None = None,
    marker_bias_sd: float = 0.0,
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
    marker_biases = generate_marker_biases(n_shared, rng, marker_bias_sd)

    out_records: list[str] = []
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
        if target_depth is not None:
            depth = target_depth
        else:
            depth = extract_depth(host_rec) or 1000

        # Calculate expected VAF, apply per-marker capture bias, then sample
        vaf = expected_vaf(host_gt, donor_gt, donor_fraction)
        vaf_biased = max(0.0, min(1.0, vaf + marker_biases[bias_idx]))
        bias_idx += 1
        ref_count, alt_count = sample_allele_counts(vaf_biased, depth, rng)

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

        line = "\t".join([
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
        ])
        out_records.append(line)

    return BlendResult(
        header=out_header,
        records=out_records,
        num_markers=num_markers,
        num_informative=num_informative,
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
