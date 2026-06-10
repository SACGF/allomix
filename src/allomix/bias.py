"""Per-marker amplification bias estimation and correction.

Implements bias correction based on the observation that capture/amplification
panels introduce systematic per-marker shifts in observed VAF relative to the
true allele frequency. These biases are consistent across samples sequenced
with the same panel (Vynck et al.).

Bias is estimated as the median deviation of observed heterozygous VAF from 0.5
at each marker across a training set of samples (positive = ALT-favoured). The
correction is integrated directly into the MLE likelihood by adjusting the
expected reference allele weight for each marker. The adjustment is applied
multiplicatively in logit space, not as a flat additive shift, so it stays
valid at informative markers whose expected VAF is far from 0.5 (issue #20);
see ``allomix.chimerism.apply_bias``.
"""

import csv
import statistics
from dataclasses import dataclass
from pathlib import Path

from allomix.genotype import MarkerData, MarkerKey, marker_key


@dataclass
class MarkerBias:
    """Per-marker bias estimate."""

    chrom: str
    pos: int
    ref: str
    alt: str
    bias: float  # median(observed_het_VAF - 0.5); positive = ALT-favoured
    n_het: int  # number of het observations used


def estimate_biases(
    marker_lists: list[list[MarkerData]],
    min_het: int = 1,
) -> dict[MarkerKey, MarkerBias]:
    """Estimate per-marker amplification bias from heterozygous observations.

    For each marker, collects VAF from all samples where the genotype is
    heterozygous (0/1). Bias is estimated as median(VAF - 0.5).

    A positive bias means the ALT allele is preferentially captured/amplified.

    Args:
        marker_lists: List of MarkerData lists, one per training sample.
        min_het: Minimum number of het observations required to estimate
            bias at a marker. Markers with fewer observations are excluded.

    Returns:
        Dict mapping marker key to MarkerBias.
    """
    # Collect het VAF deviations per marker
    het_deviations: dict[MarkerKey, list[float]] = {}
    marker_info: dict[MarkerKey, tuple[str, int, str, str]] = {}

    for markers in marker_lists:
        for m in markers:
            if m.gt != (0, 1):
                continue
            dp = m.ad_ref + m.ad_alt
            if dp <= 0:
                continue
            vaf = m.ad_alt / dp
            key = marker_key(m)
            het_deviations.setdefault(key, []).append(vaf - 0.5)
            marker_info[key] = (m.chrom, m.pos, m.ref, m.alt)

    biases: dict[MarkerKey, MarkerBias] = {}
    for key, devs in het_deviations.items():
        if len(devs) < min_het:
            continue
        chrom, pos, ref, alt = marker_info[key]
        biases[key] = MarkerBias(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            bias=statistics.median(devs),
            n_het=len(devs),
        )

    return biases


def estimate_biases_both_het(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix_lists: list[list[MarkerData]],
    min_het: int = 1,
    min_dp: int = 1,
) -> dict[MarkerKey, MarkerBias]:
    """Estimate per-marker bias from admix samples at both-het markers.

    At markers where the host and every donor are heterozygous (0/1), the true
    admixture ALT VAF is 0.5 regardless of the mixing fraction. The observed
    admix VAF there therefore gives the per-marker bias directly, from the same
    data and the same caller being analysed. This avoids needing a separate
    panel pileup, which the two-phase pipeline (GATK panel, mpileup admix) does
    not produce, and sidesteps the caller-mismatch footgun of estimating from a
    differently-called panel VCF (issue #11).

    Each admix sample contributes one observation per both-het marker; bias is
    the median of ``observed VAF - 0.5`` across observations.

    A pair's both-het markers are non-informative for that same pair (host and
    donor share the heterozygous genotype), so the resulting table only helps
    other pairs whose informative markers it covers. This is therefore a cohort
    table builder, not an inline single-run correction: pool across patients and
    apply the table with ``--bias-table``.

    Args:
        host: Markers from the host genotyping VCF.
        donors: List of marker lists, one per donor genotyping VCF.
        admix_lists: One marker list per admix sample to estimate from (pool
            across a cohort for useful marker coverage).
        min_het: Minimum observations required to keep a marker.
        min_dp: Minimum admix depth for an observation to count.

    Returns:
        Dict mapping marker key to MarkerBias (same het-site VAF convention as
        ``estimate_biases``).
    """
    host_idx = {marker_key(m): m for m in host}
    donor_idxs = [{marker_key(m): m for m in d} for d in donors]

    # Markers where host and every donor are heterozygous.
    both_het: dict[MarkerKey, tuple[str, int, str, str]] = {}
    for key, h in host_idx.items():
        if h.gt != (0, 1):
            continue
        donor_markers = [di.get(key) for di in donor_idxs]
        if any(d is None or d.gt != (0, 1) for d in donor_markers):
            continue
        both_het[key] = (h.chrom, h.pos, h.ref, h.alt)

    deviations: dict[MarkerKey, list[float]] = {}
    for admix in admix_lists:
        for m in admix:
            key = marker_key(m)
            if key not in both_het:
                continue
            dp = m.ad_ref + m.ad_alt
            if dp < min_dp or dp <= 0:
                continue
            deviations.setdefault(key, []).append(m.ad_alt / dp - 0.5)

    biases: dict[MarkerKey, MarkerBias] = {}
    for key, devs in deviations.items():
        if len(devs) < min_het:
            continue
        chrom, pos, ref, alt = both_het[key]
        biases[key] = MarkerBias(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            bias=statistics.median(devs),
            n_het=len(devs),
        )

    return biases


def save_bias_table(biases: dict[MarkerKey, MarkerBias], path: Path | str) -> None:
    """Write bias estimates to a TSV file.

    Args:
        biases: Dict mapping marker key to MarkerBias.
        path: Output file path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chrom", "pos", "ref", "alt", "bias", "n_het"])
        for key in sorted(biases.keys()):
            mb = biases[key]
            writer.writerow([mb.chrom, mb.pos, mb.ref, mb.alt, f"{mb.bias:.6f}", mb.n_het])


def load_bias_table(path: Path | str) -> dict[MarkerKey, float]:
    """Load a bias table TSV and return a dict of marker key -> bias value.

    Args:
        path: Path to bias table TSV file.

    Returns:
        Dict mapping (chrom, pos, ref, alt) to bias float.
    """
    biases: dict[MarkerKey, float] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key: MarkerKey = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            biases[key] = float(row["bias"])
    return biases


def biases_to_simple_dict(biases: dict[MarkerKey, MarkerBias]) -> dict[MarkerKey, float]:
    """Convert MarkerBias dict to a simple key -> float dict for use in estimation.

    Args:
        biases: Dict mapping marker key to MarkerBias.

    Returns:
        Dict mapping marker key to bias float value.
    """
    return {key: mb.bias for key, mb in biases.items()}
