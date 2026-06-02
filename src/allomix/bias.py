"""Per-marker amplification bias estimation and correction.

Implements bias correction based on the observation that capture/amplification
panels introduce systematic per-marker shifts in observed VAF relative to the
true allele frequency. These biases are consistent across samples sequenced
with the same panel (Vynck et al.).

Bias is estimated as the median deviation of observed heterozygous VAF from 0.5
at each marker across a training set of samples. The correction is integrated
directly into the MLE likelihood by adjusting the expected reference allele
weight for each marker.
"""

from __future__ import annotations

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
