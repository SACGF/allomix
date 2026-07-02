"""Detect the variant caller that produced a VCF, from its header (issue #42).

allomix's two-phase workflow expects host/donor genotypes from GATK joint
calling and admix allele depths from bcftools mpileup (see docs/joint_calling.md).
Mixing callers is a foot-gun in two independent ways:

1. A GATK-called admix VCF strips minority ALT reads at hom-ref blocks (the exact
   low-fraction signal chimerism detection needs), so the admix VCF should always
   come from forced ``bcftools mpileup`` at the panel sites.
2. Per-marker amplification bias is caller-specific, so a bias table estimated
   from GATK-called panel het sites and applied to bcftools mpileup admix data
   makes the estimate worse, not better.

This module reads the VCF header (it never opens a BAM) and classifies the caller
so the CLI can warn on these mismatches. Detection is best-effort: a header with
no recognised signature yields ``UNKNOWN`` and the CLI stays silent rather than
guessing.

Note that in the standard two-phase workflow the panel caller (GATK) and the
admix caller (mpileup) differ by design, so the caller-mismatch warning is scoped
to the bias-correction path (bias-table source vs admix), not a blanket
panel-vs-admix comparison.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from cyvcf2 import VCF


class Caller(Enum):
    """Variant caller behind a VCF, as far as the header reveals."""

    MPILEUP = "mpileup"
    GATK = "gatk"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class CallerInfo:
    """Detected caller plus a short human-readable note on what matched."""

    caller: Caller
    evidence: str


#: INFO IDs emitted by GATK's variant annotators but not by bcftools mpileup.
#: Used only as a fallback fingerprint when provenance header lines are stripped.
_GATK_INFO = frozenset(
    {"MQRankSum", "ReadPosRankSum", "FS", "SOR", "ExcessHet", "InbreedingCoeff", "MLEAC", "MLEAF"}
)
#: INFO IDs specific to bcftools mpileup output.
_MPILEUP_INFO = frozenset({"DP4", "I16"})

_EXPLICIT = re.compile(r"^##allomixCaller=(\S+)", re.IGNORECASE)
_INFO_ID = re.compile(r"^##INFO=<ID=([^,>]+)")
_GATK_SOURCE_TOOLS = ("HaplotypeCaller", "GenotypeGVCFs", "CombineGVCFs")


def caller_from_token(token: str | None) -> Caller:
    """Map a recorded caller token (e.g. from a bias table) to a ``Caller``."""
    if token is None:
        return Caller.UNKNOWN
    try:
        return Caller(token.strip().lower())
    except ValueError:
        return Caller.UNKNOWN


def _info_ids(lines: list[str]) -> set[str]:
    ids: set[str] = set()
    for line in lines:
        m = _INFO_ID.match(line)
        if m:
            ids.add(m.group(1))
    return ids


def detect_caller_from_header(header: str) -> CallerInfo:
    """Classify the caller from a raw VCF header string.

    Precedence, most reliable first:
      1. An explicit ``##allomixCaller=`` stamp.
      2. Provenance command / ``##source`` lines (bcftools mpileup, GATK).
      3. Caller-specific INFO fingerprints, for headers whose command lines were
         stripped (``DP4``/``I16`` => mpileup; GATK annotator IDs => GATK).
    Returns ``UNKNOWN`` when nothing matches.
    """
    lines = header.splitlines()

    # 1. Explicit allomix stamp wins when it names a caller we recognise.
    for line in lines:
        m = _EXPLICIT.match(line)
        if m:
            token = m.group(1).strip().lower()
            explicit = caller_from_token(token)
            if explicit is not Caller.UNKNOWN:
                return CallerInfo(explicit, f"##allomixCaller={token}")
            break  # recognised header key but unknown value: fall through to sniffing

    # 2. Provenance lines. A bcftools command line that ran mpileup is definitive;
    # a GATKCommandLine or GATK ##source is definitive for GATK. In the standard
    # workflow these are mutually exclusive, so the order between them is moot.
    for line in lines:
        low = line.lower()
        if low.startswith("##bcftools") and "command=" in low and "mpileup" in low:
            return CallerInfo(Caller.MPILEUP, "bcftools mpileup command line")
    for line in lines:
        if line.startswith("##GATKCommandLine"):
            return CallerInfo(Caller.GATK, "GATKCommandLine header")
        if line.startswith("##source=") and any(t in line for t in _GATK_SOURCE_TOOLS):
            return CallerInfo(Caller.GATK, line.strip())

    # 3. INFO-field fingerprints (command lines stripped).
    info_ids = _info_ids(lines)
    mpileup_hit = info_ids & _MPILEUP_INFO
    if mpileup_hit:
        return CallerInfo(Caller.MPILEUP, f"{', '.join(sorted(mpileup_hit))} INFO field(s)")
    gatk_hit = info_ids & _GATK_INFO
    if gatk_hit:
        return CallerInfo(Caller.GATK, f"{', '.join(sorted(gatk_hit))} INFO field(s)")

    return CallerInfo(Caller.UNKNOWN, "no caller signature in header")


def detect_caller(vcf_path: str | Path) -> CallerInfo:
    """Detect the caller of a VCF file from its header."""
    vcf = VCF(str(vcf_path))
    try:
        header = vcf.raw_header
    finally:
        vcf.close()
    return detect_caller_from_header(header)


__all__ = [
    "Caller",
    "CallerInfo",
    "caller_from_token",
    "detect_caller",
    "detect_caller_from_header",
]
