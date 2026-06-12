"""Shared helpers for the SRP434573 paper scripts."""

from pathlib import Path

# Committed snapshot of the joint-called SRP434573 genotype VCFs (issue #21), and
# the location a full from-scratch pipeline run writes them to.
SRP434573_COMMITTED_GENOTYPES = Path("paper/public_data/SRP434573/genotypes")
SRP434573_PIPELINE_GENOTYPES = Path("output/genotypes/SRP434573")


def resolve_srp434573_genotypes_dir() -> Path:
    """Locate the SRP434573 genotype VCFs.

    Prefer a freshly joint-called ``output/genotypes/SRP434573`` (full
    from-scratch reproduction) if it has VCFs; otherwise fall back to the
    committed snapshot under ``paper/public_data/SRP434573/genotypes`` so the
    paper builds from a fresh checkout.

    Returns:
        Path to the directory holding the per-mixture ``*.SRP434573.vcf.gz`` and
        ``*.admix.vcf.gz`` files.
    """
    if any(SRP434573_PIPELINE_GENOTYPES.glob("*.SRP434573.vcf.gz")):
        return SRP434573_PIPELINE_GENOTYPES
    return SRP434573_COMMITTED_GENOTYPES
