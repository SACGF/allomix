"""Shared helpers for the SRP434573 paper scripts."""

from pathlib import Path

# Committed snapshot of the joint-called SRP434573 genotype VCFs (issue #21), and
# the location a full from-scratch pipeline run writes them to.
SRP434573_COMMITTED_GENOTYPES = Path("paper/public_data/SRP434573/genotypes")
SRP434573_PIPELINE_GENOTYPES = Path("output/genotypes/SRP434573")

# Semi-synthetic sub-0.5% mixtures (issue #5): the two pure reference BAMs of
# each pair blended with samtools subsample, joint-called the same way and
# committed alongside the real snapshot. Generated TAU-side by
# paper/scripts/make_semisynthetic_srp434573.py.
SRP434573_COMMITTED_SYNTHETIC = Path("paper/public_data/SRP434573/genotypes_synthetic")
SRP434573_PIPELINE_SYNTHETIC = Path("output/genotypes/SRP434573_synthetic")


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


def resolve_srp434573_synthetic_dir() -> Path | None:
    """Locate the semi-synthetic mixture VCFs, or ``None`` if not generated yet.

    Prefer a freshly joint-called ``output/genotypes/SRP434573_synthetic`` if it
    has admix VCFs; otherwise fall back to the committed snapshot under
    ``paper/public_data/SRP434573/genotypes_synthetic``. Returns ``None`` when
    neither exists, so a fresh checkout (before the TAU-side generation step has
    run) builds exactly as before, without the synthetic points.

    Returns:
        Path to the directory holding the per-pair ``*.synthetic.admix.vcf.gz``
        files, or ``None``.
    """
    if any(SRP434573_PIPELINE_SYNTHETIC.glob("*.synthetic.admix.vcf.gz")):
        return SRP434573_PIPELINE_SYNTHETIC
    if any(SRP434573_COMMITTED_SYNTHETIC.glob("*.synthetic.admix.vcf.gz")):
        return SRP434573_COMMITTED_SYNTHETIC
    return None
