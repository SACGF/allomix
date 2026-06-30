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

    Prefer a freshly joint-called pipeline dir; fall back to the committed
    snapshot so the paper builds from a fresh checkout.
    """
    if any(SRP434573_PIPELINE_GENOTYPES.glob("*.SRP434573.vcf.gz")):
        return SRP434573_PIPELINE_GENOTYPES
    return SRP434573_COMMITTED_GENOTYPES


def resolve_srp434573_synthetic_dir() -> Path | None:
    """Locate the semi-synthetic mixture VCFs, or None if not generated yet.

    Prefer the freshly joint-called pipeline dir, else the committed snapshot.
    None when neither exists, so a fresh checkout (before the TAU-side generation
    step has run) builds as before, without the synthetic points.
    """
    if any(SRP434573_PIPELINE_SYNTHETIC.glob("*.synthetic.admix.vcf.gz")):
        return SRP434573_PIPELINE_SYNTHETIC
    if any(SRP434573_COMMITTED_SYNTHETIC.glob("*.synthetic.admix.vcf.gz")):
        return SRP434573_COMMITTED_SYNTHETIC
    return None
