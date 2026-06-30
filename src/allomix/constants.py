"""Constants shared across allomix modules.

A constant lives here once it is used in two or more modules. Constants used
in a single module stay local to that module. This module imports nothing from
the package, so anything can import it without risk of a circular import.
"""

# Diploid ploidy. Named so the dosage math (ref_dose = PLOIDY - alt_dose,
# weight /= PLOIDY) reads as ploidy, not a bare 2 confusable with the
# donor-count cap in the multi-donor estimator.
PLOIDY = 2

# A sequencing error changes the true base into one of the 3 other bases, so
# (assuming even spread) a miscall to one specific base, e.g. a true REF read as
# the ALT allele, has probability ``error_rate / N_OTHER_BASES``. The
# per-direction error floor of the 4-state model shared by chimerism, qc, detect,
# and simulate.
N_OTHER_BASES = len("ACGT") - 1

# Default robust-refit residual cut for the median/MAD outlier filter, in robust
# SDs. 3.5 leaves clean data essentially untouched (drops <1% of markers by
# chance) while removing copy-number / LoH-inconsistent markers. Both the
# estimator default (chimerism) and the CLI/analysis default for --robust-k.
ROBUST_K_DEFAULT = 3.5

# Pipeline defaults shared by the CLI, the genotype reader, and the estimators.
DEFAULT_ERROR_RATE = 0.01  # per-base sequencing error rate (1%)
DEFAULT_MIN_DP = 100  # minimum admixture read depth at a marker
DEFAULT_MIN_GQ = 20  # minimum host/donor genotype quality

# Confidence level for the profile-likelihood / likelihood-ratio CIs (95%).
# The estimators and the host-presence detector all build CIs at this level;
# two-sided normal quantiles derive from it as ``1 - (1 - CI_LEVEL) / 2``.
CI_LEVEL = 0.95

# VAF cutoffs for calling a genotype from read counts: hom-ref at or below
# HOM_REF_MAX_VAF, hom-alt at or above HOM_ALT_MIN_VAF, het in between. Used by
# the reference-sample consistency check (genotype) and the synthetic caller
# (simulate). Symmetric about 0.5, so the pair is defined from one number.
HOM_REF_MAX_VAF = 0.05
HOM_ALT_MIN_VAF = 1.0 - HOM_REF_MAX_VAF  # 0.95
