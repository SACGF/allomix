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

# Sex inference from non-PAR chrX heterozygosity (allomix.qc.sex), shared with
# the calibration script. Starting values were set from the 7 public SRP434573
# individuals (M1-M4, F1-F3): females showed a chrX het rate of 0.44-0.71 and
# males 0.00-0.04 after a GQ >= 20 filter. They are to be finalised from
# ``scripts/sex_calibration_summary.py`` on the internal cohort.
#
# Female model: het probability = FEMALE_X_HET_SCALE * the sample's own
# autosomal het rate (the panel's design het rate), clamped to [0.05, 0.95].
FEMALE_X_HET_SCALE = 1.0
# Male model: residual spurious-het probability on a hemizygous chromosome
# (GATK diploid calls at a male chrX site that pass GQ).
MALE_X_SPURIOUS_HET = 0.02
# |log10 likelihood ratio| (female vs male) needed for a confident call.
SEX_LR_THRESHOLD = 2.0
# Fewer usable non-PAR chrX sites than this -> UNAVAILABLE (no inference).
MIN_X_SITES = 5
# Minimum reference-sample depth at a chrX / autosomal site for it to count in
# sex inference. Deliberately not the admixture ``--min-dp`` (default 100): GATK
# downsamples per alignment start, so amplicon reference VCFs cap DP near 50
# (SRP434573 does), and the admixture depth floor would discard every site.
SEX_MIN_REF_DP = 20

# Sex-mismatch cross-check (allomix.qc.sex_mismatch, #48): the independent
# sex-chromosome estimate of the donor fraction for a sex-mismatched host/donor
# pair, compared with the autosomal MLE in qc.assess_quality.
#
# Fewer usable non-PAR chrX markers than this and the chrX copy-number readout is
# not reported (basis falls back to None).
SEXCHROM_MIN_MARKERS = 5
# Fewer non-PAR chrY depth sites than this and the chrY depth readout is skipped
# in favour of chrX (and when it is the only readout, its CI is the whole [0, 1]).
SEXCHROM_MIN_CHRY_SITES = 3
# The two estimates are concordant when their 95% CIs overlap or they differ by
# at most this many percentage points. Provisional: 2.0 pp is a placeholder set
# from the first in-silico and public-data runs and is to be finalised from the
# wetlab validation (plan: plans/sex_markers.md, Phase 3).
SEXCHROM_CONCORDANCE_PP = 2.0
