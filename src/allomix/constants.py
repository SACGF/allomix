"""Constants shared across allomix modules.

A constant lives here once it is used in two or more modules. Constants used
in a single module stay local to that module. This module imports nothing from
the package, so anything can import it without risk of a circular import.
"""

# DNA has 4 bases, so a sequencing error changes the true base into one of the
# 3 other bases. Assuming errors are spread evenly, a miscall to one specific
# base (e.g. a true REF read as the ALT allele) has probability
# ``error_rate / N_OTHER_BASES``. This is the per-direction error floor of the
# 4-state model shared by chimerism, qc, detect, and simulate.
N_OTHER_BASES = len("ACGT") - 1
