"""Sex and host/donor pair-status enums shared across the package.

These live in a leaf module (no allomix imports) because they are needed on
both sides of an import boundary: ``allomix.qc.sex`` infers them from
``genotype.MarkerData``, and ``genotype.classify_markers`` routes non-PAR chrX
markers on the pair status (#46). ``allomix.qc.sex`` re-exports both, so
``from allomix.qc.sex import Sex, PairStatus`` keeps working.
"""

from __future__ import annotations

from enum import Enum


class Sex(str, Enum):
    """Sex of a reference sample, inferred or declared.

    Member value is the wire string used in TSV/JSON output and on the CLI, so a
    ``Sex`` is a ``str``. ``UNAVAILABLE`` means the panel had no usable non-PAR
    chrX content; ``AMBIGUOUS`` means there was content and it did not resolve.
    """

    FEMALE = "F"
    MALE = "M"
    AMBIGUOUS = "ambiguous"
    UNAVAILABLE = "unavailable"

    @property
    def confident(self) -> bool:
        """True for a concrete sex (female or male)."""
        return self is Sex.FEMALE or self is Sex.MALE


class PairStatus(str, Enum):
    """Host-vs-donor(s) sex relationship, on effective sexes.

    ``MATCHED_FEMALE`` / ``MATCHED_MALE`` when every party has the same confident
    sex, ``MISMATCHED`` when the sexes differ, ``UNKNOWN`` when any party is
    ambiguous or unavailable.
    """

    MATCHED_FEMALE = "matched_female"
    MATCHED_MALE = "matched_male"
    MISMATCHED = "mismatched"
    UNKNOWN = "unknown"


__all__ = ["PairStatus", "Sex"]
