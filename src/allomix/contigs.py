"""Contig classification for marker routing.

Every marker in allomix is keyed by ``(chrom, pos, ref, alt)`` and nothing else
interprets a coordinate, with one exception: deciding whether a chrX or chrY
site is pseudoautosomal. This module owns that decision and the broader
classification of a contig as autosome, sex chromosome, mitochondrion, or
"other" (alternate haplotypes, decoys, unplaced and random scaffolds, HLA
contigs), so the rest of the package can route markers without caring about
contig naming conventions or genome build.

Pseudoautosomal regions (PAR), 1-based inclusive. The constants below are the
EXACT UNION of the GRCh37 and GRCh38 intervals: overlapping intervals are
merged, non-overlapping ones are kept separate. Per-build source values, 1-based
inclusive as NCBI publishes them in each assembly's ``par_align.gff``. The same
data is packaged in ``bioutils.par`` (fork ``davmlaw/bioutils`` branch
``add-par-regions``, datasets ``ncbi-GRCh37`` and ``ncbi-GRCh38``) in interbase
coordinates, i.e. ``[start - 1, end]``; the cross-check test converts back::

                GRCh37                    GRCh38
    X PAR1      60001-2699520             10001-2781479
    X PAR2      154931044-155260560       155701383-156030895
    Y PAR1      10001-2649520             10001-2781479
    Y PAR2      59034050-59363566         56887903-57217415

Union (what this module masks). The PAR1 intervals overlap and merge; the PAR2
intervals do not overlap on either chromosome (GRCh38 moved the chrX PAR2 by
770 kb and the chrY PAR2 by 2.1 Mb), so each PAR2 stays as two intervals::

    X    10001-2781479,  154931044-155260560,  155701383-156030895
    Y    10001-2781479,  56887903-57217415,    59034050-59363566

Why the union, and why no genome build: allomix always excludes PAR markers,
so a site that is PAR in one build and non-PAR in the other can be excluded
in both. Outside the union the PAR status of a position is the same in both
builds, so the classification is a pure function of ``(chrom, pos)`` and the
tool does not need to know which build the VCF is on.

The cost is the positions that are PAR in one build only. On chrX: 10001-60000
and 2699521-2781479 from PAR1 (50 kb and 82 kb), and the GRCh37 PAR2 interval
154931044-155260560, which is non-PAR Xq28 in GRCh38 coordinates (330 kb). The
GRCh38 chrX PAR2 interval lies beyond GRCh37's chrX length, so it costs nothing
there. On chrY: 2649521-2781479 from PAR1 (132 kb), and the GRCh38 PAR2
interval 56887903-57217415, which is non-PAR distal Yq in GRCh37 coordinates
(330 kb); the GRCh37 chrY PAR2 interval lies beyond GRCh38's chrY length. No
sample-ID panel known to us targets any of this, and the exclusions are
counted and reported so a panel that does would show it. The transcription of
the per-build values is covered by a test that cross-checks these constants
against ``bioutils.par`` when the fork is installed (see
``tests/test_contigs.py`` and ``docs/architecture.md``).
"""

from __future__ import annotations

from enum import Enum


class ContigClass(Enum):
    """Classification of a marker's contig for routing decisions.

    ``OTHER`` is any contig that is not chr1-22, X, Y, M or MT: alternate
    haplotypes, decoys, unplaced and random scaffolds, HLA contigs and the like.
    Markers on such contigs are mapping hazards and are excluded from every
    analysis.
    """

    AUTOSOME = "autosome"
    X_NONPAR = "X_nonPAR"
    X_PAR = "X_PAR"
    Y_NONPAR = "Y_nonPAR"
    Y_PAR = "Y_PAR"
    MT = "MT"
    OTHER = "other"


#: Exact union of the GRCh37 and GRCh38 PAR intervals, 1-based inclusive
#: ``(start, end)``, sorted, per normalised chromosome name. Built the same way
#: the cross-check test builds it from ``bioutils.par.get_par_map`` (after
#: converting its interbase intervals to 1-based): pool both builds' intervals
#: and merge any that overlap.
PAR_X: tuple[tuple[int, int], ...] = (
    (10001, 2781479),  # PAR1, merged
    (154931044, 155260560),  # PAR2, GRCh37
    (155701383, 156030895),  # PAR2, GRCh38
)
PAR_Y: tuple[tuple[int, int], ...] = (
    (10001, 2781479),  # PAR1, merged
    (56887903, 57217415),  # PAR2, GRCh38
    (59034050, 59363566),  # PAR2, GRCh37
)
PAR_REGIONS: dict[str, tuple[tuple[int, int], ...]] = {"X": PAR_X, "Y": PAR_Y}

_AUTOSOME_NAMES = frozenset(str(i) for i in range(1, 23))
_MT_NAMES = frozenset({"M", "MT"})


def normalize_chrom(chrom: str) -> str:
    """Strip an optional ``chr`` prefix (any case) and upper-case the name.

    ``chrX``, ``X`` and ``chrx`` all normalise to ``X``; ``chrM`` to ``M``.
    Names that are not primary contigs come back upper-cased but otherwise
    unchanged, so ``chr1_KI270706v1_random`` becomes ``1_KI270706V1_RANDOM``.
    """
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return c.upper()


def _in_par(regions: tuple[tuple[int, int], ...], pos: int) -> bool:
    return any(start <= pos <= end for start, end in regions)


def classify_contig(chrom: str, pos: int) -> ContigClass:
    """Classify a marker's contig from its name and 1-based position.

    A pure function of ``(chrom, pos)``: no genome build is consulted (see the
    module docstring for why the PAR mask does not need one).

    Args:
        chrom: Contig name as it appears in the VCF, with or without a ``chr``
            prefix, any case.
        pos: 1-based position, used only to split X and Y into PAR and non-PAR.

    Returns:
        The ``ContigClass`` for the site.
    """
    c = normalize_chrom(chrom)
    if c in _AUTOSOME_NAMES:
        return ContigClass.AUTOSOME
    if c == "X":
        return ContigClass.X_PAR if _in_par(PAR_REGIONS["X"], pos) else ContigClass.X_NONPAR
    if c == "Y":
        return ContigClass.Y_PAR if _in_par(PAR_REGIONS["Y"], pos) else ContigClass.Y_NONPAR
    if c in _MT_NAMES:
        return ContigClass.MT
    return ContigClass.OTHER


def is_sex_chrom(chrom: str) -> bool:
    """True if ``chrom`` is not an autosome (chr-prefix optional, any case).

    Kept as a name-only predicate for the identity-QC callers that need "not
    an autosome" and nothing finer. True for X, Y, M and MT, and also for any
    ``OTHER`` contig (alt, decoy, unplaced, random), since those are not
    autosomes either. Use ``classify_contig`` when PAR status matters.
    """
    return normalize_chrom(chrom) not in _AUTOSOME_NAMES
