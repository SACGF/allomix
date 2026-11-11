"""Contig classification for marker routing.

Every marker in allomix is keyed by ``(chrom, pos, ref, alt)`` and nothing else
interprets a coordinate, with one exception: deciding whether a chrX or chrY
site is pseudoautosomal. This module owns that decision and the broader
classification of a contig as autosome, sex chromosome, mitochondrion, or
"other" (alternate haplotypes, decoys, unplaced and random scaffolds, HLA
contigs), so the rest of the package can route markers without caring about
contig naming conventions or genome build.

Pseudoautosomal regions (PAR) come from ``bioutils.par`` (bioutils >= 0.6.2),
which packages the ``PAR`` rows of NCBI's ``assembly_regions.txt`` for each
build (GRCh37.p13, GRCh38.p14) in interbase coordinates: ``[start_i, end_i]``
covers the 1-based positions ``start_i + 1`` to ``end_i``, so a VCF ``POS`` is
inside when ``start_i < pos <= end_i``. The intervals, as bioutils gives them::

                GRCh37                    GRCh38
    X PAR1      60000-2699520             10000-2781479
    X PAR2      154931043-155260560       155701382-156030895
    Y PAR1      10000-2649520             10000-2781479
    Y PAR2      59034049-59363566         56887902-57217415

The mask is the union of both builds: a position is PAR if it is PAR in
either. The PAR1 intervals overlap; the PAR2 intervals do not overlap on
either chromosome (GRCh38 moved the chrX PAR2 by 770 kb and the chrY PAR2 by
2.1 Mb), so the positions between them are non-PAR.

Why the union, and why no genome build: allomix always excludes PAR markers,
so a site that is PAR in one build and non-PAR in the other can be excluded
in both. Outside the union the PAR status of a position is the same in both
builds, so the classification is a pure function of ``(chrom, pos)`` and the
tool does not need to know which build the VCF is on.

The cost is the positions that are PAR in one build only. In 1-based positions,
on chrX: 10001-60000 and 2699521-2781479 from PAR1 (50 kb and 82 kb), and the
GRCh37 PAR2 interval 154931044-155260560, which is non-PAR Xq28 in GRCh38
coordinates (330 kb). The GRCh38 chrX PAR2 interval lies beyond GRCh37's chrX
length, so it costs nothing there. On chrY: 2649521-2781479 from PAR1 (132 kb),
and the GRCh38 PAR2 interval 56887903-57217415, which is non-PAR distal Yq in
GRCh37 coordinates (330 kb); the GRCh37 chrY PAR2 interval lies beyond GRCh38's
chrY length. No sample-ID panel known to us targets any of this, and the
exclusions are counted and reported so a panel that does would show it. The
classification tests check both sides of every PAR edge, so a change to the
bioutils data fails them rather than moving the mask silently. The bioutils
data also carries T2T-CHM13v2.0; it is left out of the union because allomix
does not support T2T VCFs.
"""

from __future__ import annotations

from enum import Enum

from bioutils.par import get_par_map


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


#: Genome builds whose PAR intervals are masked.
PAR_BUILDS: tuple[str, ...] = ("GRCh37", "GRCh38")

#: Every PAR interval of ``PAR_BUILDS``, interbase ``(start_i, end_i)`` as
#: ``bioutils.par`` gives them, sorted, per normalised chromosome name ("X", "Y").
#: The builds' PAR1 intervals overlap; see the module docstring.
PAR_REGIONS: dict[str, tuple[tuple[int, int], ...]] = {
    chrom: tuple(
        sorted((s, e) for build in PAR_BUILDS for s, e in get_par_map(build)[chrom].values())
    )
    for chrom in ("X", "Y")
}

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
    return any(start_i < pos <= end_i for start_i, end_i in regions)


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
