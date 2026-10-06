"""Sex inference for reference samples and host/donor pair status.

Infers the sex of each reference sample (host, donors) from its own genotypes,
compares it with any declared sex, and summarises the pair as matched or
mismatched so later phases can route non-PAR chrX markers on it (#50, #46,
#48). Sex is never inferred for the admixture sample, which is a mixture.

Primary signal: heterozygosity at non-PAR chrX sites. A female carries two X
chromosomes and is heterozygous at about the panel's design het rate (taken
from the same sample's autosomes); a male is hemizygous, so a diploid caller
emits a het there only by error. The het count ``k`` among ``n`` usable sites
is compared by binomial likelihood ratio between those two models and called
when ``|log10 LR|`` clears ``SEX_LR_THRESHOLD``. Below ``MIN_X_SITES`` usable
sites the result is ``UNAVAILABLE`` (no usable chrX content on this panel);
with sites but no clear call it is ``AMBIGUOUS``.

Secondary signal: chrY relative depth, the lab statistic ``max(non-PAR chrY
DP) / median(autosomal DP)``, when a forced midpoint pileup of the reference
sample is supplied (``--ref-depth-vcf``; an invariant chrY site has no GATK
record, so the genotype VCF cannot provide it). It is graded with the lab's
thresholds (``chry_call``: ``M`` / ``M*`` / ``F*`` / ``F``) and used only to
resolve a chrX call that is ``AMBIGUOUS`` or ``UNAVAILABLE``, never to overrule
a confident one. Only an unstarred ``M`` resolves: chrY depth is positive
evidence for a male, but its absence is not evidence for a female until it is
known that the panel captures the chrY targets at all (a panel whose intervals
list chrY sites the capture does not pull down would read every sample as
``F``).

Declared sex (``--recipient-sex`` / ``--donor-sex``) is used two ways: it
resolves an ambiguous or unavailable inference for routing, and a confident
inference that contradicts it is a sample mix-up signal that
``qc.assess_quality`` turns into a FAIL.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import binom

from allomix.constants import (
    CHRY_REL_DEPTH_FEMALE_WEAK,
    CHRY_REL_DEPTH_MALE,
    CHRY_REL_DEPTH_MALE_WEAK,
    FEMALE_X_HET_SCALE,
    MALE_X_SPURIOUS_HET,
    MIN_X_SITES,
    SEX_LR_THRESHOLD,
)
from allomix.contigs import ContigClass, classify_contig
from allomix.genotype import MarkerData
from allomix.sex_types import PairStatus, Sex

#: One forced-pileup record: ``(chrom, pos, dp)`` for a single sample.
RegionDepth = tuple[str, int, int]

# Clamp on the female-model het probability so a sample with an extreme
# autosomal het rate (or none) still gives a finite, sensible likelihood.
_P_FEMALE_MIN = 0.05
_P_FEMALE_MAX = 0.95
# Female het probability when the sample has no usable autosomal sites to
# estimate the panel's design het rate from.
_P_FEMALE_DEFAULT = 0.5


# ``Sex`` and ``PairStatus`` are defined in ``allomix.sex_types`` (a leaf
# module) so ``genotype.classify_markers`` can route chrX markers on the pair
# status without importing this module; re-exported here for callers.


@dataclass
class SexInference:
    """Sex inference for one reference sample.

    Attributes:
        sex: Inferred sex from chrX heterozygosity (``AMBIGUOUS`` /
            ``UNAVAILABLE`` when it did not resolve).
        n_x_sites: Usable non-PAR chrX sites (called, passing depth and GQ).
        n_x_het: Heterozygous calls among them.
        x_het_rate: ``n_x_het / n_x_sites``, or None when there are no sites.
        log10_lr: log10 likelihood ratio, female model over male model, or None
            when below ``MIN_X_SITES``. Positive favours female.
        chry_rel_depth: ``max(non-PAR chrY DP) / median(autosomal DP)`` from the
            forced-pileup depth input, or None when none was given or it had no
            chrY or autosomal sites.
        n_y_sites: Usable non-PAR chrY sites seen in the genotype VCF
            (informational).
        declared: Declared sex, or None.
        source: How ``effective`` was arrived at: ``"inferred"`` (no
            declaration, or an unresolved inference with none), ``"declared"``
            (unresolved inference resolved by the declaration),
            ``"inferred+declared"`` (both agree), ``"conflict"`` (confident
            inference contradicts the declaration) or ``"unavailable"`` (no
            chrX content and nothing declared).
        chry_call: The lab's grade of ``chry_rel_depth`` (``"M"``, ``"M*"``,
            ``"F*"``, ``"F"``), or None without it.
        chry_resolved: True when ``sex`` came from the chrY depth because the
            chrX call was ambiguous or unavailable.
    """

    sex: Sex
    n_x_sites: int
    n_x_het: int
    x_het_rate: float | None
    log10_lr: float | None
    chry_rel_depth: float | None
    n_y_sites: int
    declared: Sex | None
    source: str
    chry_call: str | None = None
    chry_resolved: bool = False

    @property
    def effective(self) -> Sex:
        """Sex to use for routing.

        The inference when confident (including a conflict: the data win over
        the label, and QC fails the sample); the declaration when the inference
        is ambiguous or unavailable and a sex was declared; otherwise the
        unresolved inferred state.
        """
        if self.sex.confident:
            return self.sex
        if self.declared is not None and self.declared.confident:
            return self.declared
        return self.sex


@dataclass
class SexResult:
    """Sex inference for every reference sample plus the pair status."""

    host: SexInference
    donors: list[SexInference]
    pair: PairStatus


def parse_declared_sex(text: str | None) -> Sex | None:
    """Parse a declared sex from the CLI.

    ``f`` / ``female`` -> ``FEMALE``, ``m`` / ``male`` -> ``MALE``,
    case-insensitive and whitespace-stripped. Anything else (including None and
    blank) is None: the text is kept for display only and ignored by the logic.
    """
    if text is None:
        return None
    t = text.strip().lower()
    if t in ("f", "female"):
        return Sex.FEMALE
    if t in ("m", "male"):
        return Sex.MALE
    return None


def split_region_depths(depths: list[RegionDepth]) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(chrY non-PAR depths, autosomal depths)`` as arrays."""
    y = [dp for chrom, pos, dp in depths if classify_contig(chrom, pos) is ContigClass.Y_NONPAR]
    auto = [dp for chrom, pos, dp in depths if classify_contig(chrom, pos) is ContigClass.AUTOSOME]
    return np.array(y, dtype=float), np.array(auto, dtype=float)


def chry_depth_ratio(depths: list[RegionDepth]) -> tuple[float | None, int]:
    """The lab statistic ``max(non-PAR chrY DP) / median(autosomal DP)`` for one sample.

    Returns:
        ``(ratio, n_chry_sites)``; the ratio is None when there are no chrY
        sites, no autosomal sites, or a zero autosomal median.
    """
    y, auto = split_region_depths(depths)
    if y.size == 0 or auto.size == 0:
        return None, int(y.size)
    med = float(np.median(auto))
    if med <= 0.0:
        return None, int(y.size)
    return float(y.max() / med), int(y.size)


def grade_chry_depth(rel_depth: float | None) -> str | None:
    """Grade a chrY relative depth with the lab thresholds.

    ``M`` above ``CHRY_REL_DEPTH_MALE``, ``M*`` above
    ``CHRY_REL_DEPTH_MALE_WEAK``, ``F*`` above ``CHRY_REL_DEPTH_FEMALE_WEAK``,
    else ``F``; None when there is no depth.
    """
    if rel_depth is None:
        return None
    if rel_depth > CHRY_REL_DEPTH_MALE:
        return "M"
    if rel_depth > CHRY_REL_DEPTH_MALE_WEAK:
        return "M*"
    if rel_depth > CHRY_REL_DEPTH_FEMALE_WEAK:
        return "F*"
    return "F"


def _usable(m: MarkerData, min_dp: int, min_gq: int) -> bool:
    """Called, clean biallelic diploid GT passing the depth and GQ filters."""
    a, b = m.gt
    if a < 0 or b < 0 or a > 1 or b > 1:
        return False
    if m.dp < min_dp:
        return False
    return m.gq is None or m.gq >= min_gq


def _source(inferred: Sex, declared: Sex | None) -> str:
    if declared is None or not declared.confident:
        return "unavailable" if inferred is Sex.UNAVAILABLE else "inferred"
    if not inferred.confident:
        return "declared"
    return "inferred+declared" if inferred is declared else "conflict"


def infer_sex(
    markers: list[MarkerData],
    *,
    min_dp: int,
    min_gq: int,
    declared: Sex | None = None,
    region_depth: list[RegionDepth] | None = None,
) -> SexInference:
    """Infer the sex of one reference sample from its chrX heterozygosity.

    Only non-PAR chrX sites with a called diploid genotype, ``dp >= min_dp`` and
    ``gq >= min_gq`` (when GQ is present) are used. The female model's het
    probability is ``FEMALE_X_HET_SCALE`` times the autosomal het rate of the
    same sample under the same filters (clamped to [0.05, 0.95]; 0.5 when the
    sample has no usable autosomal sites); the male model's is
    ``MALE_X_SPURIOUS_HET``. The binomial log10 likelihood ratio of the chrX het
    count under the two models calls FEMALE at or above ``SEX_LR_THRESHOLD``,
    MALE at or below its negative, AMBIGUOUS between. Fewer than ``MIN_X_SITES``
    usable chrX sites gives UNAVAILABLE.

    When ``region_depth`` is given, its chrY relative depth is graded with the
    lab thresholds, and an unstarred ``M`` grade turns an AMBIGUOUS or
    UNAVAILABLE chrX call into MALE (``chry_resolved``). A confident chrX call is
    never changed.

    Args:
        markers: Parsed reference-sample markers (host or one donor).
        min_dp: Minimum depth at a site for it to count.
        min_gq: Minimum GQ at a site for it to count (ignored when GQ absent).
        declared: Declared sex, or None.
        region_depth: Forced-pileup ``(chrom, pos, dp)`` depths for this sample
            (``sex_mismatch.region_depths_from_vcf``), or None.

    Returns:
        A ``SexInference`` carrying the call, the counts behind it, and the
        declared sex with how the two were reconciled.
    """
    n_x = n_x_het = n_auto = n_auto_het = n_y = 0
    for m in markers:
        if not _usable(m, min_dp, min_gq):
            continue
        het = m.gt[0] != m.gt[1]
        cls = classify_contig(m.chrom, m.pos)
        if cls is ContigClass.AUTOSOME:
            n_auto += 1
            n_auto_het += het
        elif cls is ContigClass.X_NONPAR:
            n_x += 1
            n_x_het += het
        elif cls is ContigClass.Y_NONPAR:
            n_y += 1

    x_het_rate = n_x_het / n_x if n_x else None
    if n_x < MIN_X_SITES:
        inferred = Sex.UNAVAILABLE
        log10_lr: float | None = None
    else:
        p_f = FEMALE_X_HET_SCALE * (n_auto_het / n_auto) if n_auto else _P_FEMALE_DEFAULT
        p_f = min(max(p_f, _P_FEMALE_MIN), _P_FEMALE_MAX)
        p_m = MALE_X_SPURIOUS_HET
        ll_f = float(binom.logpmf(n_x_het, n_x, p_f))
        ll_m = float(binom.logpmf(n_x_het, n_x, p_m))
        log10_lr = (ll_f - ll_m) / math.log(10.0)
        if log10_lr >= SEX_LR_THRESHOLD:
            inferred = Sex.FEMALE
        elif log10_lr <= -SEX_LR_THRESHOLD:
            inferred = Sex.MALE
        else:
            inferred = Sex.AMBIGUOUS

    rel_depth = chry_depth_ratio(region_depth)[0] if region_depth is not None else None
    y_call = grade_chry_depth(rel_depth)
    chry_resolved = not inferred.confident and y_call == "M"
    if chry_resolved:
        inferred = Sex.MALE

    return SexInference(
        sex=inferred,
        n_x_sites=n_x,
        n_x_het=n_x_het,
        x_het_rate=x_het_rate,
        log10_lr=log10_lr,
        chry_rel_depth=rel_depth,
        n_y_sites=n_y,
        declared=declared,
        source=_source(inferred, declared),
        chry_call=y_call,
        chry_resolved=chry_resolved,
    )


def pair_status(host: Sex, donors: list[Sex]) -> PairStatus:
    """Classify the host/donor(s) sex relationship from effective sexes.

    ``UNKNOWN`` if any party is ambiguous or unavailable, ``MATCHED_FEMALE`` /
    ``MATCHED_MALE`` when every party has the same sex, ``MISMATCHED``
    otherwise.
    """
    sexes = [host, *donors]
    if any(not s.confident for s in sexes):
        return PairStatus.UNKNOWN
    if all(s is Sex.FEMALE for s in sexes):
        return PairStatus.MATCHED_FEMALE
    if all(s is Sex.MALE for s in sexes):
        return PairStatus.MATCHED_MALE
    return PairStatus.MISMATCHED


def assess_sex(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    *,
    min_dp: int,
    min_gq: int,
    declared_host: Sex | None = None,
    declared_donors: list[Sex | None] | None = None,
    host_region_depth: list[RegionDepth] | None = None,
    donor_region_depths: list[list[RegionDepth] | None] | None = None,
) -> SexResult:
    """Infer sex for the host and every donor and summarise the pair.

    ``declared_donors`` and ``donor_region_depths`` are aligned with ``donors``
    (one entry per donor, None for none); None means nothing for any donor.
    """
    decl = declared_donors or [None] * len(donors)
    depths = donor_region_depths or [None] * len(donors)
    host_inf = infer_sex(
        host,
        min_dp=min_dp,
        min_gq=min_gq,
        declared=declared_host,
        region_depth=host_region_depth,
    )
    donor_infs = [
        infer_sex(d, min_dp=min_dp, min_gq=min_gq, declared=decl[i], region_depth=depths[i])
        for i, d in enumerate(donors)
    ]
    pair = pair_status(host_inf.effective, [d.effective for d in donor_infs])
    return SexResult(host=host_inf, donors=donor_infs, pair=pair)


__all__ = [
    "PairStatus",
    "RegionDepth",
    "Sex",
    "SexInference",
    "SexResult",
    "assess_sex",
    "chry_depth_ratio",
    "grade_chry_depth",
    "infer_sex",
    "pair_status",
    "parse_declared_sex",
    "split_region_depths",
]
