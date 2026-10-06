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

chrY relative depth is reserved as a secondary signal (``chry_rel_depth``,
``n_y_sites``) and is not yet computed: an invariant chrY site has no GATK
record, so the depth has to come from a pileup the pipeline does not yet make
(plan Phase 3).

Declared sex (``--recipient-sex`` / ``--donor-sex``) is used two ways: it
resolves an ambiguous or unavailable inference for routing, and a confident
inference that contradicts it is a sample mix-up signal that
``qc.assess_quality`` turns into a FAIL.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.stats import binom

from allomix.constants import (
    FEMALE_X_HET_SCALE,
    MALE_X_SPURIOUS_HET,
    MIN_X_SITES,
    SEX_LR_THRESHOLD,
)
from allomix.contigs import ContigClass, classify_contig
from allomix.genotype import MarkerData
from allomix.sex_types import PairStatus, Sex

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
        chry_rel_depth: Reserved for the chrY depth secondary; always None in
            this release.
        n_y_sites: Usable non-PAR chrY sites seen (informational).
        declared: Declared sex, or None.
        source: How ``effective`` was arrived at: ``"inferred"`` (no
            declaration, or an unresolved inference with none), ``"declared"``
            (unresolved inference resolved by the declaration),
            ``"inferred+declared"`` (both agree), ``"conflict"`` (confident
            inference contradicts the declaration) or ``"unavailable"`` (no
            chrX content and nothing declared).
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

    Args:
        markers: Parsed reference-sample markers (host or one donor).
        min_dp: Minimum depth at a site for it to count.
        min_gq: Minimum GQ at a site for it to count (ignored when GQ absent).
        declared: Declared sex, or None.

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

    return SexInference(
        sex=inferred,
        n_x_sites=n_x,
        n_x_het=n_x_het,
        x_het_rate=x_het_rate,
        log10_lr=log10_lr,
        chry_rel_depth=None,
        n_y_sites=n_y,
        declared=declared,
        source=_source(inferred, declared),
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
) -> SexResult:
    """Infer sex for the host and every donor and summarise the pair.

    ``declared_donors`` is aligned with ``donors`` (one entry per donor, None for
    no declaration); None means nothing declared for any donor.
    """
    decl = declared_donors or [None] * len(donors)
    host_inf = infer_sex(host, min_dp=min_dp, min_gq=min_gq, declared=declared_host)
    donor_infs = [
        infer_sex(d, min_dp=min_dp, min_gq=min_gq, declared=decl[i]) for i, d in enumerate(donors)
    ]
    pair = pair_status(host_inf.effective, [d.effective for d in donor_infs])
    return SexResult(host=host_inf, donors=donor_infs, pair=pair)


__all__ = [
    "PairStatus",
    "Sex",
    "SexInference",
    "SexResult",
    "assess_sex",
    "infer_sex",
    "pair_status",
    "parse_declared_sex",
]
