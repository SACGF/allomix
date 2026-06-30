"""Relatedness estimation and sample-swap detection for QC.

Two allele-frequency-free identity checks from the genotype and allele-depth data
the chimerism pipeline already parses:

1. ``relatedness_coefficient`` estimates kinship between two reference samples
   (host/donor or donor/donor) with a somalier-style robust coefficient over
   shared autosomal markers. Reported in output and, when a lab declares an
   expected relationship, compared against it to flag mislabelled or unexpectedly
   related reference samples.

2. ``admix_consistency`` checks the admixture against host and donor(s). Where
   host and all donors share the same homozygous genotype, the admixture must
   carry that homozygote up to sequencing error, whatever the mixing fraction.
   Excess minority-allele reads there mean a third genome (sample swap or
   wrong-patient VCF).

Both use only autosomal markers; sex-chromosome dosage is unreliable in
sex-mismatched transplants. The robust coefficient is adapted from somalier
(Pedersen & Quinlan, MIT licensed); standard kinship estimation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy.stats import binom, norm

from allomix.constants import CI_LEVEL, DEFAULT_ERROR_RATE
from allomix.genotype import MarkerData, is_sex_chrom, marker_key


class Relatedness(str, Enum):
    """A host-vs-donor relatedness class, estimated from data or declared by a lab.

    Member value is the lowercase wire string used on the CLI / batch CSV and in messages, so a
    ``Relatedness`` is a ``str``. Members carry an ``ordinal`` on a 0 (identical) ..
    4 (unrelated) scale for degree-distance comparisons, the lower edge of the
    coefficient band that maps to them (``min_coef``), and an English ``label``.
    Coefficient scale: ~1.0 identical, ~0.5 first-degree, ~0.25 second-degree,
    ~0.125 third-degree, ~0 unrelated.

    Two members are context-specific:

    - ``IDENTICAL`` is only ever a *detected* outcome (a duplicate, or an
      identical-twin/syngeneic donor that cannot be monitored). It is not a valid
      declaration, so ``from_declared`` rejects it.
    - ``RELATED`` is only ever a *declared* catch-all ("any related degree"). It has
      no coefficient band, so ``from_coef`` never returns it.
    """

    # name = (wire value, ordinal, lower coefficient edge, label)
    IDENTICAL = ("identical", 0, 0.70, "identical / duplicate")
    FIRST_DEGREE = ("first-degree", 1, 0.35, "first-degree (parent/child/sibling)")
    SECOND_DEGREE = ("second-degree", 2, 0.17, "second-degree (half-sib/uncle/grandparent)")
    THIRD_DEGREE = ("third-degree", 3, 0.08, "third-degree (cousin)")
    UNRELATED = ("unrelated", 4, float("-inf"), "unrelated")
    RELATED = ("related", 5, float("nan"), "related (any degree)")

    def __new__(cls, value: str, ordinal: int, min_coef: float, label: str) -> Relatedness:
        member = str.__new__(cls, value)
        member._value_ = value
        member.ordinal = ordinal
        member.min_coef = min_coef
        member.label = label
        return member

    @classmethod
    def from_coef(cls, coef: float) -> Relatedness:
        """Map a relatedness coefficient to the highest band it clears.

        Always a concrete degree (identical..unrelated), never the ``RELATED``
        catch-all, whose ``min_coef`` is NaN so the ``>`` test is never true.
        """
        for member in cls:  # ordered identical..unrelated, descending edges
            if coef > member.min_coef:
                return member
        return cls.UNRELATED

    @classmethod
    def declared_values(cls) -> tuple[str, ...]:
        """Wire strings accepted as a declaration: every class except identical."""
        return tuple(m.value for m in cls if m is not cls.IDENTICAL)

    @classmethod
    def from_declared(cls, value: str | Relatedness | None) -> Relatedness | None:
        """Normalize a declared relatedness from the CLI / CSV to a member.

        Returns None for no expectation (None, blank, or "NA"); the matching member
        otherwise. Raises ValueError on an unrecognised value (likely a typo) or on
        "identical", a detected-only outcome that is never a valid declaration, so
        either is a hard error rather than silently turning the check off.
        """
        if value is None:
            return None
        if isinstance(value, cls):
            member = value
        else:
            text = value.strip().lower()
            if text in _NO_EXPECTATION_VALUES:
                return None
            try:
                member = cls(text)
            except ValueError:
                raise ValueError(
                    f"unrecognised expected relatedness {value!r}; choose from "
                    f"{', '.join(cls.declared_values())} or NA"
                ) from None
        if member is cls.IDENTICAL:
            raise ValueError(
                "'identical' is not a valid declared relatedness; it is a "
                "detected-only outcome (duplicate / identical-twin donor)"
            )
        return member


#: Inputs that mean "no expectation declared": no verdict, and no error.
_NO_EXPECTATION_VALUES = {"", "na"}

# When a declaration and the estimate sit on opposite sides of the
# related/unrelated boundary, only a *close* relationship (second-degree or nearer)
# is a hard FAIL. Third-degree (cousin) sits within sampling noise of the boundary,
# so such a crossing is REVIEW: keeps the swap/mislabel signal without failing
# legitimate distant kin that estimate just over the line.
BOUNDARY_FAIL_MAX_DEGREE = Relatedness.SECOND_DEGREE

# Minimum heterozygous sites (in the scarcer sample) before a coefficient is
# trusted at all; below this, ``coefficient`` is None.
MIN_HET_SITES = 5
# Categorical confidence from the usable het-site count.
CONF_HIGH_HETS = 40
CONF_MED_HETS = 20

# Admix-vs-(host+donor) consensus-homozygote swap check.
# Per-site significance for a single consensus-hom marker carrying excess
# minority-allele reads.
SITE_ALPHA = 1e-3
# Minimum consensus-homozygous markers before the overall swap p-value is
# meaningful; below this the check is reported but not acted on.
MIN_CONSENSUS = 20

# Two-sided normal critical value at CI_LEVEL (z_0.975 ~= 1.9600 for 95%) for the
# Wald CI: upper tail 1 - (1 - CI_LEVEL) / 2.
_Z_TWO_SIDED = float(norm.ppf(1.0 - (1.0 - CI_LEVEL) / 2.0))


@dataclass
class RelatednessResult:
    """Robust relatedness estimate between two samples over shared autosomes.

    Attributes:
        coefficient: Robust relatedness coefficient, or None when too few shared
            heterozygous sites (< ``MIN_HET_SITES``) to estimate it.
        ci_low: Lower bound of the approximate 95% CI (None if no coefficient).
        ci_high: Upper bound of the approximate 95% CI (None if no coefficient).
        confidence: "low" / "med" / "high" from the usable het-site count.
        relationship: English label for the detected relatedness (``degree.label``).
        degree: Detected ``Relatedness`` (always a concrete degree, never the
            ``RELATED`` catch-all), or None if no estimate.
        n_sites: Shared, clean, autosomal biallelic markers compared.
        shared_hets: Sites heterozygous in both samples.
        ibs0: Opposite-homozygote sites (0 vs 2 alt dose).
    """

    a_name: str
    b_name: str
    coefficient: float | None
    ci_low: float | None
    ci_high: float | None
    confidence: str
    relationship: str
    degree: Relatedness | None
    n_sites: int
    het_a: int
    het_b: int
    shared_hets: int
    ibs0: int

    @property
    def pair(self) -> str:
        return f"{self.a_name} vs {self.b_name}"


@dataclass
class AdmixConsistencyResult:
    """Admixture-vs-(host+donor) consensus-homozygote consistency check.

    Attributes:
        n_consensus_hom: Markers where host and all donors share one homozygous
            genotype and the admixture has usable depth.
        n_discordant: Consensus-hom markers where the admixture carries
            significantly more minority-allele reads than sequencing error
            alone explains.
        discordant_fraction: ``n_discordant / n_consensus_hom`` (0.0 if none).
        swap_pval: P(>= n_discordant discordant sites) under no swap, i.e. a
            Binomial(n_consensus_hom, SITE_ALPHA) tail. Small means the
            admixture carries a third genome.
    """

    n_consensus_hom: int
    n_discordant: int
    discordant_fraction: float
    swap_pval: float


@dataclass
class RelatednessVerdict:
    """Outcome of comparing an estimated relatedness against a declared one."""

    pair: str
    declared: str
    detected: str
    status: str  # "PASS", "REVIEW", or "FAIL"
    message: str  # one-line explanation suitable for a QC warning


def _clean_dose(gt: tuple[int, int]) -> int | None:
    """Alt-allele dose {0,1,2} for a clean biallelic diploid GT, else None.

    Returns None for missing (-1) or multi-allelic (allele index > 1) calls so
    they are skipped rather than miscoded.
    """
    a, b = gt
    if a < 0 or b < 0 or a > 1 or b > 1:
        return None
    return a + b


def _confidence(n_eff: int) -> str:
    """Categorical confidence from the usable het-site count."""
    if n_eff >= CONF_HIGH_HETS:
        return "high"
    if n_eff >= CONF_MED_HETS:
        return "med"
    return "low"


def _coef_ci(coef: float, n_eff: int) -> tuple[float, float]:
    """Approximate 95% CI for the coefficient.

    Treats the coefficient as a proportion over ``n_eff`` = min(het_a, het_b)
    sites and forms a Wald interval. Deliberately simple: captures the dominant
    effect (the interval widens as the het-site count shrinks) without claiming
    exactness.
    """
    p = min(max(coef, 0.0), 1.0)
    se = math.sqrt(max(p * (1.0 - p), 1e-6) / n_eff)
    return coef - _Z_TWO_SIDED * se, coef + _Z_TWO_SIDED * se


def relatedness_coefficient(
    a: list[MarkerData],
    b: list[MarkerData],
    a_name: str,
    b_name: str,
) -> RelatednessResult:
    """Estimate robust relatedness between two samples over shared autosomes.

    Uses the somalier-style coefficient::

        relatedness = (shared_hets - 2 * ibs0) / min(het_a, het_b)

    where ``shared_hets`` counts markers heterozygous in both samples and
    ``ibs0`` counts opposite homozygotes (0 vs 2 alt dose). Only shared,
    autosomal, clean biallelic genotypes are used. Allele frequencies are not
    needed, which suits an arbitrary panel.

    ``coefficient`` is None when fewer than ``MIN_HET_SITES`` heterozygous sites
    are shared, in which case ``degree`` is None and ``relationship`` is
    "undetermined".
    """
    b_by_key = {marker_key(m): m for m in b}

    het_a = het_b = shared_hets = ibs0 = n_sites = 0
    for ma in a:
        if is_sex_chrom(ma.chrom):
            continue
        mb = b_by_key.get(marker_key(ma))
        if mb is None:
            continue
        da = _clean_dose(ma.gt)
        db = _clean_dose(mb.gt)
        if da is None or db is None:
            continue
        n_sites += 1
        a_het = da == 1
        b_het = db == 1
        if a_het:
            het_a += 1
        if b_het:
            het_b += 1
        if a_het and b_het:
            shared_hets += 1
        elif (da == 0 and db == 2) or (da == 2 and db == 0):
            ibs0 += 1

    n_eff = min(het_a, het_b)
    if n_eff < MIN_HET_SITES:
        return RelatednessResult(
            a_name=a_name,
            b_name=b_name,
            coefficient=None,
            ci_low=None,
            ci_high=None,
            confidence=_confidence(n_eff),
            relationship="undetermined",
            degree=None,
            n_sites=n_sites,
            het_a=het_a,
            het_b=het_b,
            shared_hets=shared_hets,
            ibs0=ibs0,
        )

    coef = (shared_hets - 2 * ibs0) / n_eff
    degree = Relatedness.from_coef(coef)
    ci_low, ci_high = _coef_ci(coef, n_eff)
    return RelatednessResult(
        a_name=a_name,
        b_name=b_name,
        coefficient=coef,
        ci_low=ci_low,
        ci_high=ci_high,
        confidence=_confidence(n_eff),
        relationship=degree.label,
        degree=degree,
        n_sites=n_sites,
        het_a=het_a,
        het_b=het_b,
        shared_hets=shared_hets,
        ibs0=ibs0,
    )


def admix_consistency(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    error_rate: float = DEFAULT_ERROR_RATE,
    min_dp: int = 1,
) -> AdmixConsistencyResult:
    """Check the admixture against host+donor at consensus-homozygous markers.

    At markers where host and every donor share the same homozygous genotype,
    the admixture is a mixture of identical homozygotes, so it must show that
    homozygote up to sequencing error, whatever the mixing fraction. Each such
    marker is tested with a binomial tail on the minority-allele read count
    against ``error_rate``; a marker is discordant when that p-value is below
    ``SITE_ALPHA``. The overall ``swap_pval`` is the Binomial tail for seeing at
    least ``n_discordant`` discordant markers by chance under no swap.

    This is complementary to the MLE goodness-of-fit, which only uses
    informative markers and never tests these consensus sites.

    ``error_rate`` is the per-site probability of a minority-allele read under no
    swap. Using the full symmetric rate (not ``error_rate/3``) is the conservative
    choice: it raises the bar for calling a site discordant.
    """
    donor_maps = [{marker_key(m): m for m in d} for d in donors]
    admix_map = {marker_key(m): m for m in admix}

    n_consensus = 0
    n_discordant = 0
    for mh in host:
        if is_sex_chrom(mh.chrom):
            continue
        dose_h = _clean_dose(mh.gt)
        if dose_h is None or dose_h == 1:
            continue  # host must be homozygous
        key = marker_key(mh)
        # All donors must share the same homozygous genotype.
        consensus = True
        for dm in donor_maps:
            md = dm.get(key)
            if md is None or _clean_dose(md.gt) != dose_h:
                consensus = False
                break
        if not consensus:
            continue
        ma = admix_map.get(key)
        if ma is None or ma.dp < min_dp or ma.dp <= 0:
            continue
        # Minority allele is the one absent from the consensus homozygote.
        minor_reads = ma.ad_alt if dose_h == 0 else ma.ad_ref
        n_consensus += 1
        # P(>= minor_reads) under Binomial(dp, error_rate).
        p_site = float(binom.sf(minor_reads - 1, ma.dp, error_rate))
        if p_site < SITE_ALPHA:
            n_discordant += 1

    if n_consensus == 0:
        return AdmixConsistencyResult(0, 0, 0.0, 1.0)

    swap_pval = float(binom.sf(n_discordant - 1, n_consensus, SITE_ALPHA))
    return AdmixConsistencyResult(
        n_consensus_hom=n_consensus,
        n_discordant=n_discordant,
        discordant_fraction=n_discordant / n_consensus,
        swap_pval=swap_pval,
    )


def _verdict_status(
    declared_relatedness: Relatedness, detected_relatedness: Relatedness, tolerance: int
) -> str:
    """Decide PASS / REVIEW / FAIL for a usable declaration and detected estimate.

    Asymmetric, following the realistic failure modes:

    - Losing relatedness (close relationship declared, detected unrelated) is the
      random-swap signature and a FAIL. "Close" means second-degree or nearer; a
      third-degree (cousin) crossing is REVIEW, since cousins routinely estimate
      just over the unrelated line at panel marker counts.
    - Gaining relatedness by accident is implausible except for sample reuse, which
      reads as identical/duplicate. A detected *identical* is therefore a FAIL
      (reuse, or an identical-twin donor that cannot be monitored); a moderate
      unexpected relationship (e.g. unrelated declared, first-degree detected) is
      only REVIEW (not a swap signature, and a related donor still yields a usable
      estimate).
    - Within the related class, degree distance <= tolerance is PASS, else REVIEW.

    ``declared_relatedness`` and ``detected_relatedness`` share one ``Relatedness``
    scale, so the degree distance is just their ``ordinal`` difference.
    """
    # Detected identical: sample reuse or an identical-twin (syngeneic) donor;
    # either way no host/donor differences to measure, so FAIL. (IDENTICAL is not an
    # accepted declaration, so this never matches the declared side.)
    if detected_relatedness == Relatedness.IDENTICAL:
        return "FAIL"

    if declared_relatedness == Relatedness.RELATED:
        # Any related degree satisfies "related"; no relationship is REVIEW, not FAIL.
        return "REVIEW" if detected_relatedness == Relatedness.UNRELATED else "PASS"

    if declared_relatedness == Relatedness.UNRELATED:
        # A non-identical unexpected relationship is not a random-swap signature.
        return "PASS" if detected_relatedness == Relatedness.UNRELATED else "REVIEW"

    # Declared a specific degree (first/second/third).
    if detected_relatedness == Relatedness.UNRELATED:
        too_close = declared_relatedness.ordinal <= BOUNDARY_FAIL_MAX_DEGREE.ordinal
        return "FAIL" if too_close else "REVIEW"
    distance = abs(declared_relatedness.ordinal - detected_relatedness.ordinal)
    return "PASS" if distance <= tolerance else "REVIEW"


def evaluate_expected(
    result: RelatednessResult,
    declared: str | Relatedness | None,
    tolerance: int = 1,
) -> RelatednessVerdict | None:
    """Compare an estimated relatedness against a declared expectation.

    Returns None when there is no declaration (None / blank / "NA"). Raises
    ValueError on an unrecognised value (likely a typo) rather than silently
    skipping the check.

    Verdict rules (asymmetric, by realistic failure mode):
        - Losing relatedness: a close relationship declared (second-degree or
          nearer) but detected unrelated is a FAIL (random-swap signature). A
          third-degree (cousin) or "related" catch-all crossing to unrelated is
          only REVIEW, because cousins routinely estimate just over the line.
        - Gaining relatedness: a detected identical/duplicate is a FAIL (sample
          reuse, or an identical-twin donor; either way unmeasurable). A moderate
          unexpected relationship (e.g. unrelated declared, first-degree detected)
          is only REVIEW.
        - Within the related class, degree distance ``d <= tolerance`` is PASS; a
          larger gap is REVIEW.
        - A declaration we cannot check (too few het sites, no coefficient) is
          REVIEW.

    Note: a duplicate (identical reference pair) is also flagged unconditionally in
    ``assess_quality``, independent of any declaration.

    ``tolerance`` is the allowed degree distance for a PASS (default 1).

    Raises:
        ValueError: if ``declared`` is a non-blank, non-"NA" string that is not a
            recognised relationship.
    """
    declared_relatedness = Relatedness.from_declared(declared)
    if declared_relatedness is None:
        return None

    if result.degree is None or result.coefficient is None:
        return RelatednessVerdict(
            pair=result.pair,
            declared=declared_relatedness.value,
            detected="undetermined",
            status="REVIEW",
            message=(
                f"relatedness check inconclusive for {result.pair}: declared "
                f"{declared_relatedness.value}, too few shared het sites to estimate "
                f"(het_a={result.het_a}, het_b={result.het_b})"
            ),
        )

    status = _verdict_status(declared_relatedness, result.degree, tolerance)
    coef_str = f"r={result.coefficient:.2f}, {result.n_sites} markers"
    return RelatednessVerdict(
        pair=result.pair,
        declared=declared_relatedness.value,
        detected=result.relationship,
        status=status,
        message=(
            f"relatedness check {status} for {result.pair}: declared "
            f"{declared_relatedness.value}, detected {result.relationship} ({coef_str})"
        ),
    )


__all__ = [
    "RelatednessResult",
    "AdmixConsistencyResult",
    "RelatednessVerdict",
    "Relatedness",
    "relatedness_coefficient",
    "admix_consistency",
    "evaluate_expected",
]
