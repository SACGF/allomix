"""Relatedness estimation and sample-swap detection for QC.

Two allele-frequency-free identity checks computed from the same genotype and
allele-depth data the chimerism pipeline already parses:

1. ``relatedness_coefficient`` estimates the kinship between two reference
   samples (host and a donor, or donor and donor) with a somalier-style robust
   coefficient over shared autosomal markers. It is reported in output and,
   when a lab declares an expected relationship, compared against it to flag
   mislabelled or unexpectedly related reference samples.

2. ``admix_consistency`` checks the admixture sample against the supplied host
   and donor(s). At markers where host and all donors share the same homozygous
   genotype, the admixture must carry that homozygote up to sequencing error,
   regardless of mixing fraction. Excess minority-allele reads there mean the
   admixture carries a third genome (a sample swap or wrong-patient VCF).

Both are panel-agnostic and use only autosomal markers; sex-chromosome dosage is
unreliable in sex-mismatched transplants, consistent with the rest of the tool.

The robust coefficient is adapted from somalier (Pedersen & Quinlan, MIT
licensed). The math is standard kinship estimation.
"""

import math
from dataclasses import dataclass

from scipy.stats import binom, norm

from allomix.genotype import MarkerData, is_sex_chrom, marker_key

# ---------------------------------------------------------------------------
# Tunable constants
# ---------------------------------------------------------------------------

# Coefficient bands -> degree index. A coefficient strictly above a band's lower
# edge (and within the next one up) maps to that degree. Scale: ~1.0 identical,
# ~0.5 first-degree, ~0.25 second-degree, ~0.125 third-degree, ~0 unrelated.
IDENTICAL_MIN = 0.70
FIRST_DEGREE_MIN = 0.35
SECOND_DEGREE_MIN = 0.17
THIRD_DEGREE_MIN = 0.08

# Degree indices.
DEGREE_IDENTICAL = 0
DEGREE_FIRST = 1
DEGREE_SECOND = 2
DEGREE_THIRD = 3
DEGREE_UNRELATED = 4

DEGREE_LABELS = {
    DEGREE_IDENTICAL: "identical / duplicate",
    DEGREE_FIRST: "first-degree (parent/child/sibling)",
    DEGREE_SECOND: "second-degree (half-sib/uncle/grandparent)",
    DEGREE_THIRD: "third-degree (cousin)",
    DEGREE_UNRELATED: "unrelated",
}

# Declared expected-relationship strings accepted on input, mapped to a degree.
# "related" is a catch-all matching any related class (degrees 1-3); handled
# separately in the comparison since it is not a single degree. "identical" is
# deliberately NOT accepted as a declaration. The only case where host and donor
# are genuinely identical is an identical-twin (syngeneic) donor, and then there
# are no host/donor genetic differences to measure, so genotype-based chimerism
# does not apply and there is nothing useful to declare. Identical is still a
# *detected* outcome (flagged as a duplicate / unmeasurable pair).
DECLARED_DEGREE = {
    "first-degree": DEGREE_FIRST,
    "second-degree": DEGREE_SECOND,
    "third-degree": DEGREE_THIRD,
    "unrelated": DEGREE_UNRELATED,
}
RELATED_CATCH_ALL = "related"
#: Values accepted from the CLI / batch CSV. "NA" (or blank) means no expectation.
VALID_DECLARATIONS = (*DECLARED_DEGREE.keys(), RELATED_CATCH_ALL)

# When a declaration and the estimate sit on opposite sides of the
# related/unrelated boundary, only a *close* relationship (this degree or nearer:
# identical / first / second) makes it a hard FAIL. Third-degree (cousin) sits
# within sampling noise of the boundary, so such a crossing is a REVIEW, not a
# FAIL. This keeps the swap/mislabel signal (e.g. declared sibling, detected
# unrelated) without failing legitimate distant kin that estimate just over the
# line.
BOUNDARY_FAIL_MAX_DEGREE = DEGREE_SECOND

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

# Two-sided 95% normal critical value (z_0.975 ~= 1.9600) for the Wald CI.
_Z_95 = float(norm.ppf(0.975))


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------


@dataclass
class RelatednessResult:
    """Robust relatedness estimate between two samples over shared autosomes.

    Attributes:
        a_name: First sample label.
        b_name: Second sample label.
        coefficient: Robust relatedness coefficient, or None when there are too
            few shared heterozygous sites (< ``MIN_HET_SITES``) to estimate it.
        ci_low: Lower bound of the approximate 95% CI (None if no coefficient).
        ci_high: Upper bound of the approximate 95% CI (None if no coefficient).
        confidence: "low" / "med" / "high" from the usable het-site count.
        relationship: English relationship label for ``degree``.
        degree: Degree index 0-4 (see ``DEGREE_LABELS``), or None if no estimate.
        n_sites: Shared, clean, autosomal biallelic markers compared.
        het_a: Heterozygous sites in sample A.
        het_b: Heterozygous sites in sample B.
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
    degree: int | None
    n_sites: int
    het_a: int
    het_b: int
    shared_hets: int
    ibs0: int

    @property
    def pair(self) -> str:
        """Human-readable pair label, e.g. ``host vs donor1``."""
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
    """Outcome of comparing an estimated relatedness against a declared one.

    Attributes:
        pair: Pair label the verdict applies to.
        declared: The declared expected relationship.
        detected: The detected English relationship.
        status: "PASS", "REVIEW", or "FAIL".
        message: A one-line explanation suitable for a QC warning.
    """

    pair: str
    declared: str
    detected: str
    status: str
    message: str


# ---------------------------------------------------------------------------
# Relatedness coefficient
# ---------------------------------------------------------------------------


def _clean_dose(gt: tuple[int, int]) -> int | None:
    """Alt-allele dose {0,1,2} for a clean biallelic diploid GT, else None.

    Returns None for missing (-1) or multi-allelic (allele index > 1) calls so
    they are skipped rather than miscoded.
    """
    a, b = gt
    if a < 0 or b < 0 or a > 1 or b > 1:
        return None
    return a + b


def _degree_from_coef(coef: float) -> int:
    """Map a relatedness coefficient to a degree index (see DEGREE_LABELS)."""
    if coef > IDENTICAL_MIN:
        return DEGREE_IDENTICAL
    if coef > FIRST_DEGREE_MIN:
        return DEGREE_FIRST
    if coef > SECOND_DEGREE_MIN:
        return DEGREE_SECOND
    if coef > THIRD_DEGREE_MIN:
        return DEGREE_THIRD
    return DEGREE_UNRELATED


def _confidence(n_eff: int) -> str:
    """Categorical confidence from the usable het-site count."""
    if n_eff >= CONF_HIGH_HETS:
        return "high"
    if n_eff >= CONF_MED_HETS:
        return "med"
    return "low"


def _coef_ci(coef: float, n_eff: int) -> tuple[float, float]:
    """Approximate 95% CI for the coefficient.

    Treats the coefficient as a proportion estimated over ``n_eff`` =
    min(het_a, het_b) sites and forms a Wald interval. This is a deliberately
    simple analytic interval: it captures the dominant effect (the interval
    widens as the het-site count shrinks) without claiming exactness.
    """
    p = min(max(coef, 0.0), 1.0)
    se = math.sqrt(max(p * (1.0 - p), 1e-6) / n_eff)
    return coef - _Z_95 * se, coef + _Z_95 * se


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

    Args:
        a: First sample's markers.
        b: Second sample's markers.
        a_name: Label for the first sample.
        b_name: Label for the second sample.

    Returns:
        A ``RelatednessResult``. ``coefficient`` is None when fewer than
        ``MIN_HET_SITES`` heterozygous sites are shared, in which case ``degree``
        is None and ``relationship`` is "undetermined".
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
    degree = _degree_from_coef(coef)
    ci_low, ci_high = _coef_ci(coef, n_eff)
    return RelatednessResult(
        a_name=a_name,
        b_name=b_name,
        coefficient=coef,
        ci_low=ci_low,
        ci_high=ci_high,
        confidence=_confidence(n_eff),
        relationship=DEGREE_LABELS[degree],
        degree=degree,
        n_sites=n_sites,
        het_a=het_a,
        het_b=het_b,
        shared_hets=shared_hets,
        ibs0=ibs0,
    )


# ---------------------------------------------------------------------------
# Admixture consensus-homozygote consistency / swap check
# ---------------------------------------------------------------------------


def admix_consistency(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    error_rate: float = 0.01,
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

    Args:
        host: Parsed host markers.
        donors: One parsed marker list per donor.
        admix: Parsed admixture markers (raw allele depths).
        error_rate: Per-site probability of a minority-allele read under no
            swap. Using the full symmetric rate (rather than ``error_rate/3``)
            is the conservative choice: it raises the bar for calling a site
            discordant.
        min_dp: Minimum admixture depth for a marker to be tested.

    Returns:
        An ``AdmixConsistencyResult``.
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


# ---------------------------------------------------------------------------
# Expected-vs-detected comparison
# ---------------------------------------------------------------------------


#: Inputs that mean "no expectation declared": no verdict, and no error.
NO_EXPECTATION_VALUES = {"", "na"}


def _normalise_declaration(declared: str | None) -> str | None:
    """Lowercase a declaration; None for no-expectation; raise on a typo.

    Returns None for None / blank / "NA" (a deliberate "no expectation"), and the
    lowercased value when it is a recognised declaration. Raises ValueError for
    anything else, so a typo (e.g. "frist-degree") is a hard error rather than
    silently turning the relatedness check off.
    """
    if declared is None:
        return None
    d = declared.strip().lower()
    if d in NO_EXPECTATION_VALUES:
        return None
    if d not in VALID_DECLARATIONS:
        raise ValueError(
            f"unrecognised expected relatedness {declared!r}; choose from "
            f"{', '.join(VALID_DECLARATIONS)} or NA"
        )
    return d


def _verdict_status(declared: str, degree: int, tolerance: int) -> str:
    """Decide PASS / REVIEW / FAIL for a usable declaration and detected degree.

    The verdict is asymmetric, following the realistic failure modes:

    - Losing relatedness (a close relationship declared, but detected unrelated)
      is the random-swap signature and a FAIL. "Close" means second-degree or
      nearer; a third-degree (cousin) crossing is a REVIEW, since cousins
      routinely estimate just over the unrelated line at panel marker counts.
    - Gaining relatedness by accident is implausible except for sample reuse,
      which reads as identical/duplicate. A detected *identical* result is
      therefore a FAIL (reuse, or an identical-twin donor that cannot be
      monitored), but a moderate unexpected relationship (e.g. unrelated
      declared, first-degree detected) is only a REVIEW: it is not a swap
      signature and a related donor still yields a usable estimate.
    - Within the related class, a degree distance <= tolerance is PASS,
      otherwise REVIEW.
    """
    detected_identical = degree == DEGREE_IDENTICAL
    detected_unrelated = degree == DEGREE_UNRELATED

    # Detected identical: sample reuse, or an identical-twin (syngeneic) donor.
    # Either way the pair has no host/donor differences to measure, so FAIL.
    # ("identical" is not an accepted declaration, so this never matches one.)
    if detected_identical:
        return "FAIL"

    if declared == RELATED_CATCH_ALL:
        # Any related degree satisfies "related"; finding no relationship is
        # worth a look (REVIEW) but not a hard fail.
        return "REVIEW" if detected_unrelated else "PASS"

    if declared == "unrelated":
        # A non-identical unexpected relationship is not a random-swap signature.
        return "REVIEW" if not detected_unrelated else "PASS"

    decl_deg = DECLARED_DEGREE[declared]
    if detected_unrelated:
        return "FAIL" if decl_deg <= BOUNDARY_FAIL_MAX_DEGREE else "REVIEW"
    return "PASS" if abs(decl_deg - degree) <= tolerance else "REVIEW"


def evaluate_expected(
    result: RelatednessResult,
    declared: str | None,
    tolerance: int = 1,
) -> RelatednessVerdict | None:
    """Compare an estimated relatedness against a declared expectation.

    Returns None when there is no declaration (None / blank / "NA"), so callers
    that pass no expectation get no verdict. Raises ValueError on an
    unrecognised value (likely a typo), rather than silently skipping the check.

    Verdict rules (asymmetric, by realistic failure mode):
        - Losing relatedness: a close relationship declared (second-degree or
          nearer) but detected unrelated is a FAIL (the random-swap signature).
          A third-degree (cousin) or "related" catch-all crossing to unrelated
          is only REVIEW, because cousins routinely estimate just over the line.
        - Gaining relatedness: a detected identical/duplicate is a FAIL (sample
          reuse, or an identical-twin donor; either way unmeasurable). A moderate
          unexpected relationship (e.g. unrelated declared, first-degree
          detected) is only REVIEW.
        - Within the related class, degree distance ``d <= tolerance`` is PASS;
          a larger gap is REVIEW.
        - A declaration we cannot check (too few het sites, no coefficient) is
          REVIEW.

    Note: a duplicate (identical reference pair) is also flagged unconditionally
    in ``assess_quality``, independent of any declaration.

    Args:
        result: The estimated relatedness for the pair.
        declared: Declared expected relationship, or None / "NA" for none.
        tolerance: Allowed degree distance for a PASS (default 1).

    Returns:
        A ``RelatednessVerdict`` or None.

    Raises:
        ValueError: if ``declared`` is a non-blank, non-"NA" string that is not a
            recognised relationship.
    """
    declared_norm = _normalise_declaration(declared)
    if declared_norm is None:
        return None

    if result.degree is None or result.coefficient is None:
        return RelatednessVerdict(
            pair=result.pair,
            declared=declared_norm,
            detected="undetermined",
            status="REVIEW",
            message=(
                f"relatedness check inconclusive for {result.pair}: declared "
                f"{declared_norm}, too few shared het sites to estimate "
                f"(het_a={result.het_a}, het_b={result.het_b})"
            ),
        )

    status = _verdict_status(declared_norm, result.degree, tolerance)
    coef_str = f"r={result.coefficient:.2f}, {result.n_sites} markers"
    return RelatednessVerdict(
        pair=result.pair,
        declared=declared_norm,
        detected=result.relationship,
        status=status,
        message=(
            f"relatedness check {status} for {result.pair}: declared "
            f"{declared_norm}, detected {result.relationship} ({coef_str})"
        ),
    )


__all__ = [
    "RelatednessResult",
    "AdmixConsistencyResult",
    "RelatednessVerdict",
    "relatedness_coefficient",
    "admix_consistency",
    "evaluate_expected",
    "VALID_DECLARATIONS",
    "DEGREE_LABELS",
]
