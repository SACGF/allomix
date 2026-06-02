"""Per-site empirical sequencing error rate estimation.

Estimates per-marker, per-direction substitution rates from a training cohort
of joint-called VCFs. At hom-ref calls the observed ALT-read rate estimates
``P(observe ALT | true REF)`` (called ``e_refalt``); at hom-alt calls the
observed REF-read rate estimates ``P(observe REF | true ALT)`` (``e_altref``).

The two are not generally equal: oxidation damage, strand bias and flanking
context all produce direction-specific error rates. The output table is
consumed by ``chimerism.estimate_single_donor_bb`` and
``chimerism.estimate_multi_donor`` via the ``marker_errors`` parameter, and by
the host-presence detector (planned `src/allomix/detect.py`) as the per-site
background.

This module mirrors ``allomix.bias``: estimator pooled across reads, save/load
TSV with NA for missing per-direction entries, and a runtime loader that
applies a configurable floor so a zero observed rate cannot make a single stray
read produce infinite log-likelihood penalties.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from allomix.genotype import MarkerData, MarkerKey, marker_key

# Default floor applied at load time. Per-site rates with fewer observations
# than ``min_reads`` are excluded by the estimator, but even at the cutoff a
# site with zero observed errors would give ``e = 0`` and make any single stray
# read in a downstream sample produce ``-inf`` log-likelihood. The floor caps
# how confident we let any single site become; 1e-5 sits well below realistic
# clean per-site rates (~1e-4 to 5e-4) so it does not blunt the useful tail.
DEFAULT_ERROR_FLOOR = 1e-5


@dataclass
class MarkerError:
    """Per-marker, per-direction empirical error rates.

    ``e_refalt`` and ``e_altref`` may each be ``None`` independently if no
    confident observations were available in that direction. The runtime
    dispatch (in ``chimerism.total_log_likelihood_bb``) treats either ``None``
    as a fall-through to the global symmetric error model.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    e_refalt: float | None  # ALT-read rate at hom-ref calls
    e_altref: float | None  # REF-read rate at hom-alt calls
    n_reads_homref: int
    n_reads_homalt: int


def estimate_error_rates(
    marker_lists: list[list[MarkerData]],
    min_reads: int = 1000,
    max_vaf_homref: float = 0.10,
    min_vaf_homalt: float = 0.90,
) -> dict[MarkerKey, MarkerError]:
    """Estimate per-marker, per-direction error rates from training samples.

    Reads are pooled across the cohort, not averaged across samples, so a
    high-depth sample contributes proportionally more weight. This is the
    correct MLE when the per-direction error rate is shared across samples at
    a site.

    Args:
        marker_lists: List of MarkerData lists, one per training sample. Apply
            ``min_gq`` at parse time (e.g. ``parse_vcf(..., min_gq=20)``) to
            exclude low-confidence calls.
        min_reads: Minimum total reads required *per direction* to retain a
            site's estimate. Sites below threshold get ``None`` for that
            direction's rate (the runtime falls through to ``--error-rate``).
        max_vaf_homref: Drop hom-ref observations where ``ad_alt/dp`` exceeds
            this threshold. Protects against miscalled hets and contamination
            inflating the rate. Default 0.10, well above realistic error rates.
        min_vaf_homalt: Drop hom-alt observations where ``ad_alt/dp`` falls
            below this threshold (symmetric to ``max_vaf_homref``).

    Returns:
        Dict mapping (chrom, pos, ref, alt) to ``MarkerError``. Sites with no
        usable observations in either direction are omitted.
    """
    # Per-site accumulators
    n_alt_homref: dict[MarkerKey, int] = {}
    n_tot_homref: dict[MarkerKey, int] = {}
    n_ref_homalt: dict[MarkerKey, int] = {}
    n_tot_homalt: dict[MarkerKey, int] = {}
    info: dict[MarkerKey, tuple[str, int, str, str]] = {}

    for markers in marker_lists:
        for m in markers:
            dp = m.ad_ref + m.ad_alt
            if dp <= 0:
                continue
            key = marker_key(m)
            info.setdefault(key, (m.chrom, m.pos, m.ref, m.alt))
            vaf = m.ad_alt / dp
            if m.gt == (0, 0):
                if vaf > max_vaf_homref:
                    continue
                n_alt_homref[key] = n_alt_homref.get(key, 0) + m.ad_alt
                n_tot_homref[key] = n_tot_homref.get(key, 0) + dp
            elif m.gt == (1, 1):
                if vaf < min_vaf_homalt:
                    continue
                n_ref_homalt[key] = n_ref_homalt.get(key, 0) + m.ad_ref
                n_tot_homalt[key] = n_tot_homalt.get(key, 0) + dp
            # Hets are used by bias estimation, not error estimation.

    out: dict[MarkerKey, MarkerError] = {}
    for key, (chrom, pos, ref, alt) in info.items():
        nh_tot = n_tot_homref.get(key, 0)
        na_tot = n_tot_homalt.get(key, 0)
        e_ra: float | None = None
        e_ar: float | None = None
        if nh_tot >= min_reads:
            e_ra = n_alt_homref[key] / nh_tot
        if na_tot >= min_reads:
            e_ar = n_ref_homalt[key] / na_tot
        if e_ra is None and e_ar is None:
            continue
        out[key] = MarkerError(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            e_refalt=e_ra,
            e_altref=e_ar,
            n_reads_homref=nh_tot,
            n_reads_homalt=na_tot,
        )
    return out


def save_error_table(
    errors: dict[MarkerKey, MarkerError], path: Path | str,
) -> None:
    """Write error-rate estimates to a TSV file.

    Format:

        chrom  pos  ref  alt  e_refalt  e_altref  n_reads_homref  n_reads_homalt

    Missing per-direction rates are written as ``NA``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "chrom", "pos", "ref", "alt", "e_refalt", "e_altref",
                "n_reads_homref", "n_reads_homalt",
            ]
        )
        for key in sorted(errors.keys()):
            me = errors[key]
            w.writerow(
                [
                    me.chrom,
                    me.pos,
                    me.ref,
                    me.alt,
                    "NA" if me.e_refalt is None else f"{me.e_refalt:.6e}",
                    "NA" if me.e_altref is None else f"{me.e_altref:.6e}",
                    me.n_reads_homref,
                    me.n_reads_homalt,
                ]
            )


def load_error_table(
    path: Path | str,
    error_floor: float = DEFAULT_ERROR_FLOOR,
) -> dict[MarkerKey, tuple[float | None, float | None]]:
    """Load an error-rate table.

    Args:
        path: Path to a TSV produced by ``save_error_table``.
        error_floor: Lower bound applied to each non-``None`` per-direction
            rate. A site with zero observed errors would otherwise drive the
            log-likelihood to ``-inf`` on a single stray read; the floor keeps
            the likelihood finite while sitting well below realistic clean
            per-site rates. Set to 0 to disable.

    Returns:
        Dict mapping (chrom, pos, ref, alt) to ``(e_refalt, e_altref)``. Each
        entry of the tuple is ``None`` if the table stored ``NA`` in that
        column.
    """
    out: dict[MarkerKey, tuple[float | None, float | None]] = {}
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key: MarkerKey = (
                row["chrom"], int(row["pos"]), row["ref"], row["alt"],
            )
            e_ra = (
                None if row["e_refalt"] == "NA"
                else max(float(row["e_refalt"]), error_floor)
            )
            e_ar = (
                None if row["e_altref"] == "NA"
                else max(float(row["e_altref"]), error_floor)
            )
            out[key] = (e_ra, e_ar)
    return out


def errors_to_simple_dict(
    errors: dict[MarkerKey, MarkerError],
    error_floor: float = DEFAULT_ERROR_FLOOR,
) -> dict[MarkerKey, tuple[float | None, float | None]]:
    """Convert ``MarkerError`` dict to the ``(e_refalt, e_altref)`` form
    expected by the estimators, applying the same floor as the loader.
    """
    out: dict[MarkerKey, tuple[float | None, float | None]] = {}
    for key, me in errors.items():
        e_ra = None if me.e_refalt is None else max(me.e_refalt, error_floor)
        e_ar = None if me.e_altref is None else max(me.e_altref, error_floor)
        out[key] = (e_ra, e_ar)
    return out
