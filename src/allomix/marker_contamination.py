"""Per-marker co-pooled contamination correction (Step 30, issue #30).

On co-pooled panels (index hopping on a patterned flowcell) a donor-homozygous
marker carries extra reads on the host (donor-absent) allele from co-pooled
genomes that happen to carry that allele. The per-sample contamination scalar
(``allomix.contamination``) does not localize this: it taxes every donor-absent
marker by the average, while the inflation actually lands only on the markers a
co-pooled genome carries, scaling with how many co-pooled individuals carry the
host allele there. Left uncorrected this puts a positive floor under the MLE at
true-zero host and inflates every low-dilution estimate upward.

The correction predicts the per-marker contamination on donor-hom markers from
the co-pooled carrier dose and subtracts it from the host-allele count before the
MLE. The host signal is identical at every donor-hom marker (independent of
carrier count); only the contamination scales with it. So per marker:

    host_allele_reads_corrected = host_allele_reads - slope * n_carriers * depth

where ``slope`` is the per-carrier contamination rate and ``n_carriers`` is the
capped number of co-pooled individuals carrying the host allele at that site.
Only the dose-dependent part is subtracted; the flat error floor stays with the
per-site error model (``allomix.error_rates``, issue #23), so it is not
double-counted.

Two pieces are measured per run, not hardcoded (see ``claude/step30_design.md``):

- GATE (per flowcell): at consensus-homozygous sites (host and every donor
  homozygous for the same allele, so the minor allele is pure background), fit
  the minor-allele fraction against the co-pooled carrier dose. A significant
  positive slope means this flowcell has real dose-dependent contamination. A
  clean flowcell has a flat (non-significant) slope and self-selects out: the
  correction becomes a no-op. The consensus-hom slope predicts the
  informative-marker slope at r=0.92, so the gate transfers reliably.
- MAGNITUDE (per patient): the consensus-hom slope does NOT transfer 1:1 to the
  informative markers (it over-predicts on clean mixtures, under-predicts on
  dirty ones). So once gated in, the correction slope is calibrated on the
  informative donor-hom markers themselves (host-allele fraction against carrier
  dose, weighted, pooled across a patient's serial timepoints to beat per-sample
  noise).

This module mirrors ``allomix.bias`` / ``allomix.error_rates``: an estimator that
builds a table from a training cohort, save/load of that table, and a lightweight
runtime object (``ContaminationCorrection``) carried on a
``allomix.likelihood.PanelCalibration`` and consumed by
``allomix.chimerism.estimate_single_donor_bb``. Only donor-homozygous markers
(Vynck types 0 and 1) are corrected; donor-het markers are left alone (smaller,
near-balanced effect).
"""

import math
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from allomix.constants import DEFAULT_MIN_DP, DEFAULT_MIN_GQ
from allomix.genotype import (
    InformativeMarker,
    MarkerData,
    MarkerKey,
    classify_markers,
    marker_key,
)

# Cap on the carrier dose. The host allele is common, so the dose-response is
# measurable but saturates; a hard cap stops a high-frequency site dominating the
# regression. Validated at 5 on SRP434573 (carrier COUNT, not allele copies; see
# claude/step30_design.md open question 1).
DEFAULT_DOSE_CAP = 5

# Gate defaults. The correction applies only when the consensus-hom slope is
# significantly positive (alpha) AND exceeds ``min_slope`` (minimum per-carrier
# effect worth correcting, in fraction-of-depth units). min_slope defaults to 0,
# so the gate is significance-only unless a caller sets it.
DEFAULT_GATE_ALPHA = 0.05
DEFAULT_GATE_MIN_SLOPE = 0.0


@dataclass(frozen=True)
class ContaminationCorrection:
    """Runtime per-marker contamination correction consumed during estimation.

    Carried on ``allomix.likelihood.PanelCalibration.contamination_correction``
    and applied by ``apply_contamination_correction`` before the single-donor
    MLE. When ``gated`` is False or ``slope`` is non-positive the correction is a
    no-op, so loading a table built on a clean flowcell leaves the estimate
    unchanged.

    Attributes:
        carriers: Per-marker raw cohort carrier counts ``(n_ref, n_alt)``: cohort
            individuals carrying the REF and ALT allele. The host's and donor's own
            contributions are removed at correction time, so the cohort is expected
            to include the host and donor (who were on the flowcell).
        slope: Per-carrier contamination rate (fraction of depth), subtracted as
            ``slope * dose * depth`` from the host-allele count. 0 on a gated-out
            (clean) flowcell.
        dose_cap: Maximum carrier dose used (see ``DEFAULT_DOSE_CAP``).
        gated: Whether the flowcell's consensus-hom dose response was significant.
            False leaves the estimate byte-identical to no correction.
        gate_slope: Consensus-hom slope the gate decision was made on (provenance).
        gate_p_value: One-sided p-value of ``gate_slope`` (provenance).
        n_consensus: Consensus-hom observations behind the gate (provenance).
        n_informative: Informative donor-hom observations behind ``slope``.
    """

    carriers: dict[MarkerKey, tuple[int, int]] = field(default_factory=dict)
    slope: float = 0.0
    dose_cap: int = DEFAULT_DOSE_CAP
    gated: bool = False
    gate_slope: float = 0.0
    gate_p_value: float = 1.0
    n_consensus: int = 0
    n_informative: int = 0

    def carrier_dose(self, m: InformativeMarker, host_allele: int) -> int:
        """Co-pooled carrier dose of ``host_allele`` (0=REF, 1=ALT) at marker ``m``.

        The raw cohort count of carriers of that allele, minus the host's and the
        (first) donor's own carriage, clamped to ``[0, dose_cap]``. For a donor-hom
        marker the host carries the host allele and the donor does not, so this
        removes exactly the one host contribution, leaving the co-pooled count.
        """
        n_ref, n_alt = self.carriers.get((m.chrom, m.pos, m.ref, m.alt), (0, 0))
        total = n_alt if host_allele == 1 else n_ref
        if host_allele in m.host_gt:
            total -= 1
        if host_allele in m.donor_gts[0]:
            total -= 1
        return max(0, min(self.dose_cap, total))


def apply_contamination_correction(
    markers: list[InformativeMarker],
    correction: ContaminationCorrection | None,
) -> list[InformativeMarker]:
    """Subtract dose-predicted contamination from donor-hom host-allele counts.

    Returns the input list unchanged (the same object) when there is nothing to
    do, so the default estimation path stays byte-identical: ``correction`` is
    None, the flowcell gated out, or the slope is non-positive. Otherwise returns
    a new list where each donor-homozygous marker (Vynck type 0 or 1) has
    ``slope * dose * depth`` reads removed from its host allele and its depth
    rebuilt from the two allele counts. Donor-het and non-donor-hom markers pass
    through untouched.
    """
    if correction is None or not correction.gated or correction.slope <= 0.0:
        return markers

    out: list[InformativeMarker] = []
    for m in markers:
        # Host allele = the allele the host carries that the donor lacks. For a
        # donor-hom marker that is REF (type 0: host 0/0, donor 1/1) or ALT
        # (type 1: host 1/1, donor 0/0). Contamination lands on those reads.
        if m.marker_type == 0:
            host_allele, reads = 0, m.admix_ad_ref
        elif m.marker_type == 1:
            host_allele, reads = 1, m.admix_ad_alt
        else:
            out.append(m)
            continue

        dose = correction.carrier_dose(m, host_allele)
        if dose <= 0:
            out.append(m)
            continue

        contam = correction.slope * dose * m.admix_dp
        new_reads = max(0, int(round(reads - contam)))
        if host_allele == 0:
            out.append(replace(m, admix_ad_ref=new_reads, admix_dp=new_reads + m.admix_ad_alt))
        else:
            out.append(replace(m, admix_ad_alt=new_reads, admix_dp=m.admix_ad_ref + new_reads))
    return out


@dataclass(frozen=True)
class SlopeFit:
    """Weighted-least-squares slope fit, with a one-sided significance test."""

    slope: float
    intercept: float
    se: float
    p_value: float  # one-sided P(slope > 0)
    n: int


def _wls_slope_fit(xs: list[float], ys: list[float], ws: list[float]) -> SlopeFit:
    """Weighted least squares of ``y ~ a + b x`` with a one-sided slope test.

    The single-group case of ``_grouped_slope_fit``; kept as the simple entry
    point. Weights are the per-observation depths (a proxy for inverse variance),
    the slope SE is ``sqrt(s2 / Sxx)`` on ``n - 2`` degrees of freedom, and the
    p-value is one-sided for a positive slope (the only direction contamination
    can take).
    """
    return _grouped_slope_fit([(xs, ys, ws)])


def _grouped_slope_fit(
    groups: list[tuple[list[float], list[float], list[float]]],
) -> SlopeFit:
    """One shared slope across groups, each with its own intercept.

    Used to pool a patient's serial timepoints: the host-allele fraction at
    donor-hom markers carries a per-timepoint host level (the group intercept)
    plus a shared dose-dependent contamination slope. Centering each group on its
    own weighted means before pooling removes the host level, so the slope is not
    biased by timepoints with different host fractions. A single group reduces to
    ordinary WLS.

    The slope is ``sum_g Sxy_g / sum_g Sxx_g`` (each centered within its group),
    its SE is ``sqrt(s2 / sum Sxx)`` with ``s2`` the pooled weighted residual
    variance on ``N - G - 1`` degrees of freedom (G group intercepts + 1 slope),
    and the p-value is one-sided for a positive slope. Returns a null fit when
    there are too few points or no spread in x.

    Args:
        groups: One ``(xs, ys, ws)`` triple per group.

    Returns:
        A ``SlopeFit`` (its ``intercept`` is the pooled weighted mean of y, for
        reference only; per-group intercepts are not returned).
    """
    sxy = 0.0
    sxx = 0.0
    n_total = 0
    n_groups = 0
    resid_ss = 0.0
    sum_wy = 0.0
    sum_w = 0.0
    arrays: list[tuple[np.ndarray, np.ndarray, np.ndarray, float, float]] = []
    for xs, ys, ws in groups:
        x = np.asarray(xs, dtype=float)
        y = np.asarray(ys, dtype=float)
        w = np.asarray(ws, dtype=float)
        if x.size < 2 or w.sum() <= 0.0:
            continue
        W = float(w.sum())
        xm = float((w * x).sum() / W)
        ym = float((w * y).sum() / W)
        g_sxx = float((w * (x - xm) ** 2).sum())
        if g_sxx <= 0.0:
            continue
        sxy += float((w * (x - xm) * (y - ym)).sum())
        sxx += g_sxx
        n_total += int(x.size)
        n_groups += 1
        sum_wy += float((w * y).sum())
        sum_w += W
        arrays.append((x, y, w, xm, ym))

    if sxx <= 0.0 or n_total == 0:
        return SlopeFit(0.0, 0.0, float("inf"), 1.0, n_total)

    b = sxy / sxx
    intercept = sum_wy / sum_w if sum_w > 0 else 0.0
    dof = n_total - n_groups - 1
    if dof <= 0:
        return SlopeFit(b, intercept, float("inf"), 1.0, n_total)

    for x, y, w, xm, ym in arrays:
        resid = y - ym - b * (x - xm)
        resid_ss += float((w * resid**2).sum())
    s2 = resid_ss / dof
    var_b = s2 / sxx
    if var_b <= 0.0:
        p = 0.0 if b > 0.0 else 1.0
        return SlopeFit(b, intercept, 0.0, p, n_total)

    se = math.sqrt(var_b)
    p = float(student_t.sf(b / se, dof))
    return SlopeFit(b, intercept, se, p, n_total)


def estimate_carrier_counts(
    genotype_lists: list[list[MarkerData]],
) -> dict[MarkerKey, tuple[int, int]]:
    """Count co-pooled carriers of each allele over a cohort.

    For every marker, counts how many cohort individuals carry the REF allele
    and how many carry the ALT allele (one count per individual per allele,
    regardless of het/hom). This is the dose lookup the correction subtracts
    against, the deployment-invariant part of the table, from the same
    joint-called cohort genotypes that feed the #23 error table.

    Returns a dict mapping marker key to ``(n_ref_carriers, n_alt_carriers)``.
    """
    n_ref: dict[MarkerKey, int] = {}
    n_alt: dict[MarkerKey, int] = {}
    for markers in genotype_lists:
        for m in markers:
            key = marker_key(m)
            if 0 in m.gt:
                n_ref[key] = n_ref.get(key, 0) + 1
            if 1 in m.gt:
                n_alt[key] = n_alt.get(key, 0) + 1
    keys = set(n_ref) | set(n_alt)
    return {k: (n_ref.get(k, 0), n_alt.get(k, 0)) for k in keys}


def _consensus_hom_points(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix: list[MarkerData],
    carriers: dict[MarkerKey, tuple[int, int]],
    dose_cap: int,
    min_dp: int,
) -> tuple[list[float], list[float], list[float]]:
    """Dose / minor-fraction / weight points at consensus-homozygous sites.

    A consensus-hom site is one where the host and every donor are homozygous for
    the same allele, so the admixture is a mixture of identical homozygotes and
    the minor allele cannot come from a contributor. Its dose is the co-pooled
    carrier count of the minor allele (host and donors do not carry it, so no
    exclusion is needed). Used for the per-flowcell gate.
    """
    hgt = {marker_key(m): m.gt for m in host}
    dgts = [{marker_key(m): m.gt for m in d} for d in donors]
    xs: list[float] = []
    ys: list[float] = []
    ws: list[float] = []
    for m in admix:
        key = marker_key(m)
        gh = hgt.get(key)
        if gh is None or gh not in ((0, 0), (1, 1)):
            continue
        consensus = True
        for dg in dgts:
            g = dg.get(key)
            if g is None or g != gh:
                consensus = False
                break
        if not consensus:
            continue
        dp = m.ad_ref + m.ad_alt
        if dp < min_dp or dp <= 0:
            continue
        if gh == (0, 0):
            minor_allele, minor_reads = 1, m.ad_alt
        else:
            minor_allele, minor_reads = 0, m.ad_ref
        n_ref, n_alt = carriers.get(key, (0, 0))
        total = n_alt if minor_allele == 1 else n_ref
        dose = max(0, min(dose_cap, total))
        xs.append(float(dose))
        ys.append(minor_reads / dp)
        ws.append(float(dp))
    return xs, ys, ws


def _informative_points(
    informative: list[InformativeMarker],
    correction: ContaminationCorrection,
) -> tuple[list[float], list[float], list[float]]:
    """Dose / host-fraction / weight points at donor-hom informative markers.

    The host-allele fraction at a donor-hom marker is host signal (constant
    across these markers) plus contamination (scaling with carrier dose). The
    weighted slope of that fraction on the dose is the contamination per carrier,
    isolating it from the constant host term in the intercept. Uses the same
    carrier dose (host/donor excluded) the correction applies.
    """
    xs: list[float] = []
    ys: list[float] = []
    ws: list[float] = []
    for m in informative:
        if m.marker_type == 0:
            host_allele, reads = 0, m.admix_ad_ref
        elif m.marker_type == 1:
            host_allele, reads = 1, m.admix_ad_alt
        else:
            continue
        if m.admix_dp <= 0:
            continue
        dose = correction.carrier_dose(m, host_allele)
        xs.append(float(dose))
        ys.append(reads / m.admix_dp)
        ws.append(float(m.admix_dp))
    return xs, ys, ws


def estimate_contamination_table(
    host: list[MarkerData],
    donors: list[list[MarkerData]],
    admix_lists: list[list[MarkerData]],
    cohort_genotypes: list[list[MarkerData]],
    *,
    min_dp: int = DEFAULT_MIN_DP,
    min_gq: int = DEFAULT_MIN_GQ,
    dose_cap: int = DEFAULT_DOSE_CAP,
    alpha: float = DEFAULT_GATE_ALPHA,
    min_slope: float = DEFAULT_GATE_MIN_SLOPE,
) -> ContaminationCorrection:
    """Build a per-marker contamination correction for one patient on a flowcell.

    Estimates the co-pooled carrier counts from the cohort genotypes, gates on
    the per-flowcell consensus-hom dose response, and (when gated in) calibrates
    the correction slope on the patient's own informative donor-hom markers,
    pooling across the supplied serial timepoints to beat per-sample noise.

    Args:
        host: Parsed host markers.
        donors: One parsed marker list per donor.
        admix_lists: One parsed marker list per admix timepoint for this patient
            (parse with ``min_dp=0``; filtering is applied here).
        cohort_genotypes: One parsed marker list per co-pooled cohort individual
            (the flowcell's joint-called genotypes), used for carrier counts.
            Expected to include the host and donors.
        min_dp: Minimum admix depth for a marker to contribute.
        min_gq: Minimum host/donor genotype quality for classification.
        dose_cap: Carrier-dose cap (see ``DEFAULT_DOSE_CAP``).
        alpha: Significance level for the consensus-hom slope gate.
        min_slope: Minimum consensus-hom slope (per carrier) worth correcting.

    Returns:
        A ``ContaminationCorrection``. ``gated`` is False (and ``slope`` 0) when
        the consensus-hom dose response is not significantly positive, which
        leaves estimation unchanged.
    """
    carriers = estimate_carrier_counts(cohort_genotypes)

    # Each timepoint is its own group, so the per-timepoint host level is absorbed
    # by a group intercept and only the shared dose-dependent contamination slope
    # is pooled (see ``_grouped_slope_fit``). Pooling raw across timepoints would
    # let their differing host fractions leak into the slope.
    gate_groups: list[tuple[list[float], list[float], list[float]]] = []
    info_groups: list[tuple[list[float], list[float], list[float]]] = []
    probe = ContaminationCorrection(carriers=carriers, dose_cap=dose_cap)
    for admix in admix_lists:
        gate_groups.append(_consensus_hom_points(host, donors, admix, carriers, dose_cap, min_dp))
        genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
        info_groups.append(_informative_points(genotypes.informative, probe))

    # Gate: per-flowcell consensus-hom dose response. Significant (guards against
    # low-depth noise) and above a minimum effect (guards against a trivially
    # small but high-depth-significant slope).
    gate = _grouped_slope_fit(gate_groups)
    gated = (gate.p_value < alpha) and (gate.slope > min_slope)

    # Magnitude: informative donor-hom slope, calibrated on the markers we
    # actually correct (carrier dose excludes host/donor). Computed regardless of
    # the gate so its provenance is recorded; only applied when gated in.
    mag = _grouped_slope_fit(info_groups)

    slope = max(0.0, mag.slope) if gated else 0.0
    return ContaminationCorrection(
        carriers=carriers,
        slope=slope,
        dose_cap=dose_cap,
        gated=gated,
        gate_slope=gate.slope,
        gate_p_value=gate.p_value,
        n_consensus=gate.n,
        n_informative=mag.n,
    )


_HEADER_KEYS = (
    "slope",
    "dose_cap",
    "gated",
    "gate_slope",
    "gate_p_value",
    "n_consensus",
    "n_informative",
)


def save_contamination_table(correction: ContaminationCorrection, path: Path | str) -> None:
    """Write a contamination correction to a TSV file.

    The scalar gate/slope fields are written as ``# key\\tvalue`` comment lines,
    followed by a ``chrom pos ref alt n_ref n_alt`` carrier-count table. The
    comment block is parsed back by ``load_contamination_table``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# allomix contamination table (Step 30, issue #30)\n")
        fh.write(f"# slope\t{correction.slope:.8e}\n")
        fh.write(f"# dose_cap\t{correction.dose_cap}\n")
        fh.write(f"# gated\t{1 if correction.gated else 0}\n")
        fh.write(f"# gate_slope\t{correction.gate_slope:.8e}\n")
        fh.write(f"# gate_p_value\t{correction.gate_p_value:.8e}\n")
        fh.write(f"# n_consensus\t{correction.n_consensus}\n")
        fh.write(f"# n_informative\t{correction.n_informative}\n")
        fh.write("chrom\tpos\tref\talt\tn_ref\tn_alt\n")
        for key in sorted(correction.carriers.keys()):
            chrom, pos, ref, alt = key
            n_ref, n_alt = correction.carriers[key]
            fh.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{n_ref}\t{n_alt}\n")


def load_contamination_table(path: Path | str) -> ContaminationCorrection:
    """Load a contamination correction written by ``save_contamination_table``."""
    meta: dict[str, str] = {}
    carriers: dict[MarkerKey, tuple[int, int]] = {}
    with open(path, encoding="utf-8") as fh:
        header_seen = False
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split("\t")
                if len(parts) == 2:
                    meta[parts[0].strip()] = parts[1].strip()
                continue
            fields = line.split("\t")
            if not header_seen:
                # The column header row (chrom pos ref alt n_ref n_alt).
                header_seen = True
                continue
            chrom, pos, ref, alt, n_ref, n_alt = fields
            carriers[(chrom, int(pos), ref, alt)] = (int(n_ref), int(n_alt))

    return ContaminationCorrection(
        carriers=carriers,
        slope=float(meta.get("slope", 0.0)),
        dose_cap=int(meta.get("dose_cap", DEFAULT_DOSE_CAP)),
        gated=meta.get("gated", "0") == "1",
        gate_slope=float(meta.get("gate_slope", 0.0)),
        gate_p_value=float(meta.get("gate_p_value", 1.0)),
        n_consensus=int(meta.get("n_consensus", 0)),
        n_informative=int(meta.get("n_informative", 0)),
    )


__all__ = [
    "ContaminationCorrection",
    "SlopeFit",
    "apply_contamination_correction",
    "estimate_carrier_counts",
    "estimate_contamination_table",
    "load_contamination_table",
    "save_contamination_table",
]
