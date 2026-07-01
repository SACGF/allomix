"""Per-marker panel inclusion QC (issue #37).

Turns the per-marker characterisation from ``scripts/measure_panel_bias.py``
(``panel_stats_per_marker.tsv``) into a keep/drop verdict per marker, with an
auditable reason for every dropped marker and a panel-level summary. This
automates step 4 of ``docs/panel_guide.md``, which was previously a manual pass
over the statistics.

The cutoffs are starting points from the panel guide, not validated thresholds.
Every one is a CLI flag, and a sample-ID panel already validated for identity
work should need little pruning.

Criteria (see ``PanelQCThresholds``):

- ``low_call_rate`` - ``call_rate`` below ``min_call_rate``; the marker fails to
  genotype too often to add signal.
- ``extreme_bias`` - ``|median_bias|`` beyond the panel's extreme tail (a
  multiple of the panel ``p95_abs_bias``, or an absolute cap). Moderate bias is
  left to the bias-correction step, so only the extreme tail is excluded.
- ``het_deficit`` - ``het_ratio_vs_hwe`` below ``min_het_ratio`` (a het deficit
  versus Hardy-Weinberg), an allele-dropout signal that bias correction does not
  fix.
- ``unstable_depth`` - high ``depth_cv``. Secondary: it never drops a marker on
  its own, it is only recorded as an extra reason on a marker already dropped by
  one of the criteria above.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

# Default cutoffs. Starting points from docs/panel_guide.md step 4, NOT validated
# thresholds; every one is a CLI flag on the panel-qc subcommand.
DEFAULT_MIN_CALL_RATE = 0.90  # drop markers genotyped in <90% of samples
DEFAULT_MIN_HET_RATIO = 0.80  # het deficit vs HWE below this flags allele dropout
DEFAULT_BIAS_P95_MULT = 2.0  # |bias| beyond this multiple of panel p95 is extreme
DEFAULT_MAX_DEPTH_CV = 1.0  # secondary: annotates an already-dropped marker

# Reason codes recorded per dropped marker (one column entry each).
REASON_LOW_CALL_RATE = "low_call_rate"
REASON_EXTREME_BIAS = "extreme_bias"
REASON_HET_DEFICIT = "het_deficit"
REASON_UNSTABLE_DEPTH = "unstable_depth"

# Primary reasons drive keep/drop; the secondary reason only annotates.
PRIMARY_REASONS = (REASON_LOW_CALL_RATE, REASON_EXTREME_BIAS, REASON_HET_DEFICIT)
ALL_REASONS = (*PRIMARY_REASONS, REASON_UNSTABLE_DEPTH)


@dataclass
class PanelQCThresholds:
    """Tunable per-marker inclusion cutoffs. All values are starting points."""

    min_call_rate: float = DEFAULT_MIN_CALL_RATE
    min_het_ratio: float = DEFAULT_MIN_HET_RATIO
    bias_p95_mult: float = DEFAULT_BIAS_P95_MULT
    max_abs_bias: float | None = None
    max_depth_cv: float = DEFAULT_MAX_DEPTH_CV


@dataclass
class MarkerVerdict:
    """A single marker's keep/drop decision and the reasons behind it."""

    marker: str | None  # chrom:pos:ref:alt when the TSV carries coordinates
    marker_index: int
    row: dict[str, str]  # original TSV row, preserved for the verdict output
    keep: bool
    reasons: list[str] = field(default_factory=list)


def _percentile(sorted_vals: list[float], pct: float) -> float:
    """Percentile from a pre-sorted list.

    Matches ``scripts/measure_panel_bias.py`` so the panel ``p95_abs_bias``
    recomputed here agrees with the value that script reports.
    """
    if not sorted_vals:
        return float("nan")
    idx = min(len(sorted_vals) - 1, int(pct / 100.0 * len(sorted_vals)))
    return sorted_vals[idx]


def _to_float(value: str | None) -> float | None:
    """Parse a TSV cell to float, treating empty/``NA`` as missing."""
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.upper() == "NA":
        return None
    return float(value)


def load_marker_stats(path: Path | str) -> list[dict[str, str]]:
    """Load ``panel_stats_per_marker.tsv`` rows as dicts (columns by name)."""
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"no marker rows in {path}")
    return rows


def panel_p95_abs_bias(rows: list[dict[str, str]]) -> float:
    """95th percentile of ``|median_bias|`` across markers, or NaN if none.

    The extreme-bias criterion is relative to this panel-wide tail, so it is
    computed from the same rows rather than read from the facts CSV, keeping
    panel-qc self-contained.
    """
    abs_biases = sorted(
        abs(b) for b in (_to_float(r.get("median_bias")) for r in rows) if b is not None
    )
    return _percentile(abs_biases, 95)


def assign_verdict(
    row: dict[str, str],
    index: int,
    thresholds: PanelQCThresholds,
    p95_abs_bias: float,
) -> MarkerVerdict:
    """Apply the cutoffs to one marker row and return its verdict.

    A missing (``NA``) statistic never triggers its criterion: the marker is not
    dropped for a metric that could not be measured.
    """
    call_rate = _to_float(row.get("call_rate"))
    het_ratio = _to_float(row.get("het_ratio_vs_hwe"))
    median_bias = _to_float(row.get("median_bias"))
    depth_cv = _to_float(row.get("depth_cv"))

    # Absolute cap overrides the panel-relative multiple when supplied.
    if thresholds.max_abs_bias is not None:
        bias_cap: float = thresholds.max_abs_bias
    elif not math.isnan(p95_abs_bias):
        bias_cap = thresholds.bias_p95_mult * p95_abs_bias
    else:
        bias_cap = float("nan")

    reasons: list[str] = []
    if call_rate is not None and call_rate < thresholds.min_call_rate:
        reasons.append(REASON_LOW_CALL_RATE)
    if median_bias is not None and not math.isnan(bias_cap) and abs(median_bias) > bias_cap:
        reasons.append(REASON_EXTREME_BIAS)
    if het_ratio is not None and het_ratio < thresholds.min_het_ratio:
        reasons.append(REASON_HET_DEFICIT)

    # Secondary: unstable depth only annotates a marker already dropped by a
    # primary criterion, it never drops one on its own (issue #37).
    drop = bool(reasons)
    if drop and depth_cv is not None and depth_cv > thresholds.max_depth_cv:
        reasons.append(REASON_UNSTABLE_DEPTH)

    marker = row.get("marker") or None
    return MarkerVerdict(
        marker=marker,
        marker_index=index,
        row=row,
        keep=not drop,
        reasons=reasons,
    )


def assign_verdicts(
    rows: list[dict[str, str]], thresholds: PanelQCThresholds
) -> tuple[list[MarkerVerdict], float]:
    """Assign verdicts to every marker. Returns (verdicts, panel p95_abs_bias)."""
    p95 = panel_p95_abs_bias(rows)
    verdicts = [assign_verdict(row, i, thresholds, p95) for i, row in enumerate(rows)]
    return verdicts, p95


def summarize(verdicts: list[MarkerVerdict]) -> dict[str, int]:
    """Panel-level counts: totals plus dropped-marker counts per reason.

    A dropped marker with several reasons is counted once under each of them, so
    the per-reason counts can sum to more than ``n_drop``.
    """
    summary: dict[str, int] = {
        "n_markers": len(verdicts),
        "n_keep": sum(1 for v in verdicts if v.keep),
        "n_drop": sum(1 for v in verdicts if not v.keep),
    }
    for reason in ALL_REASONS:
        summary[reason] = sum(1 for v in verdicts if reason in v.reasons)
    return summary


def write_verdict_tsv(verdicts: list[MarkerVerdict], path: Path | str) -> None:
    """Write the per-marker verdict TSV: original columns plus verdict/reasons."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Preserve the input column order, then append the two decision columns.
    base_cols = list(verdicts[0].row.keys())
    fieldnames = [*base_cols, "verdict", "reasons"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for v in verdicts:
            out = dict(v.row)
            out["verdict"] = "keep" if v.keep else "drop"
            out["reasons"] = ",".join(v.reasons) if v.reasons else "."
            writer.writerow(out)


def _marker_to_bed_row(marker: str) -> tuple[str, int, int]:
    """Parse a ``chrom:pos:ref:alt`` marker into a 0-based BED interval."""
    chrom, pos_str, _ref, _alt = marker.split(":", 3)
    pos = int(pos_str)  # 1-based VCF position
    return chrom, pos - 1, pos


def write_sites_bed(verdicts: list[MarkerVerdict], path: Path | str, *, which: str) -> int:
    """Write a sites BED of dropped (``which='drop'``) or kept markers.

    Returns the number of intervals written. Requires the ``marker`` coordinate
    column: raises ``ValueError`` if any selected marker lacks coordinates, since
    a BED that silently omits markers cannot be fed to ``allomix detect`` safely.
    """
    if which not in ("drop", "keep"):
        raise ValueError(f"which must be 'drop' or 'keep', got {which!r}")
    selected = [v for v in verdicts if (v.keep if which == "keep" else not v.keep)]
    missing = [v.marker_index for v in selected if not v.marker]
    if missing:
        raise ValueError(
            "cannot write a sites BED: the per-marker TSV has no 'marker' "
            "coordinate column (marker_index "
            f"{missing[0]} and {len(missing) - 1} more). Re-run "
            "scripts/measure_panel_bias.py to emit coordinates."
        )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(_marker_to_bed_row(v.marker) for v in selected)  # type: ignore[arg-type]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for chrom, start, end in rows:
            writer.writerow([chrom, start, end])
    return len(rows)
