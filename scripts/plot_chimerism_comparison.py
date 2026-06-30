#!/usr/bin/env python3
"""Plot allomix whole-blood chimerism against flow-sorted lineage values.

Validation plot (run from `scripts/`, not part of the installed allomix package).

Takes one or two allomix `batch.tsv` files (as produced by
`scripts/run_csv_batch.py`) and draws a per-sample forest plot. Each sample
shows the NGS donor estimate with its confidence interval, overlaid on the
flow cytometry lineage values (CD45 / CD3 / CD13) parsed from a copied column.

The y axis shows donor %, but is log-spaced by distance from 100% donor and
inverted so 100% is at the top. The clinically interesting action is the
low-level signal near full donor chimerism, which a plain linear (or plain log)
donor axis compresses; this spacing keeps it readable. The spread of the
CD3/CD13 subsets shows where the true whole-blood value can sit depending on the
cell differential. QC-FAIL samples are skipped; QC-REVIEW samples are drawn but
circled in red so a confident-looking estimate is not read as clean.

Convention: all percentages are % DONOR (allomix `donor_pct` and the flow
values).

Usage:
    # Single run (run2), flow column already merged into the batch.tsv:
    python scripts/plot_chimerism_comparison.py \
        output/validation_run2/batch.tsv \
        --flow-column "Chimerism result TP2" \
        --output output/chimerism_comparison.png

    # Compare runs (primary drawn rightmost/filled; compares listed oldest first;
    # labels are left to right):
    python scripts/plot_chimerism_comparison.py \
        output/validation_run3/batch.tsv \
        --compare-tsv output/validation_run1/batch.tsv output/validation_run2/batch.tsv \
        --labels run1 run2 run3 \
        --output output/run1_vs_run2_vs_run3.png
"""

import argparse
import csv
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

# Markers per flow lineage. CD45 is the pan-leukocyte value (the nominal
# whole-blood comparator); CD3 (T cells) and CD13 (myeloid) are subsets whose
# spread brackets the true whole-blood value.
LINEAGE_MARKERS = {"CD45": "s", "CD3": "^", "CD13": "o"}
_LINEAGE_RE = re.compile(r"([A-Za-z0-9]+)\s+([\d.]+)\s*%")


def parse_lineages(text: str) -> dict[str, float]:
    """Parse a flow lineage string into {lineage: donor_pct}.

    Handles the mixed `;` and `,` delimiters seen in the source sheet, e.g.
    "CD45 100%; CD3 93.19%, CD13 30.58%".

    Args:
        text: Raw cell value.

    Returns:
        Mapping of lineage name to donor percentage. Empty if nothing parses.
    """
    if not text:
        return {}
    return {name: float(val) for name, val in _LINEAGE_RE.findall(text)}


def _qc_status(row: dict, gof: float | None) -> str:
    """Three-state QC status for a batch.tsv row.

    Uses the ``qc_status`` column when present (allomix now emits PASS / REVIEW /
    FAIL). For older files that only have the binary ``qc_pass`` column, derive
    REVIEW from a failing goodness-of-fit so a poor-fit sample is still flagged.

    Args:
        row: One batch.tsv row.
        gof: Parsed goodness-of-fit p-value, or None.

    Returns:
        One of "PASS", "REVIEW", "FAIL".
    """
    status = row.get("qc_status")
    if status:
        return status.strip().upper()
    if row.get("qc_pass", "PASS") != "PASS":
        return "FAIL"
    if gof is not None and gof < 0.01:
        return "REVIEW"
    return "PASS"


def read_batch(path: Path, flow_column: str | None, donor_column: str = "Donor") -> dict[str, dict]:
    """Read a batch.tsv into {sample: row dict} with parsed numeric fields.

    Args:
        path: Path to a batch.tsv.
        flow_column: Name of the column holding flow lineage strings, or None.
        donor_column: Name of the column describing the donor type.

    Returns:
        Mapping of sample name to a dict with donor_pct, ci_lo, ci_hi, qc_pass,
        donor (type string), lineages (dict), and optional lod_pct.
    """
    out: dict[str, dict] = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            gof = None
            if row.get("gof_pval") not in (None, "", "NA"):
                gof = float(row["gof_pval"])
            rec: dict = {
                "donor_pct": float(row["donor_pct"]),
                "ci_lo": float(row["ci_lo"]),
                "ci_hi": float(row["ci_hi"]),
                "n_informative": int(row.get("n_informative", 0) or 0),
                # n_total_markers is present only on batch.tsv files produced after
                # the total-marker-count change; 0 means "not recorded" and the
                # label falls back to showing the informative count alone.
                "n_total_markers": int(row.get("n_total_markers", 0) or 0),
                "mean_depth": float(row.get("mean_depth", 0) or 0),
                "qc_status": _qc_status(row, gof),
                "donor": (row.get(donor_column) or "").strip(),
                "lineages": {},
            }
            rec["qc_pass"] = rec["qc_status"] != "FAIL"
            # lob_pct / lod_pct are present only on batch.tsv files produced
            # after the per-sample LOD change; tolerate their absence.
            if row.get("lod_pct") not in (None, "", "NA"):
                rec["lod_pct"] = float(row["lod_pct"])
            # Host-presence p-value, present only on pileup-mode batch.tsv files.
            if row.get("host_present_p") not in (None, "", "NA"):
                rec["host_present_p"] = float(row["host_present_p"])
            # Host-presence fraction estimate and CI (fractions in the file; stored
            # as host %, the same distance-from-100 axis as the donor points).
            if row.get("host_f_est") not in (None, "", "NA"):
                rec["host_f_est"] = float(row["host_f_est"]) * 100.0
                rec["host_f_ci_lo"] = float(row["host_f_ci_lo"]) * 100.0
                rec["host_f_ci_hi"] = float(row["host_f_ci_hi"]) * 100.0
            if flow_column and row.get(flow_column):
                rec["lineages"] = parse_lineages(row[flow_column])
            out[sample] = rec
    return out


def host(donor_pct: float) -> float:
    """Convert a donor percentage to host percentage."""
    return 100.0 - donor_pct


def _marker_count(rec: dict) -> str:
    """Format a run's marker count as informative/total for the x-axis label.

    The total marker count is one of the things that changes between runs (panel
    coverage, locus dropout), so showing the informative count as a fraction of
    the total makes the panel difference legible. Older batch.tsv files without a
    total fall back to the informative count alone.

    Args:
        rec: One run's batch row.

    Returns:
        e.g. "34/76" when a total is recorded, else "34".
    """
    n_inf = rec["n_informative"]
    n_total = rec.get("n_total_markers", 0)
    return f"{n_inf}/{n_total}" if n_total else f"{n_inf}"


def short_label(name: str, field: int | None, code: bool) -> str:
    """Shorten a long sample name for display.

    Args:
        name: Full sample name, e.g. "14_MO_IDH_APM5_REDACTED_REDACTED".
        field: If set, split on "_" and take this 0-based index.
        code: If True, take the last all-uppercase-letters token (robust to the
            patient code sitting at field 3 or 4). Wins over ``field``.

    Returns:
        The shortened label, or the full name if nothing suitable is found.
    """
    if code:
        tokens = [t for t in name.split("_") if t.isalpha() and t.isupper()]
        return tokens[-1] if tokens else name
    if field is not None:
        tokens = name.split("_")
        return tokens[field] if field < len(tokens) else name
    return name


# Run colours, assigned by display position (leftmost run first): run1 orange,
# run2 blue, run3 green, ... The primary run (the rightmost) is drawn filled and
# the compare runs hollow, so the fill (not the colour) marks the primary.
RUN_PALETTE = ["#ff7f0e", "#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
LINEAGE_COLOR = "#888888"
REVIEW_COLOR = "#d62728"  # ring around QC-REVIEW primary points
PRESENCE_DET = "#1b7837"  # host-presence detected (p <= 0.05)
PRESENCE_NOT = "#999999"  # host-presence not detected


def _fmt_p(p: float) -> str:
    """Compact presence p-value: bound on underflow, decimal near 1, else sci."""
    if p <= 1e-300:
        return "<1e-300"
    return f"{p:.2f}" if p >= 0.01 else f"{p:.0e}"


def plot(
    runs: list[dict[str, dict]],
    labels: list[str],
    floor: float,
    title: str,
    output: Path,
    anonymize: bool,
    label_field: int | None,
    label_code: bool,
    explicit_donor: set[str],
    sort: str,
    show_lod: bool,
) -> None:
    """Draw the forest plot and write it to ``output``.

    The y axis shows donor %, but is log-spaced by distance from 100% donor
    (internally 100 - donor%) and inverted so 100% sits at the top. A plain log
    of donor% would compress everything near 100% where the low-level signal
    lives; this keeps that region readable while labelling the axis in donor %.

    Each x-axis label row carries the informative/total marker count so the panel
    difference between runs is visible at a glance.

    Args:
        runs: One or more batch dicts in left-to-right display order. The last
            entry is the primary run (drawn filled); flow lineages, donor type,
            sample ordering and the explicit-donor star are all read from it.
            Earlier entries are compare runs, drawn hollow.
        labels: Legend label for each run, parallel to ``runs``.
        show_lod: If True, shade each sample's region at or below the primary
            run's limit of detection (where a signal is not a reportable
            detection). Needs an ``lod_pct`` column in the primary batch.tsv.
    """
    primary = runs[-1]
    n_runs = len(runs)
    run_colors = [RUN_PALETTE[i % len(RUN_PALETTE)] for i in range(n_runs)]

    # Sample order on the x axis. "chimerism" (sort by measured donor fraction)
    # reshuffles between runs because the estimate differs, so "name" and "tsv"
    # give a stable order that lets the runs' plots line up.
    match sort:
        case "name":
            ordered = sorted(primary)
        case "tsv":
            ordered = list(primary)
        case _:  # "chimerism"
            ordered = sorted(primary, key=lambda s: host(primary[s]["donor_pct"]))
    # QC-FAIL samples (e.g. a single informative marker) have point estimates
    # too unreliable to plot, so drop them rather than imply a real measurement.
    # QC-REVIEW samples are kept (and circled later) so they are not hidden.
    skipped = [s for s in ordered if not primary[s]["qc_pass"]]
    samples = [s for s in ordered if primary[s]["qc_pass"]]
    if skipped:
        print(f"Skipped {len(skipped)} QC-FAIL sample(s): {', '.join(skipped)}")
    n = len(samples)

    def code_for(i: int, s: str) -> str:
        return f"S{i + 1}" if anonymize else short_label(s, label_field, label_code)

    def clamp(v: float) -> float:
        return max(v, floor)

    fig, ax = plt.subplots(figsize=(max(9.0, 1.45 * n), 6.5))

    # Dodge: flow markers sit far left, then the runs spread across the slot in
    # display order so time reads left to right. A single run is centred.
    if n_runs == 1:
        flow_dx = -0.22
        offsets = [0.0]
    else:
        flow_dx = -0.40
        step = 0.16
        start = -((n_runs - 1) / 2.0) * step + 0.06
        offsets = [start + i * step for i in range(n_runs)]

    # Per-sample LOD band: shade the strip from 100% donor down to the primary
    # run's limit of detection. A point inside the band is below LOD, i.e. not a
    # reportable detection (statistically indistinguishable from full donor).
    review_drawn = False
    lod_drawn = False
    if show_lod:
        top_y = floor * 0.7
        for x, sample in enumerate(samples):
            lod = primary[sample].get("lod_pct")
            if lod is None:
                continue
            ax.fill_between(
                [x - 0.46, x + 0.46],
                clamp(lod),
                top_y,
                color="0.6",
                alpha=0.13,
                lw=0,
                zorder=0,
            )
            lod_drawn = True

    for x, sample in enumerate(samples):
        # Flow lineage markers, offset left of centre, each at its own donor
        # value (100% sits at the top). A thin line shows the lineage spread.
        lin = primary[sample]["lineages"]
        if lin:
            lx = x + flow_dx
            host_vals = [clamp(host(v)) for v in lin.values()]
            if max(host_vals) > min(host_vals):
                ax.plot(
                    [lx, lx],
                    [min(host_vals), max(host_vals)],
                    color=LINEAGE_COLOR,
                    lw=1.0,
                    zorder=1,
                )
            for name, donor in lin.items():
                ax.scatter(
                    lx,
                    clamp(host(donor)),
                    marker=LINEAGE_MARKERS.get(name, "d"),
                    s=42,
                    facecolors="none",
                    edgecolors=LINEAGE_COLOR,
                    zorder=3,
                )

        # NGS estimate per run, with CI error bars. The primary is always drawn
        # (samples are filtered to its QC pass). A compare run that failed or
        # found no informative markers has no point to plot (its 0 still shows
        # in the x-axis label row).
        for k, run in enumerate(runs):
            rec = run.get(sample)
            if rec is None:
                continue
            is_primary = k == n_runs - 1
            if not is_primary and (not rec["qc_pass"] or rec["n_informative"] == 0):
                continue
            _ngs_point(ax, x + offsets[k], rec, run_colors[k], clamp, hollow=not is_primary)
            # Circle a primary-run point flagged for QC review (e.g. poor model
            # fit), so a confident-looking estimate is not read as clean.
            if is_primary and rec.get("qc_status") == "REVIEW":
                ax.scatter(
                    x + offsets[k],
                    clamp(host(rec["donor_pct"])),
                    s=200,
                    facecolors="none",
                    edgecolors=REVIEW_COLOR,
                    linewidths=1.6,
                    zorder=5,
                )
                review_drawn = True

            # Host-presence estimate as its own dot + CI, dodged just right of the
            # primary point. Green when host is detected (p <= 0.05), grey when
            # not; a not-detected CI runs up to the 100%-donor line (host = 0).
            # Near full donor this estimate is tighter than the donor MLE (it uses
            # the donor-homozygous markers directly), so it can be a confident
            # small detection where the donor CI is uninformative. The clean
            # full-donor samples have a degenerate zero estimate (nothing to draw).
            if is_primary and rec.get("host_f_ci_hi", 0.0) > 0:
                detected = rec.get("host_present_p", 1.0) <= 0.05
                _host_point(ax, x + offsets[k] + 0.12, rec, clamp, detected)

    # Reference lines at clinically relevant donor fractions (99%, 95%).
    for donor_thr in (99.0, 95.0):
        ax.axhline(host(donor_thr), color="0.8", ls="--", lw=0.8, zorder=0)
        ax.text(
            n - 0.4,
            host(donor_thr),
            f" {donor_thr:g}% donor",
            va="bottom",
            ha="right",
            color="0.5",
            fontsize=8,
        )

    ax.set_yscale("log")
    ax.set_ylim(105, floor * 0.7)  # inverted: 100% donor at top, 0% at bottom
    ax.set_yticks([0.1, 1, 10, 100])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{100 - v:g}"))
    ax.set_ylabel("Donor fraction (%), log-spaced by distance from 100%")
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_title(title)

    # Per-point x-axis labels. One row per run, colour-matched to the runs and
    # stacked in display order (leftmost run on top, primary at the bottom). The
    # patient code (with the explicit-donor star) is shown once, on the primary
    # row; the other rows carry just that run's marker count and depth.
    ax.set_xticks([])
    trans = ax.get_xaxis_transform()
    row0_y = -0.02
    row_step = 0.05
    # Row y-positions below the axis: the run rows, then (optionally) the
    # host-presence p row, then donor type, then key.
    has_presence = any("host_present_p" in primary.get(s, {}) for s in samples)
    last_run_y = row0_y - row_step * (n_runs - 1)
    presence_y = last_run_y - 0.05
    donor_y = (presence_y - 0.05) if has_presence else (last_run_y - 0.065)
    xlabel_y = donor_y - 0.135
    for i, s in enumerate(samples):
        code = code_for(i, s)
        for k, run in enumerate(runs):
            rec = run.get(s)
            md = f"M:{_marker_count(rec)} D:{rec['mean_depth']:.0f}x" if rec else "NA"
            is_primary = k == n_runs - 1
            if is_primary:
                # Flag explicit-donor samples with a star on the code (keyed in
                # the title); it marks the primary run.
                mark = "★" if any(tok and tok in s for tok in explicit_donor) else ""
                txt = f"{code}{mark} {md}"
            else:
                txt = md
            ax.text(
                i,
                row0_y - row_step * k,
                txt,
                transform=trans,
                ha="center",
                va="top",
                fontsize=8,
                # A single run keeps the default black label for readability.
                color=run_colors[k] if n_runs > 1 else "black",
            )
        # Host-presence p-value row (primary run only), coloured by the call.
        if has_presence:
            p = primary.get(s, {}).get("host_present_p")
            if p is not None:
                ax.text(
                    i,
                    presence_y,
                    _fmt_p(p),
                    transform=trans,
                    ha="center",
                    va="top",
                    fontsize=7.5,
                    color=PRESENCE_DET if p <= 0.05 else PRESENCE_NOT,
                )

        # Donor type (consistent across runs, so read once from primary).
        donor = primary[s]["donor"]
        if donor:
            ax.text(
                i,
                donor_y,
                "\n".join(textwrap.wrap(donor, 16)),
                transform=trans,
                ha="center",
                va="top",
                fontsize=6.5,
                style="italic",
                color="0.4",
            )

    # Key for the tick-label format, below the rows.
    xlabel = "Sample: M = informative/total markers, D = mean depth"
    if has_presence:
        xlabel += "   |   p = host-presence test (green = host detected, grey = not)"
    ax.text((n - 1) / 2.0, xlabel_y, xlabel, transform=trans, ha="center", va="top", fontsize=9)

    # Legend order matches the left-to-right draw order: flow, then the runs.
    handles = [
        Line2D([], [], color=LINEAGE_COLOR, marker="s", mfc="none", lw=0, label="Flow CD45"),
        Line2D([], [], color=LINEAGE_COLOR, marker="^", mfc="none", lw=0, label="Flow CD3"),
        Line2D([], [], color=LINEAGE_COLOR, marker="o", mfc="none", lw=0, label="Flow CD13"),
    ]
    for k in range(n_runs):
        is_primary = k == n_runs - 1
        handles.append(
            Line2D(
                [],
                [],
                color=run_colors[k],
                marker="o",
                mfc=run_colors[k] if is_primary else "none",
                lw=0,
                label=f"NGS {labels[k]} (CI)",
            )
        )
    if review_drawn:
        handles.append(
            Line2D(
                [],
                [],
                color=REVIEW_COLOR,
                marker="o",
                mfc="none",
                mew=1.6,
                ms=11,
                lw=0,
                label="QC review (e.g. poor fit)",
            )
        )
    if lod_drawn:
        handles.append(
            Patch(facecolor="0.6", alpha=0.13, label=f"≤ {labels[-1]} LOD (not detected)")
        )
    if has_presence:
        handles.append(
            Line2D(
                [],
                [],
                color=PRESENCE_DET,
                marker="D",
                mfc=PRESENCE_DET,
                ms=6,
                lw=0,
                label="host presence est. (detected)",
            )
        )
        handles.append(
            Line2D(
                [],
                [],
                color=PRESENCE_NOT,
                marker="D",
                mfc=PRESENCE_NOT,
                ms=6,
                lw=0,
                label="host presence est. (not detected)",
            )
        )
    # Legend below the x-axis label rows, laid out horizontally.
    ax.legend(
        handles=handles,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, xlabel_y - 0.05),
        ncol=min(len(handles), 4),
        framealpha=0.9,
    )

    fig.tight_layout()
    # Reserve room for the custom label rows + donor type + key + the legend now
    # sitting below them (tight_layout ignores all of these); more rows need more.
    fig.subplots_adjust(bottom=min(0.58, 0.34 + 0.05 * n_runs))
    fig.savefig(output, dpi=150)
    print(f"Wrote {output}")


def _ngs_point(ax, x, rec, color, clamp, hollow=False) -> None:
    """Plot one NGS estimate as a point with asymmetric CI error bars."""
    y = clamp(host(rec["donor_pct"]))
    y_lo = clamp(host(rec["ci_hi"]))  # upper donor CI -> smaller distance
    y_hi = clamp(host(rec["ci_lo"]))  # lower donor CI -> larger distance
    ax.errorbar(
        x,
        y,
        yerr=[[max(0.0, y - y_lo)], [max(0.0, y_hi - y)]],
        fmt="o",
        color=color,
        mfc="none" if hollow else color,
        capsize=3,
        ms=7,
        zorder=4,
    )


def _host_point(ax, x, rec, clamp, detected: bool) -> None:
    """Plot the host-presence estimate as a dot + asymmetric CI.

    The estimate and CI are host fraction (% distance from 100% donor), the same
    axis as the donor points, so this dot reads as a second, more sensitive
    measurement of the same quantity. A diamond marker keeps it distinct from the
    round per-run NGS dots. Green when host is detected (p <= 0.05), grey when
    not; a not-detected CI reaches the floor (the 100%-donor line, host = 0).
    """
    color = PRESENCE_DET if detected else PRESENCE_NOT
    y = clamp(rec["host_f_est"])
    y_lo = clamp(rec["host_f_ci_lo"])
    y_hi = clamp(rec["host_f_ci_hi"])
    ax.errorbar(
        x,
        y,
        yerr=[[max(0.0, y - y_lo)], [max(0.0, y_hi - y)]],
        fmt="D",
        color=color,
        mfc=color,
        ms=5,
        capsize=3,
        zorder=5,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "batch_tsv",
        type=Path,
        help="Primary batch.tsv (drawn filled, rightmost; flow/donor/star read from it)",
    )
    parser.add_argument(
        "--compare-tsv",
        type=Path,
        nargs="+",
        default=None,
        help="One or more batch.tsv files to overlay (hollow). Listed left to "
        "right, the primary is drawn to their right, so pass them oldest first.",
    )
    parser.add_argument(
        "--flow-column",
        default="Chimerism result TP2",
        help="Column name holding flow lineage strings",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        metavar="LABEL",
        help="Legend label per run in left-to-right display order (compare runs "
        "first, then primary). Defaults to run1, run2, ... Must match the run count.",
    )
    parser.add_argument(
        "--floor", type=float, default=0.02, help="Host%% floor for the log axis (default 0.02)"
    )
    parser.add_argument("--title", default="Whole-blood NGS chimerism vs flow lineages")
    parser.add_argument(
        "--anonymize",
        action="store_true",
        help="Replace sample names with S1, S2, ... on the x axis",
    )
    parser.add_argument(
        "--label-field",
        type=int,
        default=None,
        metavar="N",
        help="Shorten labels: split sample name on '_' and take this 0-based field",
    )
    parser.add_argument(
        "--label-code",
        action="store_true",
        help="Shorten labels to the patient code (last all-uppercase token); robust "
        "to the code sitting at field 3 or 4. Overrides --label-field.",
    )
    parser.add_argument(
        "--explicit-donor",
        default="",
        metavar="TOK1,TOK2",
        help="Comma-separated sample tokens (e.g. patient codes) that had an explicit "
        "donor genotype; their primary-run points are drawn as a star. Pair with a "
        "matching symbol in --title.",
    )
    parser.add_argument(
        "--sort",
        choices=("name", "tsv", "chimerism"),
        default="name",
        help="X-axis sample order: 'name' (alphabetical, default), 'tsv' (file "
        "order), or 'chimerism' (by measured donor fraction). Use 'name' or 'tsv' "
        "for a stable order that matches across runs.",
    )
    parser.add_argument(
        "--hide-lod",
        action="store_true",
        help="Do not shade the per-sample LOD band (drawn by default when the "
        "primary batch.tsv has an lod_pct column).",
    )
    parser.add_argument("--output", type=Path, default=Path("output/chimerism_comparison.png"))
    args = parser.parse_args()
    explicit_donor = {t.strip() for t in args.explicit_donor.split(",") if t.strip()}

    # Display order, left to right: compare runs as given, then the primary
    # (filled) on the right.
    primary = read_batch(args.batch_tsv, args.flow_column)
    compares = [read_batch(p, args.flow_column) for p in (args.compare_tsv or [])]
    runs = [*compares, primary]
    if args.labels is None:
        labels = [f"run{i + 1}" for i in range(len(runs))]
    elif len(args.labels) != len(runs):
        parser.error(
            f"--labels expects {len(runs)} label(s) (one per run, left to right), "
            f"got {len(args.labels)}"
        )
    else:
        labels = args.labels
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot(
        runs,
        labels,
        args.floor,
        args.title,
        args.output,
        args.anonymize,
        args.label_field,
        args.label_code,
        explicit_donor,
        args.sort,
        not args.hide_lod,
    )


if __name__ == "__main__":
    main()
