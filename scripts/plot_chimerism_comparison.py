#!/usr/bin/env python3
"""Plot allomix whole-blood chimerism against flow-sorted lineage values.

Internal SA Path validation plot (not part of the allomix package).

Takes one or two allomix `batch.tsv` files (as produced by
`scripts/run_xls_batch.py`) and draws a per-sample forest plot. Each sample
shows the NGS donor estimate with its confidence interval, overlaid on the
flow cytometry lineage values (CD45 / CD3 / CD13) parsed from a copied column.

The y axis shows donor %, but is log-spaced by distance from 100% donor and
inverted so 100% is at the top. The clinically interesting action is the
low-level signal near full donor chimerism, which a plain linear (or plain log)
donor axis compresses; this spacing keeps it readable. The spread of the
CD3/CD13 subsets shows where the true whole-blood value can sit depending on the
cell differential. QC-FAIL samples are skipped.

Convention: all percentages are % DONOR (allomix `donor_pct` and the flow
values).

Usage:
    # Single run (run2), flow column already merged into the batch.tsv:
    python scripts/plot_chimerism_comparison.py \
        output/validation_run2/batch.tsv \
        --flow-column "Chimerism result TP2" \
        --output output/chimerism_comparison.png

    # Compare two runs (e.g. run1 SID/tp0-donor vs run2 full-panel/explicit):
    python scripts/plot_chimerism_comparison.py \
        output/validation_run2/batch.tsv \
        --compare-tsv output/validation_run1/batch.tsv \
        --labels "run2 (full panel)" "run1 (SID only)" \
        --output output/run1_vs_run2.png
"""

from __future__ import annotations

import argparse
import csv
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
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
            rec: dict = {
                "donor_pct": float(row["donor_pct"]),
                "ci_lo": float(row["ci_lo"]),
                "ci_hi": float(row["ci_hi"]),
                "n_informative": int(row.get("n_informative", 0) or 0),
                "mean_depth": float(row.get("mean_depth", 0) or 0),
                "qc_pass": row.get("qc_pass", "PASS") == "PASS",
                "donor": (row.get(donor_column) or "").strip(),
                "lineages": {},
            }
            # lob_pct / lod_pct are present only on batch.tsv files produced
            # after the per-sample LOD change; tolerate their absence.
            if row.get("lod_pct") not in (None, "", "NA"):
                rec["lod_pct"] = float(row["lod_pct"])
            if flow_column and row.get(flow_column):
                rec["lineages"] = parse_lineages(row[flow_column])
            out[sample] = rec
    return out


def host(donor_pct: float) -> float:
    """Convert a donor percentage to host percentage."""
    return 100.0 - donor_pct


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


def plot(
    primary: dict[str, dict],
    compare: dict[str, dict] | None,
    labels: tuple[str, str],
    floor: float,
    title: str,
    output: Path,
    anonymize: bool,
    label_field: int | None,
    label_code: bool,
    explicit_donor: set[str],
) -> None:
    """Draw the forest plot and write it to ``output``.

    The y axis shows donor %, but is log-spaced by distance from 100% donor
    (internally 100 - donor%) and inverted so 100% sits at the top. A plain log
    of donor% would compress everything near 100% where the low-level signal
    lives; this keeps that region readable while labelling the axis in donor %.

    Each x-axis label carries the informative-marker count(s) so the panel
    difference between runs is visible at a glance.
    """
    # QC-FAIL samples (e.g. a single informative marker) have point estimates
    # too unreliable to plot, so drop them rather than imply a real measurement.
    ordered = sorted(primary, key=lambda s: host(primary[s]["donor_pct"]))
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

    ngs_color = "#1f77b4"
    cmp_color = "#ff7f0e"
    lineage_color = "#888888"

    # Dodge: time runs left to right, so run1 (compare) sits left of run2
    # (primary). In single-run mode the one estimate is centred.
    flow_dx = -0.34 if compare is not None else -0.22
    run2_dx = 0.12 if compare is not None else 0.0

    for x, sample in enumerate(samples):
        rec = primary[sample]

        # Flow lineage markers, offset left of centre, each at its own donor
        # value (100% sits at the top). A thin line shows the lineage spread.
        lin = rec["lineages"]
        if lin:
            lx = x + flow_dx
            host_vals = [clamp(host(v)) for v in lin.values()]
            if max(host_vals) > min(host_vals):
                ax.plot(
                    [lx, lx],
                    [min(host_vals), max(host_vals)],
                    color=lineage_color,
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
                    edgecolors=lineage_color,
                    zorder=3,
                )

        # NGS estimate(s) with CI error bars. A compare run that failed or found
        # no informative markers has no point to plot (its 0 still shows in the
        # x-axis label).
        _ngs_point(ax, x + run2_dx, rec, ngs_color, clamp)
        if compare is not None and sample in compare:
            crec = compare[sample]
            if crec["qc_pass"] and crec["n_informative"] > 0:
                _ngs_point(ax, x - 0.12, crec, cmp_color, clamp, hollow=True)

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

    # Per-point x-axis labels, colour-matched to the runs. In compare mode run1
    # is the top row and run2 a separate row below it (time runs left to right).
    ax.set_xticks([])
    trans = ax.get_xaxis_transform()
    # Row y-positions below the axis: run label(s), then donor type, then key.
    donor_y = -0.075 if compare is None else -0.135
    xlabel_y = -0.20 if compare is None else -0.27
    for i, s in enumerate(samples):
        code = code_for(i, s)
        p = primary[s]
        # Flag explicit-donor samples with a star on the sample name (keyed in
        # title). It marks the primary run, so only that row's code is starred.
        mark = "★" if any(tok and tok in s for tok in explicit_donor) else ""
        # Runs are distinguished by colour (see legend), so the run name is
        # left off the labels.
        if compare is None:
            ax.text(
                i,
                -0.02,
                f"{code}{mark} M:{p['n_informative']} D:{p['mean_depth']:.0f}x",
                transform=trans,
                ha="center",
                va="top",
                fontsize=8,
            )
        else:
            c = compare.get(s)
            run1_txt = f"M:{c['n_informative']} D:{c['mean_depth']:.0f}x" if c else "NA"
            ax.text(
                i,
                -0.02,
                f"{code} {run1_txt}",
                transform=trans,
                ha="center",
                va="top",
                fontsize=8,
                color=cmp_color,
            )
            ax.text(
                i,
                -0.07,
                f"{code}{mark} M:{p['n_informative']} D:{p['mean_depth']:.0f}x",
                transform=trans,
                ha="center",
                va="top",
                fontsize=8,
                color=ngs_color,
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
    xlabel = "Sample: M = informative markers, D = mean depth"
    ax.text((n - 1) / 2.0, xlabel_y, xlabel, transform=trans, ha="center", va="top", fontsize=9)

    # Legend order matches the left-to-right draw order: flow, then run1, run2.
    handles = [
        Line2D([], [], color=lineage_color, marker="s", mfc="none", lw=0, label="Flow CD45"),
        Line2D([], [], color=lineage_color, marker="^", mfc="none", lw=0, label="Flow CD3"),
        Line2D([], [], color=lineage_color, marker="o", mfc="none", lw=0, label="Flow CD13"),
    ]
    if compare is not None:
        handles.append(
            Line2D(
                [], [], color=cmp_color, marker="o", mfc="none", lw=0, label=f"NGS {labels[1]} (CI)"
            )
        )
    handles.append(Line2D([], [], color=ngs_color, marker="o", lw=0, label=f"NGS {labels[0]} (CI)"))
    ax.legend(handles=handles, fontsize=8, loc="lower left", framealpha=0.9)

    fig.tight_layout()
    # Reserve room for the custom label rows + donor type + key (tight_layout
    # ignores them).
    fig.subplots_adjust(bottom=0.26 if compare is None else 0.33)
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("batch_tsv", type=Path, help="Primary batch.tsv (e.g. run2)")
    parser.add_argument(
        "--compare-tsv",
        type=Path,
        default=None,
        help="Optional second batch.tsv to overlay (e.g. run1)",
    )
    parser.add_argument(
        "--flow-column",
        default="Chimerism result TP2",
        help="Column name holding flow lineage strings",
    )
    parser.add_argument(
        "--labels",
        nargs=2,
        default=["run2", "run1"],
        metavar=("PRIMARY", "COMPARE"),
        help="Legend labels for the two runs",
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
    parser.add_argument("--output", type=Path, default=Path("output/chimerism_comparison.png"))
    args = parser.parse_args()
    explicit_donor = {t.strip() for t in args.explicit_donor.split(",") if t.strip()}

    primary = read_batch(args.batch_tsv, args.flow_column)
    compare = read_batch(args.compare_tsv, args.flow_column) if args.compare_tsv else None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plot(
        primary,
        compare,
        tuple(args.labels),
        args.floor,
        args.title,
        args.output,
        args.anonymize,
        args.label_field,
        args.label_code,
        explicit_donor,
    )


if __name__ == "__main__":
    main()
