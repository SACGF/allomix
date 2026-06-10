"""Plot allomix accuracy on the SRP434573 real-data mixture dilution series.

Reads the two TSVs written by run_srp434573_allomix.py and produces:

  output/srp434573_scatter.png       known vs allomix estimated host %, linear
                                     scatter with the perfect-recovery line. Both
                                     allomix outputs are shown: the MLE host
                                     fraction (100 - donor%) and the presence-test
                                     host fraction.
  output/srp434573_logy.png          log-Y host% dot plot. No joining lines:
                                     each timepoint shows both allomix estimates
                                     (MLE = filled circle, presence = open square,
                                     each + CI) with the known fraction as a grey
                                     diamond beside them, grouped by mixture.
  output/srp434573_three_person.png  the single three-person mixture (1:3:5 of
                                     F2:M1:M2) as grouped known-vs-estimated
                                     bars, kept separate from the accuracy series.

House style: log axis with percent-formatted ticks, QC-REVIEW points circled in
red, dpi 150.
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

OUT = Path("output")
YFLOOR = 0.012  # where "not detected" (est 0%) points are drawn on the log axis
SHOW_QC_REVIEW = False  # circle QC=REVIEW points in red and add a legend entry


def _fmt_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _read(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _label(mixture: str) -> str:
    return mixture.replace("mix_", "").replace("_into_", "→")


def _grouped(rows: list[dict]) -> tuple[list[str], dict[str, list[dict]]]:
    """Group two-person rows by mixture, each sorted by known fraction (desc)."""
    by_mix: dict[str, list[dict]] = {}
    for r in rows:
        by_mix.setdefault(r["mixture"], []).append(r)
    for v in by_mix.values():
        v.sort(key=lambda r: -(_f(r["known_pct"]) or 0.0))
    return sorted(by_mix), by_mix


def _colors(mixes: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {m: cmap(i % 10) for i, m in enumerate(mixes)}


MLE_COLOR = "#1f77b4"
PRESENCE_COLOR = "#d62728"


def plot_scatter(rows: list[dict], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    lims = [0, 11.5]
    ax.plot(lims, lims, color="0.4", linestyle="--", linewidth=1.2, zorder=1,
            label="perfect recovery (y = x)")

    kx = [_f(r["known_pct"]) for r in rows]
    mle = [(_f(r["mle_pct"]) or 0.0) for r in rows]
    pres = [(_f(r.get("presence_pct")) or 0.0) for r in rows]
    ax.scatter(kx, mle, s=55, color=MLE_COLOR, edgecolor="white", linewidth=0.6,
               marker="o", zorder=3, label="MLE (100 − donor%)")
    ax.scatter(kx, pres, s=55, facecolor="none", edgecolor=PRESENCE_COLOR,
               linewidth=1.4, marker="s", zorder=3, label="presence-test")

    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_aspect("equal")
    ticks = [0, 0.5, 1, 2.5, 5, 10]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("Known host fraction", fontsize=12)
    ax.set_ylabel("allomix estimated host %", fontsize=12)
    ax.set_title("allomix vs known host fraction (SRP434573)",
                 fontsize=12.5, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.92)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_logy(rows: list[dict], out_path: Path) -> None:
    mixes, by_mix = _grouped(rows)
    colors = _colors(mixes)

    # Lay samples out left to right, grouped by mixture with a gap between groups.
    x = 0.0
    group_span: dict[str, tuple[float, float]] = {}
    any_nd = any_review = False
    fig, ax = plt.subplots(figsize=(18, 6.6))

    for gi, m in enumerate(mixes):
        start = x
        if gi % 2 == 1:
            pass  # shading added after we know the span
        for r in by_mix[m]:
            known = _f(r["known_pct"])
            est = _f(r["mle_pct"])
            lo, hi = _f(r["mle_ci_lo"]), _f(r["mle_ci_hi"])
            pest = _f(r.get("presence_pct"))
            plo, phi = _f(r.get("presence_ci_lo")), _f(r.get("presence_ci_hi"))
            # Layout per timepoint: MLE (left), known diamond (centre), presence (right).
            x_mle, x_known, x_pres = x - 0.20, x, x + 0.20
            if known is not None:
                ax.plot(x_known, known, marker="D", markersize=6, color="0.45",
                        zorder=2)
            # MLE: filled circle + CI. Not detected (est 0) -> the same circle at
            # the log-axis floor (a stand-in 0), so the shape still reads as MLE.
            if est is not None and est > 0:
                ax.errorbar(x_mle, est,
                            yerr=[[est - lo if lo is not None else 0.0],
                                  [hi - est if hi is not None else 0.0]],
                            marker="o", markersize=7.5, color=colors[m], capsize=2.5,
                            linewidth=0, elinewidth=1.3, markeredgecolor="white",
                            markeredgewidth=0.6, zorder=3)
                y_mle = est
            else:
                ax.plot(x_mle, YFLOOR, marker="o", markersize=7.5, color=colors[m],
                        markeredgecolor="white", markeredgewidth=0.6, zorder=3)
                y_mle = YFLOOR
                any_nd = True
            if SHOW_QC_REVIEW and r.get("qc") == "REVIEW":
                ax.plot(x_mle, y_mle, marker="o", markersize=12, markerfacecolor="none",
                        markeredgecolor="red", markeredgewidth=1.2, zorder=4)
                any_review = True
            # Presence: open square + CI. Not detected -> the same square at the floor.
            if pest is not None and pest > 0:
                ax.errorbar(x_pres, pest,
                            yerr=[[pest - plo if plo is not None else 0.0],
                                  [phi - pest if phi is not None else 0.0]],
                            marker="s", markersize=7, markerfacecolor="none",
                            markeredgecolor=colors[m], markeredgewidth=1.5,
                            ecolor=colors[m], capsize=2.5, linewidth=0,
                            elinewidth=1.1, zorder=3)
            else:
                ax.plot(x_pres, YFLOOR, marker="s", markersize=7, markerfacecolor="none",
                        markeredgecolor=colors[m], markeredgewidth=1.5, zorder=3)
                any_nd = True
            x += 1.4
        group_span[m] = (start, x - 1.4)
        x += 1.2  # gap between mixtures

    # Alternating background shading + group labels.
    for gi, m in enumerate(mixes):
        lo_x, hi_x = group_span[m]
        if gi % 2 == 1:
            ax.axvspan(lo_x - 0.6, hi_x + 0.85, color="0.92", zorder=0)
        ax.text((lo_x + hi_x) / 2, 0.0145, _label(m), ha="center", va="bottom",
                fontsize=9, rotation=0, transform=ax.get_xaxis_transform(),
                fontweight="bold", color=colors[m])

    ax.set_yscale("log")
    ax.set_ylim(0.008, 13)
    ax.set_xlim(-1.0, x - 0.5)
    ax.set_xticks([])
    # The floor row is a stand-in for 0 (true 0 has no place on a log axis).
    ax.axhline(YFLOOR, color="0.6", linewidth=0.8, linestyle=":", zorder=0)
    ax.yaxis.set_major_locator(FixedLocator([YFLOOR, 0.1, 0.5, 1, 2.5, 5, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, p: "0" if abs(v - YFLOOR) < 1e-6 else _fmt_pct(v, p))
    )
    ax.grid(True, axis="y", which="major", alpha=0.25)
    ax.set_ylabel("host % (log scale)", fontsize=12)
    ax.set_title(
        "allomix host-fraction estimates: MLE (filled circle, 100 − donor%) and "
        "presence-test (open square) vs known (grey diamond), per timepoint",
        fontsize=12, fontweight="bold",
    )

    handles = [
        Line2D([0], [0], marker="o", color="0.3", linestyle="none", markersize=8,
               label="MLE fraction + 95% CI"),
        Line2D([0], [0], marker="s", color="none", linestyle="none", markersize=8,
               markerfacecolor="none", markeredgecolor="0.3", markeredgewidth=1.5,
               label="presence-test fraction + 95% CI"),
        Line2D([0], [0], marker="D", color="0.45", linestyle="none", markersize=7,
               label="known fraction"),
    ]
    if any_review:
        handles.append(Line2D([0], [0], marker="o", color="none", markersize=11,
                              markerfacecolor="none", markeredgecolor="red",
                              markeredgewidth=1.3, label="QC = REVIEW"))
    ax.legend(handles=handles, fontsize=9, loc="upper left",
              bbox_to_anchor=(1.005, 1.0), framealpha=0.95)
    if any_nd:
        ax.text(1.005, 0.0, 'a marker on the "0" row\n= not detected', fontsize=8,
                transform=ax.transAxes, va="bottom", ha="left", color="0.35")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_three_person(rows: list[dict], out_path: Path) -> None:
    order = [("F2", "host"), ("M1", "donor"), ("M2", "donor")]
    by_comp = {r["component"]: r for r in rows}
    labels, known, est, elo, ehi = [], [], [], [], []
    for comp, role in order:
        r = by_comp.get(comp)
        if not r:
            continue
        labels.append(f"{comp}\n({role})")
        known.append(_f(r["known_pct"]))
        e = _f(r["est_pct"])
        est.append(e)
        lo, hi = _f(r["ci_lo"]), _f(r["ci_hi"])
        elo.append((e - lo) if (e is not None and lo is not None) else 0.0)
        ehi.append((hi - e) if (e is not None and hi is not None) else 0.0)

    x = range(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    ax.bar([i - w / 2 for i in x], known, width=w, color="0.7", label="known")
    ax.bar([i + w / 2 for i in x], est, width=w, color="#2c7fb8",
           yerr=[elo, ehi], capsize=4, label="allomix")
    for i, (k, e) in enumerate(zip(known, est)):
        ax.text(i - w / 2, k + 1.2, f"{k:.1f}%", ha="center", fontsize=9, color="0.3")
        if e is not None:
            ax.text(i + w / 2, e + 1.2, f"{e:.1f}%", ha="center", fontsize=9,
                    color="#2c7fb8", fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("fraction of sample (%)", fontsize=11)
    ax.set_ylim(0, max(known + [e for e in est if e is not None]) * 1.18)
    ax.set_title("Three-person mixture (1:3:5 of F2:M1:M2), single sample",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> int:
    two = _read(OUT / "srp434573_two_person.tsv")
    plot_scatter(two, OUT / "srp434573_scatter.png")
    plot_logy(two, OUT / "srp434573_logy.png")
    plot_three_person(_read(OUT / "srp434573_three_person.tsv"),
                      OUT / "srp434573_three_person.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
