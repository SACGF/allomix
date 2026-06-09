"""Plot allomix accuracy on the SRP434573 real-data mixture dilution series.

Reads the two TSVs written by run_srp434573_allomix.py and produces:

  output/srp434573_scatter.png       known vs allomix estimated donor %, linear
                                     scatter with the perfect-recovery line. One
                                     dot per two-person timepoint, coloured by
                                     mixture.
  output/srp434573_logy.png          log-Y donor% dot plot. No joining lines:
                                     each timepoint shows the allomix estimate
                                     (dot + CI) with the known fraction as a grey
                                     diamond right beside it, grouped by mixture.
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


def plot_scatter(rows: list[dict], out_path: Path) -> None:
    mixes, by_mix = _grouped(rows)
    colors = _colors(mixes)

    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    lims = [0, 11.5]
    ax.plot(lims, lims, color="0.4", linestyle="--", linewidth=1.2, zorder=1,
            label="perfect recovery (y = x)")

    for m in mixes:
        xs = [_f(r["known_pct"]) for r in by_mix[m]]
        ys = [(_f(r["est_pct"]) or 0.0) for r in by_mix[m]]
        ax.scatter(xs, ys, s=55, color=colors[m], edgecolor="white",
                   linewidth=0.6, zorder=3, label=_label(m))

    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_aspect("equal")
    ticks = [0, 0.5, 1, 2.5, 5, 10]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("Known donor fraction", fontsize=12)
    ax.set_ylabel("allomix estimated donor %", fontsize=12)
    ax.set_title("allomix vs known donor fraction (SRP434573)",
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
    fig, ax = plt.subplots(figsize=(15, 6.4))

    for gi, m in enumerate(mixes):
        start = x
        if gi % 2 == 1:
            pass  # shading added after we know the span
        for r in by_mix[m]:
            known = _f(r["known_pct"])
            est = _f(r["est_pct"])
            lo, hi = _f(r["ci_lo"]), _f(r["ci_hi"])
            # known reference: grey diamond just to the right of the estimate
            if known is not None:
                ax.plot(x + 0.28, known, marker="D", markersize=6, color="0.45",
                        zorder=2)
            if est is None or est <= 0:
                ax.plot(x, YFLOOR, marker="v", markersize=10, markerfacecolor="none",
                        markeredgecolor=colors[m], markeredgewidth=1.8, zorder=3)
                any_nd = True
            else:
                yerr = [[est - lo if lo is not None else 0.0],
                        [hi - est if hi is not None else 0.0]]
                ax.errorbar(x, est, yerr=yerr, marker="o", markersize=8,
                            color=colors[m], capsize=3, linewidth=0,
                            elinewidth=1.4, markeredgecolor="white",
                            markeredgewidth=0.6, zorder=3)
                if r.get("qc") == "REVIEW":
                    ax.plot(x, est, marker="o", markersize=13, markerfacecolor="none",
                            markeredgecolor="red", markeredgewidth=1.3, zorder=4)
                    any_review = True
            x += 1.0
        group_span[m] = (start, x - 1.0)
        x += 1.0  # gap between mixtures

    # Alternating background shading + group labels.
    for gi, m in enumerate(mixes):
        lo_x, hi_x = group_span[m]
        if gi % 2 == 1:
            ax.axvspan(lo_x - 0.5, hi_x + 0.78, color="0.92", zorder=0)
        ax.text((lo_x + hi_x) / 2, 0.0145, _label(m), ha="center", va="bottom",
                fontsize=9, rotation=0, transform=ax.get_xaxis_transform(),
                fontweight="bold", color=colors[m])

    ax.set_yscale("log")
    ax.set_ylim(0.008, 13)
    ax.set_xlim(-1.0, x - 0.5)
    ax.set_xticks([])
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 0.5, 1, 2.5, 5, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.grid(True, axis="y", which="major", alpha=0.25)
    ax.set_ylabel("donor % (log scale)", fontsize=12)
    ax.set_title(
        "allomix estimate (coloured dot + 95% CI) vs known fraction "
        "(grey diamond), per timepoint",
        fontsize=12.5, fontweight="bold",
    )

    handles = [
        Line2D([0], [0], marker="o", color="0.3", linestyle="none", markersize=8,
               label="allomix estimate + 95% CI"),
        Line2D([0], [0], marker="D", color="0.45", linestyle="none", markersize=7,
               label="known fraction"),
    ]
    if any_review:
        handles.append(Line2D([0], [0], marker="o", color="none", markersize=11,
                              markerfacecolor="none", markeredgecolor="red",
                              markeredgewidth=1.3, label="QC = REVIEW"))
    if any_nd:
        handles.append(Line2D([0], [0], marker="v", color="none", markersize=10,
                              markerfacecolor="none", markeredgecolor="0.4",
                              markeredgewidth=1.8, label="not detected (est 0%)"))
    ax.legend(handles=handles, fontsize=9, loc="upper left",
              bbox_to_anchor=(1.005, 1.0), framealpha=0.95)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_three_person(rows: list[dict], out_path: Path) -> None:
    order = [("F2", "donor"), ("M1", "donor"), ("M2", "host")]
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
