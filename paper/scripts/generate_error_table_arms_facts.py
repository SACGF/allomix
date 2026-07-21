"""Compare error-table arms on the public SRP434573 data (issue #49 follow-up).

The presence test's background model can come from three places, and this script
quantifies what each buys:

  flat         the built-in flat ``--error-rate`` default, no table
  per_mixture  one ``estimate-errors`` table per mixture, from that mixture's
               two reference individuals only
  pooled       one table pooled over all seven reference individuals

Arms are produced by ``paper/scripts/run_srp434573_allomix.py`` under env vars
(``ALLOMIX_NO_ERROR_TABLE``, ``ALLOMIX_FORCE_PER_MIXTURE_ERROR_TABLE``, and
``ALLOMIX_OUT_DIR`` to keep concurrent arms from writing the same paths), staged
into ``output/error_table_arms/<arm>/``. This script reads those and writes:

  output/facts/error_table_arms.csv     one row per arm (template variables)
  output/facts/fig_error_table_arms.png two panels, sensitivity and its cost

The two panels are deliberately paired. Panel A alone would be a sales pitch:
richer error tables raise low-fraction detection. Panel B is the price, the
false-signal floor on true-zero-recipient samples, shown against the
independently measured co-pooled contamination floor. Reporting A without B
would misrepresent the trade.

If an arm is missing (nobody has run it), it is skipped rather than failing, so
a fresh checkout still builds.
"""

import csv
import statistics as st
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

FRESH_ARMS_DIR = Path("output/error_table_arms")
SNAPSHOT_ARMS_DIR = Path("paper/public_data/SRP434573/error_table_arms")
BUILD_OUT = Path("output")
FACTS = Path("output/facts")
DETECT_ALPHA = 0.05


def arm_dir(arm: str) -> Path:
    """Where one arm's result TSVs live.

    A fresh three-arm run under output/error_table_arms wins for every arm.
    Otherwise:

    - ``pooled`` comes from the ordinary build output. The pipeline's default is
      the pooled table, so ``srp434573_allomix`` already writes exactly this arm;
      committing a second copy would duplicate an existing artifact and let the
      two drift apart.
    - ``flat`` and ``per_mixture`` come from the committed snapshot, because
      nothing else in the build produces them and regenerating costs ~30 min.
    """
    fresh = FRESH_ARMS_DIR / arm
    if fresh.is_dir():
        return fresh
    return BUILD_OUT if arm == "pooled" else SNAPSHOT_ARMS_DIR / arm

# Ordered worst -> best, so the arms read as a progression. A sequential ramp
# encodes that ordering; validated for CVD separation (normal-vision min dE 19.1,
# deuteranope/protanope min 17.9, both above the 15/8 floors). The obvious tab10
# categorical pick (blue/orange/green) fails protanopia at dE 1.4.
ARMS = [
    ("flat", "Flat default", "#9ecae1"),
    ("per_mixture", "Per-mixture table", "#4292c6"),
    ("pooled", "Pooled table", "#08519c"),
]


def _f(x: str | None) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _median(xs: list[float]) -> float | None:
    return st.median(xs) if xs else None


def arm_metrics(arm: str) -> dict | None:
    """Detection, accuracy and floor metrics for one arm, or None if absent."""
    d = arm_dir(arm)
    syn = _read(d / "srp434573_synthetic.tsv")
    two = _read(d / "srp434573_two_person.tsv")
    if not syn and not two:
        return None

    m: dict = {"arm": arm}

    # Semi-synthetic ladder: detection rate and signed offset per nominal level.
    fracs = sorted({_f(r["frac_pct"]) for r in syn if _f(r.get("frac_pct")) is not None})
    devs = []
    for fr in fracs:
        rows = [r for r in syn if _f(r["frac_pct"]) == fr]
        pv = [_f(r["presence_p"]) for r in rows if _f(r.get("presence_p")) is not None]
        mle = [_f(r["mle_pct"]) for r in rows if _f(r.get("mle_pct")) is not None]
        tag = f"{fr:g}".replace(".", "p")
        if pv:
            m[f"detect_{tag}"] = f"{sum(1 for p in pv if p < DETECT_ALPHA) / len(pv):.2f}"
        mm = _median(mle)
        if mm is not None:
            m[f"mle_med_{tag}"] = f"{mm:.3f}"
            devs.append(abs(mm - fr))
    if devs:
        m["ladder_mean_abs_dev_pct"] = f"{sum(devs) / len(devs):.3f}"

    # True-zero-recipient samples: the false-signal floor this arm carries.
    # These are the pure-donor endpoint rows, identified exactly as in
    # generate_srp434573_facts.zero_host_facts: blank known_pct, MLE recipient
    # below 50% (the pure-recipient endpoints sit near 100%).
    zero = [
        r for r in two
        if _f(r.get("known_pct")) is None and (_f(r.get("mle_pct")) or 0.0) < 50.0
    ]
    zmle = [_f(r["mle_pct"]) for r in zero if _f(r.get("mle_pct")) is not None]
    zp = [_f(r["presence_p"]) for r in zero if _f(r.get("presence_p")) is not None]
    if zmle:
        m["zero_n"] = str(len(zmle))
        m["zero_mle_med_pct"] = f"{_median(zmle):.3f}"
        m["zero_mle_max_pct"] = f"{max(zmle):.3f}"
    if zp:
        m["zero_presence_falsepos_n"] = str(sum(1 for p in zp if p < DETECT_ALPHA))

    # Real ladder, the lowest two rungs, where the error model actually bites.
    for known in (0.5, 1.0):
        rows = [r for r in two if _f(r.get("known_pct")) == known]
        mle = [_f(r["mle_pct"]) for r in rows if _f(r.get("mle_pct")) is not None]
        if mle:
            tag = f"{known:g}".replace(".", "p")
            m[f"real_{tag}_mle_med_pct"] = f"{_median(mle):.3f}"
    return m


def contamination_band() -> tuple[float, float] | None:
    """Independently measured co-pooled contamination floor (median, max).

    Read from the existing srp434573 facts. It is derived from the admix VCFs
    directly and does not depend on which error table was used, so it is the
    right reference line for panel B.
    """
    p = FACTS / "srp434573.csv"
    if not p.exists():
        return None
    with p.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return None
    med, mx = _f(rows[0].get("contam_floor_median_pct")), _f(rows[0].get("contam_floor_max_pct"))
    return (med, mx) if med is not None and mx is not None else None


def plot(metrics: dict[str, dict], out: Path) -> None:
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.6, 4.3))

    # Panel A: detection rate vs nominal recipient fraction.
    for arm, label, colour in ARMS:
        m = metrics.get(arm)
        if not m:
            continue
        xs, ys = [], []
        for k, v in m.items():
            if k.startswith("detect_"):
                xs.append(float(k[len("detect_") :].replace("p", ".")))
                ys.append(float(v))
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        axa.plot(
            [xs[i] for i in order], [ys[i] for i in order],
            marker="o", ms=6, lw=2, color=colour, label=label, zorder=3,
        )
    axa.set_xlabel("Nominal recipient fraction (%)")
    axa.set_ylabel("Presence-test detection rate")
    axa.set_ylim(-0.04, 1.06)
    axa.set_title("A. Low-fraction sensitivity", loc="left", fontsize=11)
    axa.grid(alpha=0.25, lw=0.6)
    axa.set_axisbelow(True)
    axa.legend(frameon=False, fontsize=9, loc="lower right")

    # Panel B: false-signal floor on true-zero-recipient samples.
    band = contamination_band()
    if band:
        axb.axhspan(band[0], band[1], color="0.85", zorder=0)
        axb.text(
            0.98, band[1], "measured contamination floor",
            transform=axb.get_yaxis_transform(), va="bottom", ha="right",
            fontsize=8, color="0.35",
        )
    labels, meds, maxs, colours = [], [], [], []
    for arm, label, colour in ARMS:
        m = metrics.get(arm)
        if not m or "zero_mle_med_pct" not in m:
            continue
        labels.append(label.replace(" ", "\n"))
        meds.append(float(m["zero_mle_med_pct"]))
        maxs.append(float(m["zero_mle_max_pct"]))
        colours.append(colour)
    xs = range(len(labels))
    axb.bar(xs, maxs, width=0.5, color=colours, zorder=2)
    axb.scatter(xs, meds, color="white", edgecolor="0.2", zorder=4, s=42)
    # Bars are the per-arm maximum, dots the median. Values are labelled directly
    # rather than legended, because a legend swatch would have to borrow one arm's
    # colour and so read as identity rather than as a statistic.
    for x, mx, md in zip(xs, maxs, meds, strict=True):
        axb.annotate(
            f"max {mx:.3f}\nmed {md:.3f}", (x, mx), textcoords="offset points",
            xytext=(0, 5), ha="center", fontsize=8, color="0.25",
        )
    axb.set_xticks(list(xs))
    axb.set_xticklabels(labels, fontsize=9)
    axb.set_ylabel("Estimated recipient % at true zero")
    axb.set_title("B. Cost: false-signal floor", loc="left", fontsize=11)
    axb.grid(alpha=0.25, lw=0.6, axis="y")
    axb.set_axisbelow(True)
    axb.set_ylim(0, max(maxs + [band[1] if band else 0]) * 1.35)

    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main() -> None:
    FACTS.mkdir(parents=True, exist_ok=True)
    metrics = {}
    for arm, _, _ in ARMS:
        m = arm_metrics(arm)
        # Say where each arm came from: pooled normally reads the ordinary build
        # output while the other two read the committed snapshot, so the sources
        # are deliberately not uniform.
        print(f"[sources] {arm}: {arm_dir(arm)}{'' if m else ' (absent, skipped)'}")
        if m:
            metrics[arm] = m
    if not metrics:
        print("No error-table arms found under output/error_table_arms/; skipping.")
        return

    keys: list[str] = []
    for m in metrics.values():
        for k in m:
            if k not in keys:
                keys.append(k)
    out_csv = FACTS / "error_table_arms.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for arm, _, _ in ARMS:
            if arm in metrics:
                w.writerow(metrics[arm])
    print(f"Wrote {out_csv} ({len(metrics)} arms)")

    # Flat single-row view for prose templating: the per-arm CSV above has one row
    # per arm, and the template namespace exposes only the first. Keys are
    # <arm>_<metric>, e.g. pooled_detect_0p2.
    headline: dict[str, str] = {}
    for arm, m in metrics.items():
        for k, v in m.items():
            if k != "arm":
                headline[f"{arm}_{k}"] = v
    band = contamination_band()
    if band:
        headline["contam_floor_median_pct"] = f"{band[0]:g}"
        headline["contam_floor_max_pct"] = f"{band[1]:g}"
    headline["n_arms"] = str(len(metrics))
    out_headline = FACTS / "error_table_arms_headline.csv"
    with out_headline.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(headline))
        w.writeheader()
        w.writerow(headline)
    print(f"Wrote {out_headline}")

    out_png = FACTS / "fig_error_table_arms.png"
    plot(metrics, out_png)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
