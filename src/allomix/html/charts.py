"""Chart builders for the HTML timeline report.

Each chart renders to an in-memory PNG and is returned as a base64 ``data:``
URI so the report stays a single self-contained file with no external image
references. Only the timeline report imports this module, so matplotlib is an
optional dependency (the ``report`` extra); the single-sample report does not
need it.

matplotlib is imported at module top level (no lazy import, per CLAUDE.md). The
CLI checks for matplotlib before dispatching the timeline-html path and emits a
clean install hint rather than an ImportError traceback.
"""

import base64
import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (must follow the Agg backend selection)

# Greyscale-safe, colour-blind-friendly series colours, reused in order. Paired
# with distinct markers and (for failed QC) hollow fills so the chart survives
# black-and-white printing, as the report styling does elsewhere.
_SERIES_COLOURS = ["#1f3a5f", "#8a5a00", "#1a7f37", "#b3261e"]
_SERIES_MARKERS = ["o", "s", "^", "D"]


def _png_data_uri(fig) -> str:
    """Serialise a matplotlib figure to a base64 PNG ``data:`` URI and close it.

    Args:
        fig: A matplotlib ``Figure``.

    Returns:
        A ``data:image/png;base64,...`` URI string.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def trend_png(
    x_labels: list[str],
    series: list[dict],
    *,
    log_scale: bool = False,
    ylabel: str = "Donor fraction (%)",
) -> str:
    """Render a chimerism trend chart as a base64 PNG ``data:`` URI.

    Draws one line per series with a shaded 95% CI band. Timepoints that failed
    QC are drawn as hollow markers so they stand out without relying on colour.

    Args:
        x_labels: Tick label per timepoint (collection date or sample name), in
            order. Defines the x axis.
        series: One dict per plotted line, each with keys ``name`` (legend
            label), ``y`` (values, percent), ``ci_lo``, ``ci_hi`` (CI bounds,
            percent; non-finite bounds drop that point's band), and ``qc_pass``
            (bool per timepoint). All lists align with ``x_labels``.
        log_scale: When True, use a logarithmic y axis (values must be > 0).
        ylabel: Y-axis label.

    Returns:
        A ``data:image/png;base64,...`` URI string.
    """
    x = list(range(len(x_labels)))
    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    for i, s in enumerate(series):
        colour = _SERIES_COLOURS[i % len(_SERIES_COLOURS)]
        marker = _SERIES_MARKERS[i % len(_SERIES_MARKERS)]
        y = s["y"]
        ax.plot(x, y, "-", color=colour, linewidth=1.5, label=s.get("name", f"series {i + 1}"))

        # CI band only across the span where both bounds are finite.
        lo = s.get("ci_lo") or []
        hi = s.get("ci_hi") or []
        band_x, band_lo, band_hi = [], [], []
        for xi, lvi, hvi in zip(x, lo, hi):
            if _is_finite(lvi) and _is_finite(hvi):
                band_x.append(xi)
                band_lo.append(lvi)
                band_hi.append(hvi)
        if len(band_x) >= 2:
            ax.fill_between(band_x, band_lo, band_hi, color=colour, alpha=0.15, linewidth=0)

        # Filled marker for a passing timepoint, hollow for a failed one.
        qc_pass = s.get("qc_pass") or [True] * len(y)
        for xi, yi, ok in zip(x, y, qc_pass):
            ax.plot(
                xi,
                yi,
                marker=marker,
                markersize=6,
                color=colour,
                markerfacecolor=(colour if ok else "white"),
                markeredgecolor=colour,
                markeredgewidth=1.4,
            )

    if log_scale:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=30, ha="right", fontsize=9)
    ax.grid(True, axis="y", linewidth=0.5, alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if len(series) > 1:
        ax.legend(fontsize=9, frameon=False)

    return _png_data_uri(fig)


def _is_finite(x: float | None) -> bool:
    """True when ``x`` is a real, finite number (not None, NaN, or inf)."""
    return x is not None and x == x and x not in (float("inf"), float("-inf"))
