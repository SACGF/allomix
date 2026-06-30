"""Formatting helpers for the HTML report.

Single source of significant-figure, NA, p-value, and CI consistency. Every
numeric formatter guards ``math.isfinite`` (and ``None``) because several result
fields default to ``float("inf")`` (``lob_fraction``, ``lod_fraction``, ``rho``)
and confidence intervals can be infinite, so an undetectable limit renders as a
stated NA glyph, never the string "inf".

These helpers produce small HTML-safe text fragments. Percentages, counts, and
p-values are plain ASCII digits, so they need no escaping. The status badge emits
a fixed span with no caller-supplied text. Sample-derived strings (names, IDs,
free-text warnings) are escaped at the section layer with ``html.escape``, not
here.
"""

import math

#: Sentinel shown wherever a value is unavailable. The em-dash GLYPH in rendered
#: output is fine: the no-em-dash project rule is about prose, not a data cell.
NA = "—"

# Status vocabulary. Each status carries a glyph so the colour coding survives
# greyscale printing and is legible to colour-blind readers (a requirement of
# the spec). Anything unrecognised falls back to a neutral badge.
_STATUS_GLYPH = {
    "PASS": "✓",  # check mark
    "REVIEW": "!",
    "FAIL": "✗",  # ballot X
}
_STATUS_CLASS = {
    "PASS": "badge-pass",
    "REVIEW": "badge-review",
    "FAIL": "badge-fail",
}


def _finite(x: float | None) -> bool:
    """True when ``x`` is a real, finite number (not None, NaN, or inf)."""
    return x is not None and math.isfinite(x)


def num(x: float | None, dp: int = 2) -> str:
    """Format a plain number to ``dp`` decimal places, or NA if not finite.

    Args:
        x: Value to format.
        dp: Decimal places.

    Returns:
        The fixed-point string, or the NA sentinel.
    """
    if not _finite(x):
        return NA
    return f"{x:.{dp}f}"


def count(n: int | None) -> str:
    """Format an integer count, or NA when None."""
    if n is None:
        return NA
    return str(int(n))


def pct(frac: float | None, dp: int = 2) -> str:
    """Format a fraction (0.0-1.0) as a percentage string, or NA.

    Args:
        frac: Fraction in [0, 1]; may be None or non-finite.
        dp: Decimal places on the percentage.

    Returns:
        For example ``"94.80%"``, or the NA sentinel when not finite.
    """
    if not _finite(frac):
        return NA
    return f"{frac * 100:.{dp}f}%"


def pct_points(frac: float | None, dp: int = 2) -> str:
    """Format a fraction as a bare percentage number (no ``%`` sign), or NA.

    Used in table cells where the column header already carries the unit.
    """
    if not _finite(frac):
        return NA
    return f"{frac * 100:.{dp}f}"


def ci(lo: float | None, hi: float | None, dp: int = 1, *, as_pct: bool = True) -> str:
    """Format a 95% confidence interval, or NA when either bound is not finite.

    Args:
        lo: Lower bound.
        hi: Upper bound.
        dp: Decimal places.
        as_pct: When True the bounds are fractions scaled to percent; when False
            they are shown as-is (already in display units).

    Returns:
        For example ``"(95% CI 94.0 to 95.6)"``, or the NA sentinel.
    """
    if not (_finite(lo) and _finite(hi)):
        return NA
    if as_pct:
        lo, hi = lo * 100, hi * 100
    return f"(95% CI {lo:.{dp}f} to {hi:.{dp}f})"


def pval(p: float | None) -> str:
    """Format a p-value in a consistent house style, or NA.

    Very small values collapse to ``"<0.001"``; everything else shows three
    decimal places. Keeps the report's p-values visually consistent regardless
    of which test produced them.

    Args:
        p: P-value, or None / non-finite.

    Returns:
        For example ``"<0.001"`` or ``"0.042"``, or the NA sentinel.
    """
    if not _finite(p):
        return NA
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def badge(status: str, *, large: bool = False) -> str:
    """Render a PASS / REVIEW / FAIL status badge as an HTML span.

    The badge pairs a glyph with the uppercase status text so it stays legible
    in greyscale and for colour-blind readers, not relying on colour alone.

    Args:
        status: One of "PASS", "REVIEW", "FAIL" (case-insensitive). An
            unrecognised value renders as a neutral badge with its own text.
        large: When True, add the ``badge-lg`` class for the headline verdict.

    Returns:
        An HTML ``<span>`` fragment.
    """
    key = status.upper()
    cls = _STATUS_CLASS.get(key, "badge-neutral")
    glyph = _STATUS_GLYPH.get(key, "•")  # bullet for unknown
    size = " badge-lg" if large else ""
    return f'<span class="badge {cls}{size}"><span class="badge-glyph">{glyph}</span> {key}</span>'
