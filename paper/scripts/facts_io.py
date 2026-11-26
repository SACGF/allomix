"""Read and write the paper's single-row facts CSVs (vibepaper typed facts).

Every facts file the paper's prose reads goes through :func:`write_fact_csv`,
which writes a vertical ``field,value,type`` CSV with vibepaper's
``write_facts`` and a provenance sidecar (``<name>.provenance.json``) naming the
allomix version that produced it.

Types come from the paper's naming convention, so producers do not type each
field by hand:

- ``*_pct`` (and ``est_<N>pct``): ``percent``, already multiplied by 100.
- ``n_*`` and ``*_n``: ``count``, a whole number of 0 or more.
- ``*_rate``, ``detect_rate_*`` and ``*_detect_0p<N>``: ``fraction`` (0 to 1).
- ``*_ratio``, ``*_fold`` and ``fold_*``: ``ratio`` (positive).

Anything else is written untyped. So is a value that is not a finite number
(a quick-build ``nan`` LoD cell, or a ``">20"`` beyond-ceiling marker), which
renders exactly as it did before. A value that has a type but is out of its
range (a negative count, a fraction above 1) raises, because that is a bug in
the producer, not something to print.

Multi-row data tables (per-pair, per-site, summary grids) are not facts files
and keep their own ``csv`` writers.
"""

import csv
import math
import re
from pathlib import Path

from vibepaper.facts import Count, Fraction, Percent, Ratio, write_facts
from vibepaper.provenance import write_provenance

_RULES = (
    (re.compile(r".*_pct$|^est_\d+pct$"), Percent),
    (re.compile(r"^n_.*|.*_n$"), Count),
    (re.compile(r".*_rate$|^detect_rate_.*|.*_detect_0p\d+$"), Fraction),
    (re.compile(r".*_ratio$|.*_fold$|^fold_.*"), Ratio),
)


def fact_type_for(name: str) -> type | None:
    """The vibepaper type class the naming convention gives ``name``, or None."""
    for pattern, cls in _RULES:
        if pattern.fullmatch(name):
            return cls
    return None


def _plain(value: object) -> object:
    """A numpy scalar as the equivalent Python value; anything else unchanged."""
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        return value.item()
    return value


def _finite_number(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def typed_facts(facts: dict) -> dict:
    """``facts`` with each value wrapped in its convention type where it applies."""
    out = {}
    for name, value in facts.items():
        value = _plain(value)
        cls = fact_type_for(name)
        number = _finite_number(value)
        if cls is not None and number is not None:
            out[name] = cls(value if isinstance(value, int) else number)
        elif isinstance(value, float) and not math.isfinite(value):
            # write_facts refuses a float NaN; keep the text it rendered as before.
            out[name] = str(value)
        else:
            out[name] = value
    return out


def write_fact_csv(path: Path | str, facts: dict, provenance: bool = True) -> Path:
    """Write ``facts`` as a typed vertical facts CSV.

    Args:
        path: Destination ``.csv``.
        facts: Field name to value, in output order.
        provenance: Also write the ``.provenance.json`` sidecar (the allomix
            version, today's date, the producing command). Off only when
            rewriting a file whose producing run is not this one.
    """
    path = write_facts(Path(path), typed_facts(facts))
    if provenance:
        write_provenance(path, tools=["allomix"])
    return path


def read_fact_csv(path: Path | str) -> dict[str, str]:
    """One facts CSV as ``{field: raw string}``, vertical or legacy horizontal."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    if rows and rows[0][:2] == ["field", "value"]:
        return {row[0]: row[1] for row in rows[1:] if row}
    if len(rows) != 2:
        raise ValueError(f"{path}: expected a vertical or a one-row facts CSV")
    return dict(zip(rows[0], rows[1], strict=True))
