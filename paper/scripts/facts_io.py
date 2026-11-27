"""The paper's facts convention: one vibepaper ``FactsWriter`` shared by every producer.

Every facts file the paper's prose reads goes through ``FACTS.write`` (alias
:func:`write_fact_csv`), which writes a vertical ``field,value,type`` CSV and a
provenance sidecar (``<name>.provenance.json``) naming the allomix version that
produced it.

Types come from the paper's naming convention, so producers do not type each
field by hand (first match wins):

- ``*_pct`` and ``est_*pct``: ``percent``, already multiplied by 100.
- ``n_*`` and ``*_n``: ``count``, a whole number of 0 or more.
- ``*_rate``, ``detect_rate_*`` and ``*_detect_0p*``: ``fraction`` (0 to 1).
- ``*_ratio``, ``*_fold`` and ``fold_*``: ``ratio`` (positive).

Anything else is written untyped, and so is a value that is not a finite
number (a ``">20"`` beyond-ceiling marker). A missing value (None or NaN, e.g.
a quick-build LoD cell) is written as an empty cell, which loads as NaN. A
value that has a type but is out of its range (a negative count, a fraction
above 1) raises, because that is a bug in the producer, not something to print.

Read a facts CSV back with ``vibepaper.facts.read_facts``, which returns the
values as a template sees them.

Multi-row data tables (per-pair, per-site, summary grids) are not facts files
and keep their own ``csv`` writers.
"""

from vibepaper.facts import Count, FactsWriter, Fraction, Percent, Ratio

FACTS = FactsWriter(
    {
        "*_pct": Percent,
        "est_*pct": Percent,
        "n_*": Count,
        "*_n": Count,
        "*_rate": Fraction,
        "detect_rate_*": Fraction,
        "*_detect_0p*": Fraction,
        "*_ratio": Ratio,
        "*_fold": Ratio,
        "fold_*": Ratio,
    },
    missing="blank",
    tools=["allomix"],
)

write_fact_csv = FACTS.write
