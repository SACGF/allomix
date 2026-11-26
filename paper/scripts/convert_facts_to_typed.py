#!/usr/bin/env python3
"""One-off: rewrite existing one-row facts CSVs as typed vertical CSVs.

The paper's producers now write typed vertical facts (``facts_io``). This
converts a facts directory produced before that change, so a slow build's
results do not have to be re-run to pick up the types. Values are unchanged:
each file goes through the same ``write_fact_csv`` the producers use. Already
vertical files are left alone, and no provenance sidecar is written, since the
run that produced the values is not this one. Multi-row data tables keep their
contents and only have CRLF line endings rewritten as LF, matching what the
producers now write.

Usage:
    python paper/scripts/convert_facts_to_typed.py [output/facts]
"""

import csv
import sys
from pathlib import Path

from facts_io import read_fact_csv, write_fact_csv


def main() -> int:
    facts_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "output/facts")
    n_converted = n_relined = 0
    for path in sorted(facts_dir.glob("*.csv")):
        raw = path.read_bytes()
        rows = list(csv.reader(raw.decode("utf-8").splitlines()))
        if not rows or rows[0][:2] == ["field", "value"]:
            continue
        if len(rows) == 2:
            write_fact_csv(path, read_fact_csv(path), provenance=False)
            n_converted += 1
            print(f"converted {path}")
        elif b"\r\n" in raw:
            path.write_bytes(raw.replace(b"\r\n", b"\n"))
            n_relined += 1
    print(f"{n_converted} facts CSV(s) converted, {n_relined} data table(s) relined in {facts_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
