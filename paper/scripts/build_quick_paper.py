"""Render a draft of the paper from a quick build's facts (preview only).

A ``quick=1`` build produces every facts file, but its low-iteration numbers
fail some of the claims asserts, and a few LoD cells are never reached and are
written blank (NaN). ``vibepaper build`` stops on both, so the quick build could
not render the paper. This script renders it anyway, so a quick build checks
that every section, figure, table and template expression still renders:

- It copies ``output/facts`` to ``output/quick_paper/facts`` and fills each
  blank fact with a type-valid placeholder (0 for a percent, fraction, count or
  untyped value; 1 for a ratio), listing each one it fills.
- It reports failed claims asserts as warnings instead of stopping the build.

The output goes to ``output/quick_paper/`` with ``_QUICK`` in the file name, and
the vibepaper build cache goes to ``output/quick_paper/build`` so it does not
replace the cache that ``vibepaper diff`` compares the real build against. The
numbers in the draft are not results: use the full build for anything that is
read for its content.

Usage (from the repo root, after a quick build has written ``output/facts``)::

    python paper/scripts/build_quick_paper.py
"""

import csv
import shutil
import sys
from pathlib import Path

import vibepaper.build
from facts_io import FACTS
from vibepaper.facts import Count, Fraction, Percent, Ratio

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FACTS_DIR = PROJECT_ROOT / "output" / "facts"
OUT_DIR = PROJECT_ROOT / "output" / "quick_paper"
QUICK_FACTS_DIR = OUT_DIR / "facts"

# A blank fact's placeholder, by the type its field name implies.
PLACEHOLDERS = {
    Percent: ("0", "percent"),
    Fraction: ("0", "fraction"),
    Count: ("0", "count"),
    Ratio: ("1", "ratio"),
}
UNTYPED_PLACEHOLDER = ("0", "")


def fill_blank_facts(facts_dir: Path) -> list[str]:
    """Fill each blank value in the vertical facts CSVs in ``facts_dir``, in place.

    Args:
        facts_dir: Directory of facts CSVs (a copy, since the files are rewritten).

    Returns:
        One ``"<file>: <field> = <placeholder>"`` line per filled value.
    """
    filled = []
    for path in sorted(facts_dir.glob("*.csv")):
        with path.open(newline="") as fh:
            rows = list(csv.reader(fh))
        if not rows or rows[0][:2] != ["field", "value"]:
            continue  # a horizontal data table, not a facts file
        has_type = len(rows[0]) > 2 and rows[0][2] == "type"
        changed = False
        for row in rows[1:]:
            if row[1].strip():
                continue
            value, type_name = PLACEHOLDERS.get(FACTS.type_for(row[0]), UNTYPED_PLACEHOLDER)
            row[1] = value
            if has_type:
                row[2] = type_name
            filled.append(f"{path.stem}: {row[0]} = {value}")
            changed = True
        if changed:
            with path.open("w", newline="") as fh:
                csv.writer(fh).writerows(rows)
    return filled


def main() -> int:
    """Copy and fill the facts, then render the paper with asserts as warnings."""
    if not FACTS_DIR.is_dir():
        sys.stderr.write(f"{FACTS_DIR} not found: run the quick build first.\n")
        return 1
    if QUICK_FACTS_DIR.exists():
        shutil.rmtree(QUICK_FACTS_DIR)
    shutil.copytree(FACTS_DIR, QUICK_FACTS_DIR)

    filled = fill_blank_facts(QUICK_FACTS_DIR)
    print(f"[quick paper] filled {len(filled)} blank fact(s) with placeholders:")
    for line in filled:
        print(f"  {line}")

    # Collect the claims-sidecar failures instead of letting the build raise on them.
    failures: list[str] = []
    sidecar_failures = vibepaper.build.sidecar_failures

    def failures_as_warnings(*args, **kwargs) -> list[str]:
        failures.extend(sidecar_failures(*args, **kwargs))
        return []

    vibepaper.build.sidecar_failures = failures_as_warnings

    config = vibepaper.build.load_config(PROJECT_ROOT / "paper.toml")
    config["name"] = f"{config['name']}_QUICK"
    config["facts_dir"] = str(QUICK_FACTS_DIR)
    config["build_dir"] = str(OUT_DIR / "build")
    vibepaper.build.run_build(config, PROJECT_ROOT, OUT_DIR, md=True, pdf_if_available=True)

    print(
        f"[quick paper] {len(failures)} claims assert(s) failed on the quick numbers "
        "(expected for a quick build; not checked here):"
    )
    for failure in failures:
        print(f"  {failure.splitlines()[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
