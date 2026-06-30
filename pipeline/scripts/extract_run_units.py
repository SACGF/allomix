"""Extract the sequencing run unit (flowcell + lane) per sample from BAMs.

This is the BAM-side half of the index-hopping check (issue #12). allomix itself
is VCF-first and never opens a BAM, so the run/flowcell identity an index-hopping
flag needs is pulled here, in the joint-calling pipeline, and written as a VCF
header fragment of ``##allomixRunUnit`` lines. The merge step stamps those onto
the admixture VCF, so the metadata travels inside the VCF (no losable sidecar)
and allomix reads it straight back from the header (``allomix.qc.runmeta``).

Index hopping (barcode swapping on patterned-flowcell Illumina machines) leaks
reads between samples clustered in the same lane of the same flowcell. So the
unit that matters is flowcell + lane: two samples sharing it can cross-contaminate
at the sequencer. If the host shares a run unit with a post-transplant admixture
sample, hopped host reads can fake a host-presence signal, which is the
clinically dangerous case. This script records each admix sample's run unit and
flags the ones that share it with the host.

The output is optional by construction: a header line is emitted only for an
admix sample whose run unit is recoverable, so when nothing useful can be pulled
(e.g. SRA-renamed reads with no PU) the fragment is empty, the merge adds
nothing, and the downstream flag stays silent rather than reporting a false
"no risk".

Run unit is read, in order of preference:

  1. ``@RG`` ``PU`` (platform unit) tags in the BAM header. PU is the canonical
     SAM field for flowcell-barcode.lane and is what Picard / most aligners
     populate. The flowcell and lane (first two dot/colon-separated fields) are
     taken; the per-sample barcode that may follow is dropped.
  2. The flowcell and lane parsed from the first read name, for BAMs whose header
     carries no PU. Illumina names are
     ``instrument:run:flowcell:lane:tile:x:y`` (7 colon fields); fields 3 and 4
     give flowcell:lane.
  3. ``unknown`` when neither is available. This is the expected result for
     SRA-derived BAMs (SRP434573), where the deposited reads are renamed
     ``<accession>.<n>`` and the header carries no PU, so the flowcell is
     unrecoverable. The flag then degrades to "cannot determine" rather than a
     false negative.

A BAM merged across lanes / flowcells can carry several run units; all distinct
ones are kept and a share is any non-empty intersection.

Output is a VCF header fragment, one line per admix sample with a recoverable
run unit:

    ##allomixRunUnit=<ID=<sample>,RunUnit=<flowcell:lane>,Source=<RG:PU|readname>,SharesRunWithHost=<true|false|unknown>>

``RunUnit`` is ";"-joined when a sample spans more than one unit.
``SharesRunWithHost`` is ``unknown`` when the host run unit is itself
unrecoverable (the comparison cannot be made). Keep the header key and field
names in sync with ``src/allomix/runmeta.py``, which parses them.
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

# An Illumina read name with at least this many colon fields carries the
# flowcell in field 3 and the lane in field 4 (instrument:run:flowcell:lane:...).
_ILLUMINA_MIN_FIELDS = 7
_UNKNOWN = "unknown"
# Header key parsed by allomix.qc.runmeta.read_run_units — keep in sync.
_HEADER_KEY = "allomixRunUnit"


def _run_units_from_header(bam: Path, samtools: str) -> set[str]:
    """Flowcell:lane run units parsed from every ``@RG`` ``PU`` tag in the header."""
    try:
        header = subprocess.run(
            [samtools, "view", "-H", str(bam)],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return set()

    units: set[str] = set()
    for line in header.splitlines():
        if not line.startswith("@RG"):
            continue
        for field in line.split("\t"):
            if field.startswith("PU:"):
                pu = field[3:]
                units.add(_normalise_unit(pu))
    units.discard("")
    return units


def _run_unit_from_first_read(bam: Path, samtools: str) -> str | None:
    """Flowcell:lane parsed from the first read name, or None if not Illumina-shaped."""
    try:
        # `samtools view | head -1` closes the pipe after one record; samtools
        # takes SIGPIPE, which is expected, so the return code is not checked.
        proc = subprocess.run(
            f"{samtools} view {bam} 2>/dev/null | head -n 1",
            shell=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return None
    line = proc.stdout.strip()
    if not line:
        return None
    qname = line.split("\t", 1)[0]
    fields = qname.split(":")
    if len(fields) >= _ILLUMINA_MIN_FIELDS:
        flowcell, lane = fields[2], fields[3]
        if flowcell and lane:
            return f"{flowcell}:{lane}"
    return None


def _normalise_unit(pu: str) -> str:
    """Reduce a PU tag to ``flowcell:lane``.

    PU is conventionally ``flowcell.lane.barcode`` (Picard) but ``.``/``:``
    separators and a missing barcode are both seen. The first two
    dot/colon-separated tokens are the flowcell and (if present) lane; the
    trailing per-sample barcode is dropped so co-loaded samples match.
    """
    tokens = pu.replace(".", ":").split(":")
    tokens = [t for t in tokens if t]
    if not tokens:
        return ""
    return ":".join(tokens[:2])


def run_units_for_bam(bam: Path, samtools: str) -> tuple[set[str], str]:
    """Return ``(run_units, source)`` for one BAM.

    Prefers the header PU tags; falls back to the first read name; otherwise
    ``({"unknown"}, "unknown")``.
    """
    header_units = _run_units_from_header(bam, samtools)
    if header_units:
        return header_units, "RG:PU"
    read_unit = _run_unit_from_first_read(bam, samtools)
    if read_unit:
        return {read_unit}, "readname"
    return {_UNKNOWN}, _UNKNOWN


def _read_csv(csv_path: Path) -> list[dict[str, str]]:
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        missing = {"sample_id", "bam_filename", "sample_type"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{csv_path} is missing column(s): {', '.join(sorted(missing))}")
        return list(reader)


def build_metadata(csv_path: Path, samtools: str) -> list[dict[str, str]]:
    """Extract run units for every sample in a patient CSV and flag host shares."""
    rows = _read_csv(csv_path)

    units_by_sample: dict[str, set[str]] = {}
    source_by_sample: dict[str, str] = {}
    type_by_sample: dict[str, str] = {}
    host_units: set[str] = set()

    for row in rows:
        sid = row["sample_id"]
        stype = row["sample_type"].strip().upper()
        units, source = run_units_for_bam(Path(row["bam_filename"]), samtools)
        units_by_sample[sid] = units
        source_by_sample[sid] = source
        type_by_sample[sid] = stype
        if stype == "HOST":
            host_units |= {u for u in units if u != _UNKNOWN}

    out: list[dict[str, str]] = []
    for sid, units in units_by_sample.items():
        stype = type_by_sample[sid]
        known = {u for u in units if u != _UNKNOWN}
        if stype == "HOST":
            shares = "NA"
        elif not host_units or not known:
            # Either the host run unit or this sample's is unrecoverable, so the
            # comparison cannot be made; report undetermined, not "no".
            shares = "unknown"
        else:
            shares = "true" if (known & host_units) else "false"
        out.append(
            {
                "sample_id": sid,
                "sample_type": stype,
                "run_units": ";".join(sorted(units)),
                "run_unit_source": source_by_sample[sid],
                "shares_run_with_host": shares,
            }
        )
    return out


def to_header_lines(rows: list[dict[str, str]]) -> list[str]:
    """Build ``##allomixRunUnit`` header lines for admix samples with a known unit.

    Only ADMIX samples are emitted (they are the samples in the admix VCF the
    fragment annotates), and only when their run unit is recoverable, so an
    all-unknown run produces an empty fragment.
    """
    lines: list[str] = []
    for row in rows:
        if row["sample_type"] != "ADMIX":
            continue
        if row["run_units"] == _UNKNOWN:
            continue
        lines.append(
            f"##{_HEADER_KEY}=<ID={row['sample_id']},"
            f"RunUnit={row['run_units']},"
            f"Source={row['run_unit_source']},"
            f"SharesRunWithHost={row['shares_run_with_host']}>"
        )
    return lines


def write_header(lines: list[str], out_path: Path) -> None:
    """Write the header fragment (possibly empty when nothing was recoverable)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for line in lines:
            f.write(line + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path, help="Patient sample CSV.")
    parser.add_argument(
        "--out", required=True, type=Path, help="Output VCF header fragment (.hdr)."
    )
    parser.add_argument("--samtools", default="samtools", help="samtools binary.")
    args = parser.parse_args(argv)

    rows = build_metadata(args.csv, args.samtools)
    lines = to_header_lines(rows)
    write_header(lines, args.out)

    n_admix = sum(1 for r in rows if r["sample_type"] == "ADMIX")
    n_share = sum(1 for r in rows if r["shares_run_with_host"] == "true")
    sys.stderr.write(
        f"{args.csv.name}: {len(lines)}/{n_admix} admix samples with a recoverable "
        f"run unit, {n_share} sharing one with the host.\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
