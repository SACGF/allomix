"""Read sequencing run-unit metadata from a VCF header (index-hopping check, issue #12).

The joint-calling pipeline stamps the admixture VCF with one ``##allomixRunUnit``
header line per admix sample, carrying that sample's sequencing run unit
(flowcell+lane) and whether it shares that unit with the host. Sharing a flowcell
lane is the index-hopping risk: barcode swapping at the sequencer can deposit
host reads into a co-loaded admix sample and fake a host-presence signal.

allomix never opens a BAM, so this header is the only path by which run-unit
identity reaches the tool. The metadata is optional: a VCF produced without it,
or one where the run unit was unrecoverable (SRA-renamed reads, no ``@RG PU``),
simply yields no entries and the index-hopping QC stays silent.

Writer side: ``pipeline/scripts/extract_run_units.py`` writes the matching
``##allomixRunUnit`` lines. Keep the header key and field names in sync.
"""

import re
from dataclasses import dataclass

from cyvcf2 import VCF

#: Header key for the per-sample run-unit lines. Custom (not the spec-reserved
#: ``##SAMPLE``) so it does not collide with that line's reserved fields.
HEADER_KEY = "allomixRunUnit"
_LINE = re.compile(r"^##" + HEADER_KEY + r"=<(.+)>\s*$")


@dataclass(frozen=True)
class RunUnitInfo:
    """Run-unit metadata for one sample, read from the VCF header.

    Attributes:
        run_unit: Flowcell:lane identifier, or None when not recorded.
        source: How it was derived upstream ("RG:PU", "readname"), or None.
        shares_run_with_host: True if this sample shares its run unit with the
            host (index-hopping risk), False if not, None when undetermined
            (e.g. the host run unit was unrecoverable).
    """

    run_unit: str | None
    source: str | None
    shares_run_with_host: bool | None


def _parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    v = value.strip().lower()
    if v == "true":
        return True
    if v == "false":
        return False
    return None


def _parse_fields(body: str) -> dict[str, str]:
    """Parse ``Key=Value,Key=Value`` from inside a structured header line.

    Our values (flowcell:lane, a source token, a boolean) never contain commas
    or spaces, so a plain comma split is sufficient; surrounding quotes are
    stripped defensively in case a writer adds them.
    """
    out: dict[str, str] = {}
    for kv in body.split(","):
        if "=" in kv:
            key, val = kv.split("=", 1)
            out[key.strip()] = val.strip().strip('"')
    return out


def read_run_units(vcf_path: str) -> dict[str, RunUnitInfo]:
    """Read all ``##allomixRunUnit`` lines from a VCF header.

    Args:
        vcf_path: Path to the (admixture) VCF.

    Returns:
        ``{sample_id: RunUnitInfo}`` for every line present. Empty when the VCF
        carries no run-unit metadata.
    """
    vcf = VCF(str(vcf_path))
    try:
        header = vcf.raw_header
    finally:
        vcf.close()

    out: dict[str, RunUnitInfo] = {}
    for line in header.splitlines():
        m = _LINE.match(line)
        if not m:
            continue
        fields = _parse_fields(m.group(1))
        sid = fields.get("ID")
        if not sid:
            continue
        out[sid] = RunUnitInfo(
            run_unit=fields.get("RunUnit"),
            source=fields.get("Source"),
            shares_run_with_host=_parse_bool(fields.get("SharesRunWithHost")),
        )
    return out


__all__ = ["HEADER_KEY", "RunUnitInfo", "read_run_units"]
