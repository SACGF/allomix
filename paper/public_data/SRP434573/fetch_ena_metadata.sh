#!/usr/bin/env bash
# Fetch the ENA run metadata and FASTQ download list for an SRA/ENA study, so
# the SRP434573 input mapping is reproducible (not a hand-fetched one-off).
#
# Produces, next to this script:
#   ena_runs.tsv        run -> sample_alias + library/platform/read-count fields
#                       (consumed by make_sample_csvs.py to build the CSVs)
#   download_fastqs.sh  wget commands for every run's FASTQ (the ENA
#                       "download all" list, regenerated from fastq_ftp)
#
# Usage:
#   ./fetch_ena_metadata.sh            # defaults to SRP434573
#   ./fetch_ena_metadata.sh SRP434573  # or any other study accession
set -euo pipefail

ACCESSION="${1:-SRP434573}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
API="https://www.ebi.ac.uk/ena/portal/api/filereport"

# Fields used by make_sample_csvs.py (sample_alias encodes the mixture design)
# plus context fields kept for provenance/QC.
FIELDS="run_accession,sample_accession,sample_alias,sample_title,library_name,library_strategy,library_source,instrument_platform,read_count,scientific_name"

echo "Fetching run metadata for ${ACCESSION} ..." >&2
curl -fsSL \
    "${API}?accession=${ACCESSION}&result=read_run&fields=${FIELDS}&format=tsv" \
    -o "${HERE}/ena_runs.tsv"
echo "  wrote ena_runs.tsv ($(( $(wc -l < "${HERE}/ena_runs.tsv") - 1 )) runs)" >&2

echo "Fetching FASTQ download list for ${ACCESSION} ..." >&2
# ENA always returns run_accession as the first column, so request it explicitly
# and take column 2 (fastq_ftp). fastq_ftp holds host-relative FTP paths
# (semicolon-separated for paired runs); prefix ftp:// and emit one `wget -nc`
# per file, matching the ENA "download all" script format.
curl -fsSL \
    "${API}?accession=${ACCESSION}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
    | tail -n +2 \
    | cut -f2 \
    | tr ';' '\n' \
    | sed '/^[[:space:]]*$/d; s#^#wget -nc ftp://#' \
    > "${HERE}/download_fastqs.sh"
echo "  wrote download_fastqs.sh ($(wc -l < "${HERE}/download_fastqs.sh") files)" >&2
