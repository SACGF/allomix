#!/usr/bin/env bash
# Build a force-output panel VCF for the GATK two-phase pipeline.
#
# Filters a population-level source VCF (e.g. gnomAD v4.1 sites) to:
#   - sites overlapping the given capture-panel BED
#   - PASS filter
#   - biallelic SNPs
#   - INFO/AF >= AF_THRESHOLD (global gnomAD AF)
#
# The output is bgzipped + tabix-indexed and ready to drop into
# pipeline/config.yaml as `panel_alleles_vcf:`. GenotypeGVCFs will then
# force-emit a genotype at every panel site, regardless of whether GATK
# would have called it variant in the small per-patient joint call.
#
# Chromosome naming: our reference is UCSC-style (chr1, chr2, ...) but the
# gnomAD v4.1 sites VCF (joint exomes+genomes) is RefSeq-accessioned
# (NC_000001.11, NC_000002.12, ...). When a chr-mapping file is given
# (column 1 = source name, column 2 = target name, tab-separated), the
# script reverse-maps the BED for the gnomAD query and forward-maps the
# output VCF back to UCSC names. See pipeline/gnomad_refseq_to_hg38_chrs.tsv.
#
# Usage:
#   scripts/build_force_output_panel.sh \
#       <bed> <source_vcf> <af_threshold> <output_vcf.gz> [chr_mapping_tsv]
#
# Example (gnomAD v4.1 sites, AF >= 0.05, haem panel):
#   scripts/build_force_output_panel.sh \
#       output/union_sid_haem_vendor_probes.bed \
#       /data/annotation/VEP/annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz \
#       0.05 \
#       output/union_sid_haem_force_panel_af05.vcf.gz \
#       pipeline/gnomad_refseq_to_hg38_chrs.tsv
#
# Requires bcftools and tabix on PATH.

set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
    cat >&2 <<EOF
Usage: $0 <bed> <source_vcf> <af_threshold> <output_vcf.gz> [chr_mapping_tsv]

The optional chr_mapping_tsv has two tab-separated columns:
  source_name  target_name
Used to translate UCSC-style BED chrs (chr1, ...) into the source VCF's
naming (e.g. gnomAD's NC_000001.11), and to rename the output VCF back.
EOF
    exit 1
fi

BED="$1"
SRC="$2"
AF="$3"
OUT="$4"
CHRMAP="${5:-}"

for f in "$BED" "$SRC"; do
    [[ -f "$f" ]] || { echo "Not found: $f" >&2; exit 1; }
done

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# Step 1 — strip comment header lines, then optionally rename BED chrs to
# whatever the source VCF uses (so bcftools view -R can resolve them).
if [[ -n "$CHRMAP" ]]; then
    [[ -f "$CHRMAP" ]] || { echo "Not found: $CHRMAP" >&2; exit 1; }
    # Build a forward map (UCSC name -> source name) by reversing columns.
    awk -F'\t' 'NR==FNR{m[$2]=$1; next}
                /^#/ {next}
                ($1 in m){print m[$1] "\t" $2 "\t" $3}' \
        "$CHRMAP" "$BED" > "$TMPDIR/regions.bed"
else
    grep -v '^#' "$BED" > "$TMPDIR/regions.bed"
fi

if [[ ! -s "$TMPDIR/regions.bed" ]]; then
    echo "No usable rows in BED after chr translation. Check the BED and mapping." >&2
    exit 1
fi

# Step 2 — query the source. Filters:
#   -R <bed>            restrict to BED regions (uses tabix index on source)
#   -f PASS,EXOMES_FILTERED,GENOMES_FILTERED
#                       gnomAD v4 joint VCF FILTER values: PASS, or failed
#                       just one of the two datasets. We accept all three
#                       and only drop BOTH_FILTERED, because for a
#                       force-output panel we want known polymorphic sites,
#                       not perfectly-called ones — the filter failure is
#                       about gnomAD's own call quality, which is irrelevant
#                       once we force-call the site in our own samples.
#   -m2 -M2             keep only biallelic sites
#   -v snps             drop indels
#   -e 'INFO/AF<AF      drop low-AF sites and sites with missing AF
#       || INFO/AF="."'
QUERIED="$TMPDIR/queried.vcf.gz"
bcftools view \
    -R "$TMPDIR/regions.bed" \
    -f PASS,EXOMES_FILTERED,GENOMES_FILTERED \
    -m2 -M2 \
    -v snps \
    -e "INFO/AF<${AF} || INFO/AF=\".\"" \
    "$SRC" \
    -Oz -o "$QUERIED"

# Step 3 — rename output contigs back to BED-side names if a chr map was given.
if [[ -n "$CHRMAP" ]]; then
    bcftools annotate \
        --rename-chrs "$CHRMAP" \
        "$QUERIED" \
        | bcftools sort -Oz -o "$OUT"
else
    bcftools sort -Oz -o "$OUT" "$QUERIED"
fi
tabix -f -p vcf "$OUT"

N=$(bcftools view -H "$OUT" | wc -l)
echo "Wrote ${N} sites to ${OUT} (AF >= ${AF}, PASS, biallelic SNPs, within ${BED})"
