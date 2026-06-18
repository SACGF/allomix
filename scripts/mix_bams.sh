#!/bin/bash
# Mix two BAMs at a target donor fraction, producing a single merged BAM
# suitable for joint calling as one admixture sample.
#
# Each input BAM is subsampled independently, then merged, with all reads
# retagged to a single read-group sample name so the variant caller treats
# the result as one sample rather than two.
#
# Usage:
#   scripts/mix_bams.sh HOST_BAM DONOR_BAM DONOR_FRACTION OUTPUT_BAM [SAMPLE_NAME] [SEED]
#
# Arguments:
#   HOST_BAM         Path to the host/majority BAM
#   DONOR_BAM        Path to the donor/minority BAM
#   DONOR_FRACTION   Target fraction of reads from DONOR_BAM (0.0 - 1.0)
#   OUTPUT_BAM       Path to write the mixed, indexed BAM to
#   SAMPLE_NAME      Optional read-group SM tag (default: MIX_F<fraction>)
#   SEED             Optional integer seed for reproducible subsampling
#                    (default: 42). Host uses SEED; donor uses SEED+1 so the
#                    two sources are sampled independently.
#
# Example:
#   scripts/mix_bams.sh host.bam donor.bam 0.05 output/mix_f005_rep1.bam MIX_F005_REP1 42

set -euo pipefail

if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 HOST_BAM DONOR_BAM DONOR_FRACTION OUTPUT_BAM [SAMPLE_NAME] [SEED]" >&2
    exit 1
fi

HOST_BAM="$1"
DONOR_BAM="$2"
DONOR_FRACTION="$3"
OUTPUT_BAM="$4"
SEED="${6:-42}"

if [ -n "${5:-}" ]; then
    SAMPLE_NAME="$5"
else
    # Default sample name: "MIX_F0_05" for fraction 0.05
    SAMPLE_NAME="MIX_F$(echo "$DONOR_FRACTION" | tr '.' '_')"
fi

# Validate fraction is a number in [0, 1]
if ! awk -v x="$DONOR_FRACTION" 'BEGIN { exit !(x+0 == x && x >= 0 && x <= 1) }'; then
    echo "Error: DONOR_FRACTION must be a number between 0 and 1 (got: $DONOR_FRACTION)" >&2
    exit 1
fi

for f in "$HOST_BAM" "$DONOR_BAM"; do
    if [ ! -f "$f" ]; then
        echo "Error: input BAM not found: $f" >&2
        exit 1
    fi
done

HOST_FRACTION=$(awk -v d="$DONOR_FRACTION" 'BEGIN { printf "%.6f", 1 - d }')
HOST_SEED="$SEED"
DONOR_SEED=$((SEED + 1))

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

HOST_SUB="$TMPDIR/host.sub.bam"
DONOR_SUB="$TMPDIR/donor.sub.bam"
MERGED="$TMPDIR/merged.bam"

mkdir -p "$(dirname "$OUTPUT_BAM")"

echo "[1/4] Subsampling host to fraction $HOST_FRACTION (seed $HOST_SEED)"
samtools view \
    --bam \
    --subsample "$HOST_FRACTION" \
    --subsample-seed "$HOST_SEED" \
    --output "$HOST_SUB" \
    "$HOST_BAM"

echo "[2/4] Subsampling donor to fraction $DONOR_FRACTION (seed $DONOR_SEED)"
samtools view \
    --bam \
    --subsample "$DONOR_FRACTION" \
    --subsample-seed "$DONOR_SEED" \
    --output "$DONOR_SUB" \
    "$DONOR_BAM"

echo "[3/4] Merging subsampled BAMs"
samtools merge \
    -f \
    -o "$MERGED" \
    "$HOST_SUB" "$DONOR_SUB"

echo "[4/4] Unifying read groups under sample '$SAMPLE_NAME'"
samtools addreplacerg \
    -m overwrite_all \
    -r "ID:${SAMPLE_NAME}" \
    -r "SM:${SAMPLE_NAME}" \
    -r "LB:${SAMPLE_NAME}" \
    -r "PL:ILLUMINA" \
    --output-fmt BAM \
    -o "$OUTPUT_BAM" \
    "$MERGED"

samtools index "$OUTPUT_BAM"

echo "Done: $OUTPUT_BAM (sample: $SAMPLE_NAME)"
