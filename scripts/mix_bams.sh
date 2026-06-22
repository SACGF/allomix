#!/bin/bash
# Mix two BAMs at a target donor (minor) fraction, producing a single merged BAM
# suitable for joint calling as one admixture sample.
#
# Each input BAM is subsampled independently, then merged, with all reads
# retagged to a single read-group sample name so the variant caller treats the
# result as one sample rather than two.
#
# Depth normalization: the subsample fractions are sized from each input's
# ON-TARGET depth (summed per-base coverage over PANEL_BED) so the realized minor
# fraction equals DONOR_FRACTION regardless of the two libraries' depth ratio.
# Subsampling each input by a fixed fraction of its own reads (the naive
# approach) makes the realized fraction drift with the depth ratio. On-target
# reads, not raw BAM counts, because off-target rates differ between samples.
# Given on-target depths Dmajor, Dminor and target minor fraction t:
#     D0      = min( Dminor / t , Dmajor / (1 - t) )   # largest non-oversampling output depth
#     f_major = (1 - t) * D0 / Dmajor
#     f_minor =       t * D0 / Dminor
# Per-site (rather than total on-target) normalization would be marginally more
# faithful but is second-order and not practical with samtools' global
# --subsample; see issue #5 for the rationale.
#
# Usage:
#   scripts/mix_bams.sh HOST_BAM DONOR_BAM DONOR_FRACTION OUTPUT_BAM PANEL_BED [SAMPLE_NAME] [SEED]
#
# Arguments:
#   HOST_BAM         Path to the host/majority (background) BAM
#   DONOR_BAM        Path to the donor/minority (titrated) BAM
#   DONOR_FRACTION   Target fraction of on-target reads from DONOR_BAM (0 < f < 1)
#   OUTPUT_BAM       Path to write the mixed, indexed BAM to
#   PANEL_BED        Capture-panel BED; on-target depth is measured over this
#   SAMPLE_NAME      Optional read-group SM tag (default: MIX_F<fraction>)
#   SEED             Optional integer seed for reproducible subsampling
#                    (default: 42). Host uses SEED; donor uses SEED+1 so the
#                    two sources are sampled independently.
#
# Example:
#   scripts/mix_bams.sh host.bam donor.bam 0.005 output/mix_f0005_rep1.bam \
#       paper/public_data/SRP434573/SRP434573.bed MIX_F0005_REP1 100

set -euo pipefail

if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
    echo "Usage: $0 HOST_BAM DONOR_BAM DONOR_FRACTION OUTPUT_BAM PANEL_BED [SAMPLE_NAME] [SEED]" >&2
    exit 1
fi

HOST_BAM="$1"
DONOR_BAM="$2"
DONOR_FRACTION="$3"
OUTPUT_BAM="$4"
PANEL_BED="$5"
SEED="${7:-42}"

if [ -n "${6:-}" ]; then
    SAMPLE_NAME="$6"
else
    # Default sample name: "MIX_F0_005" for fraction 0.005
    SAMPLE_NAME="MIX_F$(echo "$DONOR_FRACTION" | tr '.' '_')"
fi

# Validate fraction is a number strictly in (0, 1): the depth-normalized
# subsampling divides by both t and (1 - t).
if ! awk -v x="$DONOR_FRACTION" 'BEGIN { exit !(x+0 == x && x > 0 && x < 1) }'; then
    echo "Error: DONOR_FRACTION must be a number strictly between 0 and 1 (got: $DONOR_FRACTION)" >&2
    exit 1
fi

for f in "$HOST_BAM" "$DONOR_BAM" "$PANEL_BED"; do
    if [ ! -f "$f" ]; then
        echo "Error: input file not found: $f" >&2
        exit 1
    fi
done

# On-target depth (summed per-base coverage over the panel) of each pure BAM.
D_MAJOR=$(samtools bedcov "$PANEL_BED" "$HOST_BAM" | awk '{s += $NF} END { printf "%.0f", s }')
D_MINOR=$(samtools bedcov "$PANEL_BED" "$DONOR_BAM" | awk '{s += $NF} END { printf "%.0f", s }')

if [ "$D_MAJOR" -le 0 ] || [ "$D_MINOR" -le 0 ]; then
    echo "Error: zero on-target depth over $PANEL_BED (major=$D_MAJOR minor=$D_MINOR)" >&2
    exit 1
fi

# Subsample fractions that land the realized minor fraction on DONOR_FRACTION at
# output on-target depth D0 (the largest depth that over-samples neither input).
read -r MAJOR_SUBFRAC MINOR_SUBFRAC D0 < <(awk \
    -v t="$DONOR_FRACTION" -v dmin="$D_MINOR" -v dmaj="$D_MAJOR" 'BEGIN {
        d0a = dmin / t
        d0b = dmaj / (1 - t)
        d0  = (d0a < d0b) ? d0a : d0b
        fmaj = (1 - t) * d0 / dmaj
        fmin =       t * d0 / dmin
        if (fmaj > 1) fmaj = 1
        if (fmin > 1) fmin = 1
        printf "%.8f %.8f %.0f\n", fmaj, fmin, d0
    }')

HOST_SEED="$SEED"
DONOR_SEED=$((SEED + 1))

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

HOST_SUB="$TMPDIR/host.sub.bam"
DONOR_SUB="$TMPDIR/donor.sub.bam"
MERGED="$TMPDIR/merged.bam"

mkdir -p "$(dirname "$OUTPUT_BAM")"

echo "On-target depth: major=$D_MAJOR minor=$D_MINOR -> output D0=$D0 (target minor fraction $DONOR_FRACTION)"

echo "[1/4] Subsampling host (major) to fraction $MAJOR_SUBFRAC (seed $HOST_SEED)"
samtools view \
    --bam \
    --subsample "$MAJOR_SUBFRAC" \
    --subsample-seed "$HOST_SEED" \
    --output "$HOST_SUB" \
    "$HOST_BAM"

echo "[2/4] Subsampling donor (minor) to fraction $MINOR_SUBFRAC (seed $DONOR_SEED)"
samtools view \
    --bam \
    --subsample "$MINOR_SUBFRAC" \
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
