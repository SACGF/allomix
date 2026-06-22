#!/bin/bash
# Mix N BAMs at known target fractions into a single merged BAM suitable for
# joint calling as one admixture sample. Used to build semi-synthetic chimerism
# mixtures (host + 1 or 2 donors) from real reference BAMs (issue #5).
#
# Each input is subsampled independently, then merged, with all reads retagged to
# a single read-group sample name so the variant caller treats the result as one
# sample rather than several.
#
# Depth normalization: subsample fractions are sized from each input's ON-TARGET
# depth (summed per-base coverage over PANEL_BED) so the realized fraction of each
# component equals its target regardless of the inputs' depth ratios. Subsampling
# each input by a fixed fraction of its own reads (the naive approach) makes the
# realized fractions drift with the depth ratios. On-target reads, not raw BAM
# counts, because off-target rates differ between samples. For components with
# on-target depths D_i and target fractions t_i (which must sum to 1):
#     D0  = min_i( D_i / t_i )        # largest output depth that over-samples no input
#     s_i = t_i * D0 / D_i            # subsample fraction for component i
# Per-site (rather than total on-target) normalization would be marginally more
# faithful but is second-order and not practical with samtools' global
# --subsample; see issue #5 for the rationale.
#
# Usage:
#   scripts/mix_bams.sh OUTPUT_BAM PANEL_BED SAMPLE_NAME SEED BAM1 FRAC1 [BAM2 FRAC2 ...]
#
# Arguments:
#   OUTPUT_BAM    Path to write the mixed, indexed BAM to
#   PANEL_BED     Capture-panel BED; on-target depth is measured over this
#   SAMPLE_NAME   Read-group SM tag for the merged sample
#   SEED          Integer base seed; component i is subsampled with SEED+i so the
#                 sources are sampled independently and reproducibly
#   BAMi FRACi    One or more (BAM, target on-target fraction) pairs. FRACs are in
#                 (0, 1) and must sum to 1. At least 2 components are required.
#
# Examples:
#   # Two-person: host (minor) at 0.5%, donor background at 99.5%
#   scripts/mix_bams.sh out/mix.bam panel.bed SYN_F2_M1_f0005_rep1 100 \
#       donor.bam 0.995 host.bam 0.005
#   # Three-person: host 0.5%, two donors split the rest equally
#   scripts/mix_bams.sh out/mix3.bam panel.bed SYN3_F2_M1_M2_h0005_eq_rep1 100 \
#       donor1.bam 0.4975 donor2.bam 0.4975 host.bam 0.005

set -euo pipefail

if [ "$#" -lt 6 ]; then
    echo "Usage: $0 OUTPUT_BAM PANEL_BED SAMPLE_NAME SEED BAM1 FRAC1 [BAM2 FRAC2 ...]" >&2
    exit 1
fi

OUTPUT_BAM="$1"; shift
PANEL_BED="$1"; shift
SAMPLE_NAME="$1"; shift
SEED="$1"; shift

if [ $(( $# % 2 )) -ne 0 ]; then
    echo "Error: components must be given as 'BAM FRAC' pairs (odd argument count)" >&2
    exit 1
fi

BAMS=(); FRACS=()
while [ "$#" -gt 0 ]; do
    BAMS+=("$1"); FRACS+=("$2"); shift 2
done
N="${#BAMS[@]}"
if [ "$N" -lt 2 ]; then
    echo "Error: need at least 2 components" >&2
    exit 1
fi

# Each fraction strictly in (0, 1) (the depth normalization divides by t_i).
for t in "${FRACS[@]}"; do
    if ! awk -v x="$t" 'BEGIN { exit !(x+0 == x && x > 0 && x < 1) }'; then
        echo "Error: each fraction must be a number strictly between 0 and 1 (got: $t)" >&2
        exit 1
    fi
done

# Fractions must sum to 1 (within rounding tolerance).
SUM=$(printf '%s\n' "${FRACS[@]}" | awk '{ s += $1 } END { printf "%.8f", s }')
if ! awk -v s="$SUM" 'BEGIN { exit !(s > 0.999 && s < 1.001) }'; then
    echo "Error: component fractions must sum to 1 (got: $SUM)" >&2
    exit 1
fi

for f in "$PANEL_BED" "${BAMS[@]}"; do
    if [ ! -f "$f" ]; then
        echo "Error: input file not found: $f" >&2
        exit 1
    fi
done

# On-target depth (summed per-base coverage over the panel) of each input.
DEPTHS=()
for b in "${BAMS[@]}"; do
    d=$(samtools bedcov "$PANEL_BED" "$b" | awk '{ s += $NF } END { printf "%.0f", s }')
    if [ "$d" -le 0 ]; then
        echo "Error: zero on-target depth over $PANEL_BED for $b" >&2
        exit 1
    fi
    DEPTHS+=("$d")
done

# D0 = min_i(D_i / t_i): the largest output depth that over-samples no input.
D0=$(for i in $(seq 0 $((N - 1))); do printf '%s %s\n' "${DEPTHS[$i]}" "${FRACS[$i]}"; done \
    | awk 'BEGIN { d0 = -1 } { v = $1 / $2; if (d0 < 0 || v < d0) d0 = v } END { printf "%.6f", d0 }')

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$(dirname "$OUTPUT_BAM")"

echo "On-target depths: ${DEPTHS[*]}  target fractions: ${FRACS[*]}  output D0=$D0"

SUBBAMS=()
for i in $(seq 0 $((N - 1))); do
    subfrac=$(awk -v t="${FRACS[$i]}" -v d="${DEPTHS[$i]}" -v d0="$D0" \
        'BEGIN { x = t * d0 / d; if (x > 1) x = 1; printf "%.8f", x }')
    seed=$(( SEED + i ))
    sub="$TMPDIR/comp${i}.bam"
    echo "[comp $((i + 1))/$N] $(basename "${BAMS[$i]}") target=${FRACS[$i]} -> subsample $subfrac (seed $seed)"
    samtools view \
        --bam \
        --subsample "$subfrac" \
        --subsample-seed "$seed" \
        --output "$sub" \
        "${BAMS[$i]}"
    SUBBAMS+=("$sub")
done

MERGED="$TMPDIR/merged.bam"
echo "Merging $N subsampled BAMs"
samtools merge -f -o "$MERGED" "${SUBBAMS[@]}"

echo "Unifying read groups under sample '$SAMPLE_NAME'"
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

echo "Done: $OUTPUT_BAM (sample: $SAMPLE_NAME, $N components)"
