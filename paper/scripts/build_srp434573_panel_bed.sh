#!/usr/bin/env bash
# Recover the SRP434573 MIP capture panel BED from BAM coverage.
#
# Thin wrapper around the generic scripts/recover_panel_bed.py, passing the
# SRP434573-specific paths and thresholds. The SRP434573 thesis (Chu Xufeng,
# HUST 2024) describes a ~1062-autosomal-SNP MIP panel but publishes no
# coordinates; because this is a MIP/amplicon assay the panel self-recovers from
# coverage (see SACGF/allomix issue #16 and paper/public_data/SRP434573/README.md).
#
# A position is kept when covered at >=100x (MAPQ/BASEQ >=20) in at least 50 of
# the 64 runs, and adjacent kept positions are merged into one interval per
# amplicon. This recovers 1052 hg38 intervals (1025 autosomal + 27 on chrX),
# matching the issue #16 laptop probe (~1053) and the stated ~1062 SNPs.
#
# Usage (from repo root):
#   paper/scripts/build_srp434573_panel_bed.sh
#
# Override the BAM glob or output path by passing extra args, which are forwarded
# to recover_panel_bed.py (e.g. --min-samples 40 to tolerate more dropout).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

python "${REPO_ROOT}/scripts/recover_panel_bed.py" \
    --bam-glob 'output/bam/*.bam' \
    --out paper/public_data/SRP434573/SRP434573.bed \
    --min-samples 50 \
    "$@"
