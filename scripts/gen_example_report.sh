#!/usr/bin/env bash
#
# Regenerate the example HTML reports committed under docs/examples/.
#
# Runs the allomix CLI on the de-identified public SRP434573 mixtures so the
# example reports linked from the README can be rebuilt deterministically. This
# is the bash counterpart of gen_example_report.py: it shows the exact command
# lines a user would type, nothing hidden behind a Python wrapper.
#
# The generation timestamp is pinned (--report-timestamp) so the output is
# byte-stable run to run; the only other nondeterminism, the timeline trend-chart
# PNG, is fixed for fixed input data.
#
# Two examples are produced:
#
#   1. Single-sample report (the common case): the lowest real fraction in the
#      ladder, alias 1_199_F2-M1_v1. The titrated minor contributor F2 is the
#      HOST and the background M1 is the DONOR (see
#      paper/scripts/run_srp434573_allomix.py), so the monitored 0.5% quantity is
#      the HOST fraction. The MLE reads slightly low because of the
#      donor-homozygous contamination background documented in
#      claude/srp434573_figure_review_findings.md; the host-presence test still
#      gates it correctly.
#
#   2. Timeline report (secondary): the whole F2-into-M1 dilution ladder fed to
#      the timeline mode to show the trend chart. These are a titration series,
#      not serial timepoints from one patient, so the example is labelled a
#      dilution series.
#
# Usage:
#   scripts/gen_example_report.sh           # writes docs/examples/*.html
#   scripts/gen_example_report.sh --check   # regenerate, then fail if anything changed
#
# The timeline report needs matplotlib for the trend chart:
#   pip install -e ".[report]"
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$REPO_ROOT/paper/public_data/SRP434573/genotypes"
OUT_DIR="$REPO_ROOT/docs/examples"

# Pinned so the committed reports are byte-reproducible.
TIMESTAMP="2026-06-29 00:00:00 (example build)"

PANEL="$DATA/mix_F2_into_M1.SRP434573.vcf.gz"
ADMIX="$DATA/mix_F2_into_M1.admix.vcf.gz"
ERRORS="$DATA/mix_F2_into_M1.error_table.tsv"

SINGLE_1PCT_OUT="$OUT_DIR/srp434573_single_sample_1pct.html"
SINGLE_1PCT_CSV="$OUT_DIR/srp434573_single_sample_1pct.markers.csv"
SINGLE_OUT="$OUT_DIR/srp434573_single_sample.html"
SINGLE_CSV="$OUT_DIR/srp434573_single_sample.markers.csv"
TIMELINE_OUT="$OUT_DIR/srp434573_dilution_series.html"

mkdir -p "$OUT_DIR"

# 1a. Single-sample report at 1% (the headline example): the 1_99 rung, where
#     the titrated host sits at 1% and the estimate recovers it cleanly.
allomix monitor \
  --genotype-vcf "$PANEL" \
  --admix-vcf "$ADMIX" \
  --error-table "$ERRORS" \
  --host-sample F2 \
  --donor-sample M1 \
  --sample 1_99_F2-M1_v1 \
  --html "$SINGLE_1PCT_OUT" \
  --marker-csv "$SINGLE_1PCT_CSV" \
  --recipient-id "SRP434573 demo (F2 into M1, 1%)" \
  --donor-relationship unrelated \
  --report-timestamp "$TIMESTAMP"

# 1b. Single-sample report at 0.5% (near the panel's contamination floor): the
#     1_199 rung, kept to show the low-end behaviour and host-presence gating.
allomix monitor \
  --genotype-vcf "$PANEL" \
  --admix-vcf "$ADMIX" \
  --error-table "$ERRORS" \
  --host-sample F2 \
  --donor-sample M1 \
  --sample 1_199_F2-M1_v1 \
  --html "$SINGLE_OUT" \
  --marker-csv "$SINGLE_CSV" \
  --recipient-id "SRP434573 demo (F2 into M1)" \
  --donor-relationship unrelated \
  --report-timestamp "$TIMESTAMP"

# 2. Timeline report across the F2-into-M1 dilution ladder (high host fraction
#    1_9 = 10% down to low 1_199 = 0.5%). One --sample per rung, in order.
allomix timeline \
  --genotype-vcf "$PANEL" \
  --admix-vcf "$ADMIX" \
  --error-table "$ERRORS" \
  --host-sample F2 \
  --donor-sample M1 \
  --sample 1_9_F2-M1_v1 \
  --sample 1_19_F2-M1_v1 \
  --sample 1_39_F2-M1_v1 \
  --sample 1_79_F2-M1_v1 \
  --sample 1_99_F2-M1_v1 \
  --sample 1_199_F2-M1_v1 \
  --html "$TIMELINE_OUT" \
  --recipient-id "SRP434573 dilution series (F2 into M1)" \
  --report-timestamp "$TIMESTAMP"

for f in "$SINGLE_1PCT_OUT" "$SINGLE_1PCT_CSV" "$SINGLE_OUT" "$SINGLE_CSV" "$TIMELINE_OUT"; do
  echo "wrote ${f#"$REPO_ROOT"/}"
done

# --check: the reports are committed, so staleness is just an uncommitted diff.
if [[ "${1:-}" == "--check" ]]; then
  if ! git -C "$REPO_ROOT" diff --quiet -- "$OUT_DIR"; then
    echo "Example reports are stale: commit the regenerated docs/examples/" >&2
    git -C "$REPO_ROOT" --no-pager diff --stat -- "$OUT_DIR" >&2
    exit 1
  fi
fi
