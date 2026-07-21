#!/bin/bash
# Regenerate the three error-table arms for the supplementary comparison (#49).
#
#   flat         built-in flat --error-rate default, no table
#   per_mixture  each mixture's own estimate-errors table (2 individuals)
#   pooled       one table pooled over all 7 reference individuals
#
# Each arm is a full run_srp434573_allomix.py pass. They are run CONCURRENTLY,
# each into its own ALLOMIX_OUT_DIR, because the script otherwise writes fixed
# paths under output/ and the arms would overwrite each other. Expect ~30 min
# wall clock for all three (they contend during the detect subprocesses, so this
# is faster than sequential but not 3x).
#
# The paper does NOT run this. Only the flat and per_mixture arms are committed
# under paper/public_data/SRP434573/error_table_arms/; the pooled arm is the
# ordinary build output, because the pipeline defaults to the pooled table and
# srp434573_allomix already writes exactly that. The Snakemake rule
# error_table_arms_facts only summarises them, so a fresh checkout builds the
# figure in seconds. Re-run this only when the estimator or the committed
# genotype snapshot changes, then refresh the snapshot (last step below).
#
# Usage, from the repo root:
#   bash paper/scripts/run_error_table_arms.sh
set -uo pipefail

ARMS=output/error_table_arms
SNAPSHOT=paper/public_data/SRP434573/error_table_arms
mkdir -p "$ARMS"

run_arm () {  # run_arm <name> [ENV=VAL ...]
  local name=$1; shift
  echo "=== START $name $(date +%H:%M:%S) ==="
  env ALLOMIX_OUT_DIR="$ARMS/$name" "$@" \
      python3 paper/scripts/run_srp434573_allomix.py > "$ARMS/$name.log" 2>&1
  echo "=== END   $name $(date +%H:%M:%S) rc=$? ==="
}

run_arm pooled &
run_arm per_mixture ALLOMIX_FORCE_PER_MIXTURE_ERROR_TABLE=1 &
run_arm flat        ALLOMIX_NO_ERROR_TABLE=1 &
wait

echo "=== ALL ARMS DONE $(date +%H:%M:%S) ==="
for a in flat per_mixture pooled; do
  echo "  $a: $(ls "$ARMS/$a"/*.tsv 2>/dev/null | wc -l) tsv"
done

cat <<EOF

To refresh the committed snapshot the paper builds from:

    for a in flat per_mixture; do          # pooled is the build's own output
        mkdir -p $SNAPSHOT/\$a
        cp $ARMS/\$a/*.tsv $SNAPSHOT/\$a/
    done
    snakemake -s paper/Snakefile --cores 4 output/facts/error_table_arms.csv

Note generate_error_table_arms_facts.py prefers $ARMS when it exists, so the
figure will already reflect this run before you copy anything.
EOF
