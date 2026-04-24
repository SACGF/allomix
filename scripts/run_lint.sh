#!/usr/bin/env bash
# Run pylint over the project using .pylintrc at the repo root.
# Usage: scripts/run_lint.sh [extra pylint args]
set -euo pipefail

cd "$(dirname "$0")/.."

exec pylint \
    --rcfile=.pylintrc \
    src/allomix \
    tests \
    scripts \
    paper/scripts \
    "$@"
