#!/usr/bin/env bash
# Frozen from scripts/funannotate_clean_V3.sh (runner prefix configurable)
set -euo pipefail
: "${GENOME:?}" "${OUT_FASTA:?}"
FUNANNOTATE_RUNNER="${FUNANNOTATE_RUNNER:-}"
mkdir -p "$(dirname "${OUT_FASTA}")"
unset PYTHONPATH || true
export PYTHONNOUSERSITE=1
export LC_ALL="${LC_ALL:-C.UTF-8}"
export LANG="${LANG:-C.UTF-8}"
# shellcheck disable=SC2086
eval "${FUNANNOTATE_RUNNER}" funannotate clean \
  -i "${GENOME}" \
  -o "${OUT_FASTA}" \
  --minlen 15000000 \
  --pident 95 \
  --cov 95 \
  --exhaustive
