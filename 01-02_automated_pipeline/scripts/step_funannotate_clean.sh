#!/usr/bin/env bash
# Frozen from scripts/funannotate_clean_V3.sh (runner prefix configurable)
set -euo pipefail
: "${GENOME:?}" "${OUT_FASTA:?}"
FUNANNOTATE_RUNNER="${FUNANNOTATE_RUNNER:-micromamba run -n funannotate_clean_env}"
mkdir -p "$(dirname "${OUT_FASTA}")"
# shellcheck disable=SC2086
eval "${FUNANNOTATE_RUNNER}" funannotate clean \
  -i "${GENOME}" \
  -o "${OUT_FASTA}" \
  --minlen 15000000 \
  --pident 95 \
  --cov 95 \
  --exhaustive
