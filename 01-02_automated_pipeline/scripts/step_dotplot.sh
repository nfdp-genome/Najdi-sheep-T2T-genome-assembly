#!/usr/bin/env bash
# Frozen from scripts/dotplot_references_LR.T1.5.sh (nucmer + DotPrep)
set -euo pipefail
: "${REFERENCE_FASTA:?}" "${ASSEMBLY_FASTA:?}" "${OUTPREFIX:?}" "${DOTPREP:?}" "${THREADS:?}"
mkdir -p "$(dirname "${OUTPREFIX}")"
nucmer \
  --mum \
  -l 100 \
  -c 500 \
  -t "${THREADS}" \
  "${REFERENCE_FASTA}" \
  "${ASSEMBLY_FASTA}" \
  -p "${OUTPREFIX}"

python3 "${DOTPREP}" \
  --delta "${OUTPREFIX}.delta" \
  --out   "${OUTPREFIX}"
