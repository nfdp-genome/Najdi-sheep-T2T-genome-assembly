#!/usr/bin/env bash
# Frozen hifiasm argv from scripts/hifiasm_assembly_telomere_flag_noDual-scaf_v1.0.sh
set -euo pipefail
: "${HIFIASM_OUTDIR:?}" "${HIFIASM_PREFIX:?}" "${THREADS:?}" "${HIFIASM_HIFI:?}" "${HIFIASM_ONT:?}"
mkdir -p "${HIFIASM_OUTDIR}"
cd "${HIFIASM_OUTDIR}"
LOG="${HIFIASM_PREFIX}.hifiasm.log"
hifiasm \
  -o "${HIFIASM_PREFIX}.asm" \
  -t "${THREADS}" \
  --telo-m TTAGGG \
  --ul "${HIFIASM_ONT}" \
  "${HIFIASM_HIFI}" \
  2> "${LOG}"
