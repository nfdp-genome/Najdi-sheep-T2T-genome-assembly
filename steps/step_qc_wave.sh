#!/usr/bin/env bash
# Frozen BUSCO + QUAST from scripts/run_quast_busco_v1.1.sh
set -euo pipefail
: "${FA:?}" "${LABEL:?}" "${L1_PATH:?}" "${L2_PATH:?}" "${BUSCO_OUT_ROOT:?}" "${LOG_ROOT:?}" "${CPUS:?}"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "${BUSCO_OUT_ROOT}" "${LOG_ROOT}"

OUTNAME1="busco_${LABEL}_mammalia"
busco \
  -i "${FA}" \
  -o "${OUTNAME1}" \
  -m genome \
  -l "${L1_PATH}" \
  -c "${CPUS}" \
  --out_path "${BUSCO_OUT_ROOT}"

OUTNAME2="busco_${LABEL}_laurasiatheria"
busco \
  -i "${FA}" \
  -o "${OUTNAME2}" \
  -m genome \
  -l "${L2_PATH}" \
  -c "${CPUS}" \
  --out_path "${BUSCO_OUT_ROOT}"

QUAST_OUT_ROOT="${QUAST_OUT_ROOT:-$(dirname "${BUSCO_OUT_ROOT}")/quast}"
mkdir -p "${QUAST_OUT_ROOT}"
qdir="${QUAST_OUT_ROOT}/quast_${LABEL}"
qlog="${LOG_ROOT}/quast_${LABEL}.log"
mkdir -p "${qdir}"
quast --large -t "${CPUS}" -o "${qdir}" "${FA}" &> "${qlog}"
