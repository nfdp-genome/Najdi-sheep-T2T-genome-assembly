#!/usr/bin/env bash
# Frozen from scripts/gapfilling_hifiasm_hic_ont_ntlink_v1.0.sh
set -euo pipefail
: "${THREADS:?}" "${NTLINK_BIN:?}" "${ASM_IN:?}" "${ONT_FASTQ:?}" "${WORK_DIR:?}" "${TAG:?}" "${POLISHER:?}" "${FINAL_OUT:?}"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
TARGET_FA="${TAG}_${POLISHER}.fasta"
READS_FAQ="ont.fastq.gz"
ln -sf "${ASM_IN}" "${TARGET_FA}"
ln -sf "${ONT_FASTQ}" "${READS_FAQ}"
PREFIX="${TAG}_${POLISHER}_ntlink"
"${NTLINK_BIN}" scaffold gap_fill \
    target="${TARGET_FA}" \
    reads="${READS_FAQ}" \
    prefix="${PREFIX}" \
    t="${THREADS}" \
    k=32 \
    w=250 \
    > "${PREFIX}.ntlink.log" 2> "${PREFIX}.ntlink.err"
FINAL_CANDIDATE=$(ls -1 "${PREFIX}"*.ntLink.scaffolds.fa 2>/dev/null | head -n 1 || true)
[[ -n "${FINAL_CANDIDATE}" && -s "${FINAL_CANDIDATE}" ]] || { echo "[ERROR] ntLink output missing for ${PREFIX}" >&2; exit 1; }
mkdir -p "$(dirname "${FINAL_OUT}")"
cp -f "${FINAL_CANDIDATE}" "${FINAL_OUT}"
