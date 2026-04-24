#!/usr/bin/env bash
# Frozen from scripts/polish_racon_hifiasm_hic_with_ont_V1.0.sh (ROUNDS=2)
set -euo pipefail
: "${THREADS:?}" "${ONT_READS:?}" "${ASM_IN:?}" "${WORK_DIR:?}" "${TAG:?}" "${FINAL_OUT:?}"
ROUNDS=2
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
if [[ ! -e "ont.fastq.gz" ]]; then
    ln -s "${ONT_READS}" "ont.fastq.gz"
fi
if [[ ! -e "${TAG}.fasta" ]]; then
    ln -s "${ASM_IN}" "${TAG}.fasta"
fi
CURRENT_ASM="${TAG}.fasta"
ROUND=1
while (( ROUND <= ROUNDS )); do
    PAF_OUT="round${ROUND}.paf"
    ASM_OUT="${TAG}_racon_r${ROUND}.fasta"
    minimap2 -x map-ont -t "${THREADS}" "${CURRENT_ASM}" "ont.fastq.gz" > "${PAF_OUT}"
    racon -t "${THREADS}" "ont.fastq.gz" "${PAF_OUT}" "${CURRENT_ASM}" > "${ASM_OUT}"
    CURRENT_ASM="${ASM_OUT}"
    ROUND=$((ROUND + 1))
done
mkdir -p "$(dirname "${FINAL_OUT}")"
cp -f "${CURRENT_ASM}" "${FINAL_OUT}"
