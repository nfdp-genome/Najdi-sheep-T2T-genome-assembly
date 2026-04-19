#!/usr/bin/env bash
# Frozen from scripts/polish_medaka_hifiasm_hic_with_ont_V1.0.sh
set -euo pipefail
: "${THREADS:?}" "${ONT_READS:?}" "${ASM_IN:?}" "${WORK_DIR:?}" "${TAG:?}" "${FINAL_OUT:?}"
MEDAKA_MODEL="${MEDAKA_MODEL:-r941_min_sup_g507}"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
if [[ ! -e "ont.fastq.gz" ]]; then
    ln -s "${ONT_READS}" "ont.fastq.gz"
fi
if [[ ! -e "${TAG}.fasta" ]]; then
    ln -s "${ASM_IN}" "${TAG}.fasta"
fi
REF_FA="${TAG}.fasta"
minimap2 -ax map-ont -t "${THREADS}" "${REF_FA}" "ont.fastq.gz" \
    | samtools sort -@ "${THREADS}" -o align.bam
samtools index align.bam
OUT_DIR="medaka_r1"
rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"
medaka_consensus \
    -i "ont.fastq.gz" \
    -d "${REF_FA}" \
    -o "${OUT_DIR}" \
    -t "${THREADS}" \
    -m "${MEDAKA_MODEL}"
[[ -s "${OUT_DIR}/consensus.fasta" ]] || { echo "[ERROR] medaka consensus missing" >&2; exit 1; }
mkdir -p "$(dirname "${FINAL_OUT}")"
cp -f "${OUT_DIR}/consensus.fasta" "${FINAL_OUT}"
