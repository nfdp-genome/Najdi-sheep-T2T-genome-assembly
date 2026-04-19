#!/usr/bin/env bash
# Frozen from scripts/gapfilling_hifiasm_hic_ont_tgsgapcloser_v1.1.sh
set -euo pipefail
: "${THREADS:?}" "${ASM_IN:?}" "${ONT_FASTQ:?}" "${WORK_DIR:?}" "${TAG:?}" "${FINAL_OUT:?}"
TGSGC_BIN="$(command -v tgsgapcloser)"
RACON_BIN="$(command -v racon)"
mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"
if [[ ! -e "${TAG}.fasta" ]]; then
  ln -s "${ASM_IN}" "${TAG}.fasta"
fi
if [[ ! -e "ont.fastq.gz" ]]; then
  ln -s "${ONT_FASTQ}" "ont.fastq.gz"
fi
ONT_FASTA="ont_reads.fasta"
if [[ ! -s "${ONT_FASTA}" ]]; then
  zcat "ont.fastq.gz" | awk 'NR%4==1{print ">" substr($0,2)} NR%4==2{print}' > "${ONT_FASTA}"
fi
OUT_PREFIX="${TAG}_tgsgapcloser"
"${TGSGC_BIN}" \
  --scaff  "${TAG}.fasta" \
  --reads  "${ONT_FASTA}" \
  --output "${OUT_PREFIX}" \
  --racon  "${RACON_BIN}" \
  --thread "${THREADS}" \
  > "${OUT_PREFIX}.log" 2> "${OUT_PREFIX}.err"
FINAL1="${OUT_PREFIX}.scaff_seq"
FINAL2="${OUT_PREFIX}.contig"
if [[ -s "${FINAL1}" ]]; then
  FINAL="${FINAL1}"
elif [[ -s "${FINAL2}" ]]; then
  FINAL="${FINAL2}"
else
  echo "[ERROR] tgsgapcloser output missing" >&2
  exit 1
fi
mkdir -p "$(dirname "${FINAL_OUT}")"
cp -f "${FINAL}" "${FINAL_OUT}"
