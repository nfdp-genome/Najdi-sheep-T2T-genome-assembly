#!/usr/bin/env bash
# Hi-C mapping step: maps Hi-C reads to contig FASTA, produces sorted+indexed BAM
set -euo pipefail
: "${HIC_WORKDIR:?}" "${ASM_FASTA:?}" "${HIC_R1:?}" "${HIC_R2:?}" "${HIC_BAM_OUT:?}" "${THREADS:?}"

mkdir -p "${HIC_WORKDIR}"
cd "${HIC_WORKDIR}"

echo "[INFO] Indexing assembly: ${ASM_FASTA}"
if [[ ! -f "${ASM_FASTA}.bwt" ]]; then
  bwa index "${ASM_FASTA}"
else
  echo "[INFO] BWA index exists — skipping."
fi

echo "[INFO] Mapping Hi-C reads..."
bwa mem -5SP -t "${THREADS}" "${ASM_FASTA}" "${HIC_R1}" "${HIC_R2}" \
  | samtools view -@ "${THREADS}" -b - \
  | samtools sort -@ "${THREADS}" -o "${HIC_BAM_OUT}" -

[[ -s "${HIC_BAM_OUT}" ]] || { echo "[ERROR] Empty sorted BAM"; exit 1; }

echo "[INFO] Indexing BAM..."
samtools index "${HIC_BAM_OUT}"

[[ -f "${HIC_BAM_OUT}.bai" ]] || { echo "[ERROR] Missing BAM index"; exit 1; }

echo "[INFO] DONE: ${HIC_BAM_OUT}"
