#!/usr/bin/env bash
# Frozen SALSA + bedtools from scripts/hiC_Salsa_scaffolding_v1.1.sh
set -euo pipefail
: "${SALSA_WORKDIR:?}" "${REF_FASTA:?}" "${SORTED_BAM:?}" "${PREFIX:?}" "${SALSA_RUN:?}" "${FINAL_SCAFFOLD_FA:?}"
THREADS="${THREADS:-40}"
mkdir -p "${SALSA_WORKDIR}"
cd "${SALSA_WORKDIR}"

ASM_FASTA_BASENAME="$(basename "${REF_FASTA}")"
if [[ ! -e "${ASM_FASTA_BASENAME}" ]]; then
  ln -s "${REF_FASTA}" "${ASM_FASTA_BASENAME}"
fi
if [[ ! -e "${ASM_FASTA_BASENAME}.fai" ]]; then
  ln -s "${REF_FASTA}.fai" "${ASM_FASTA_BASENAME}.fai"
fi

LENGTHS_FILE="${ASM_FASTA_BASENAME%.fasta}.contig_lengths.txt"
if [[ ! -f "${LENGTHS_FILE}" ]]; then
  cut -f1,2 "${ASM_FASTA_BASENAME}.fai" > "${LENGTHS_FILE}"
fi

BED_RAW="${PREFIX}.bed"
BED_SORTED="${PREFIX}.namesorted.bed"
if [[ ! -s "${BED_SORTED}" ]]; then
  bedtools bamtobed -i "${SORTED_BAM}" > "${BED_RAW}"
  sort -k4,4 "${BED_RAW}" > "${BED_SORTED}"
  rm -f "${BED_RAW}"
fi

OUTDIR="salsa_out"
mkdir -p "${OUTDIR}"

"${SALSA_RUN}" \
    -a "${ASM_FASTA_BASENAME}" \
    -l "${LENGTHS_FILE}" \
    -b "${BED_SORTED}" \
    -e DNASE \
    -i 3 \
    -m yes \
    -o "${OUTDIR}"

FINAL_FASTA=""
if [[ -s "${OUTDIR}/scaffolds_FINAL.fasta" ]]; then
    FINAL_FASTA="${OUTDIR}/scaffolds_FINAL.fasta"
elif [[ -s "${OUTDIR}/scaffolds.fasta" ]]; then
    FINAL_FASTA="${OUTDIR}/scaffolds.fasta"
fi
[[ -n "${FINAL_FASTA}" ]] || { echo "[ERROR] SALSA produced no scaffolds.fasta" >&2; exit 1; }
mkdir -p "$(dirname "${FINAL_SCAFFOLD_FA}")"
cp -f "${FINAL_FASTA}" "${FINAL_SCAFFOLD_FA}"
