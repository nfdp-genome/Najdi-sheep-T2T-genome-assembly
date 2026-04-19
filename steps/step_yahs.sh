#!/usr/bin/env bash
# Frozen argv from scripts/hiC_Yahs_scaffolding_v1.0.sh
set -euo pipefail
: "${YAHS_WORKDIR:?}" "${ASM_FASTA:?}" "${HIC_BAM:?}" "${YAHS_PREFIX:?}" "${FINAL_SCAFFOLD_FA:?}"
mkdir -p "${YAHS_WORKDIR}"
cd "${YAHS_WORKDIR}"
ASM_BASENAME="$(basename "${ASM_FASTA}")"
if [[ ! -e "${ASM_BASENAME}" ]]; then
  ln -s "${ASM_FASTA}" "${ASM_BASENAME}"
fi
yahs "${ASM_BASENAME}" "${HIC_BAM}" -o "${YAHS_PREFIX}"
FINAL_FA="${YAHS_PREFIX}_scaffolds_final.fa"
[[ -s "${FINAL_FA}" ]] || { echo "[ERROR] YaHS output missing: ${FINAL_FA}" >&2; exit 1; }
mkdir -p "$(dirname "${FINAL_SCAFFOLD_FA}")"
cp -f "${FINAL_FA}" "${FINAL_SCAFFOLD_FA}"
