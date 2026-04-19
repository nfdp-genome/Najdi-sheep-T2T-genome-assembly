#!/usr/bin/env bash
# GFA→FASTA awk from scripts/hifiasm_assembly_quality_v1.3.sh (gfa_to_fa core only)
set -euo pipefail
: "${GFA_IN:?}" "${FASTA_OUT:?}"
mkdir -p "$(dirname "${FASTA_OUT}")"
[[ -s "${GFA_IN}" ]] || { echo "[ERROR] Missing GFA: ${GFA_IN}" >&2; exit 1; }
awk '$1=="S"{print ">"$2"\n"$3}' "${GFA_IN}" > "${FASTA_OUT}.tmp"
mv -f "${FASTA_OUT}.tmp" "${FASTA_OUT}"
