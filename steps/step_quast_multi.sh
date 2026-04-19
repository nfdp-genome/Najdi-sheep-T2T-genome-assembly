#!/usr/bin/env bash
# Multi-sample QUAST — genome paths from QUAST_FASTA_LIST (one per line)
set -euo pipefail
: "${OUT_DIR:?}" "${THREADS:?}" "${QUAST_FASTA_LIST:?}"
mkdir -p "${OUT_DIR}"
mapfile -t GENOMES < "${QUAST_FASTA_LIST}"
[[ ${#GENOMES[@]} -gt 0 ]] || { echo "[ERROR] empty QUAST_FASTA_LIST" >&2; exit 1; }

cmd=(quast.py "${GENOMES[@]}" --large --threads "${THREADS}" --output-dir "${OUT_DIR}")

if [[ -n "${QUAST_LABELS_FILE:-}" && -s "${QUAST_LABELS_FILE}" ]]; then
  mapfile -t LABELS < "${QUAST_LABELS_FILE}"
  if [[ ${#LABELS[@]} -eq ${#GENOMES[@]} ]]; then
    cmd+=(--labels "${LABELS[@]}")
  else
    echo "[WARN] label count != genome count; omitting --labels" >&2
  fi
fi

"${cmd[@]}"
