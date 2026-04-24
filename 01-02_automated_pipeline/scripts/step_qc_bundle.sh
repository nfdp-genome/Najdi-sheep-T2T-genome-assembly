#!/usr/bin/env bash
# QC wave (BUSCO+QUAST) then optional dotplot vs reference
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "${DIR}/step_qc_wave.sh"
if [[ -n "${REFERENCE_FASTA:-}" && -n "${DOTPREP:-}" && -n "${DOTPLOT_OUTPREFIX:-}" ]]; then
  export ASSEMBLY_FASTA="${FA}"
  export OUTPREFIX="${DOTPLOT_OUTPREFIX}"
  bash "${DIR}/step_dotplot.sh"
fi
