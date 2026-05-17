#!/usr/bin/env bash
#SBATCH --job-name=step_nfcore_hic_contactmap
#SBATCH --partition=batch
#SBATCH --constraint=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=96:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

MODULES="nextflow singularity"

CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      [[ $# -ge 2 ]] || { echo "[ERROR] --config requires a YAML path" >&2; exit 1; }
      CONFIG="$2"; shift 2
      ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
  esac
done

[[ -n "${CONFIG}" && -s "${CONFIG}" ]] || { echo "[ERROR] Missing --config YAML" >&2; exit 1; }

python3 -c "import yaml" >/dev/null 2>&1 || { echo "[ERROR] PyYAML missing" >&2; exit 1; }

eval "$(python3 - "$CONFIG" <<'PYEOF'
import sys, yaml, shlex
cfg = yaml.safe_load(open(sys.argv[1])) or {}
params = cfg.get('params', {}) or {}
for key, value in params.items():
    if isinstance(value, list):
        value = '|'.join(str(x) for x in value)
    elif value is None:
        value = ''
    else:
        value = str(value)
    print(f"{key}={shlex.quote(value)}")
print('log_file=' + shlex.quote(str(cfg.get('log_file', './logs/step_nfcore_hic_contactmap.log'))))
PYEOF
)"

mkdir -p "$(dirname "$log_file")"
exec > >(tee -a "$log_file") 2>&1

if [[ -n "${MODULES}" ]]; then
  source /etc/profile.d/modules.sh 2>/dev/null || true
  if declare -F module >/dev/null 2>&1 || command -v module >/dev/null 2>&1; then
    module purge || true
    for module_name in ${MODULES}; do module load "$module_name"; done
  else
    echo "[WARN] Environment modules are unavailable; expecting tools on PATH."
  fi
fi

for required in TAG SAMPLE FASTA HIC_R1 HIC_R2 OUT_BASE WORK_BASE; do
  [[ -n "${!required:-}" ]] || { echo "[ERROR] Missing required config param: ${required}" >&2; exit 1; }
done

[[ -f "${FASTA}" ]] || { echo "[ERROR] Missing FASTA: ${FASTA}" >&2; exit 1; }
[[ -f "${HIC_R1}" ]] || { echo "[ERROR] Missing HIC_R1: ${HIC_R1}" >&2; exit 1; }
[[ -f "${HIC_R2}" ]] || { echo "[ERROR] Missing HIC_R2: ${HIC_R2}" >&2; exit 1; }

NXFCORE_HIC_REV="${NXFCORE_HIC_REV:-2.1.0}"
NXFCORE_PROFILE="${NXFCORE_PROFILE:-kaust}"
MIN_CIS_DIST="${MIN_CIS_DIST:-1000}"
SPLIT_FASTQ="${SPLIT_FASTQ:-true}"
RUN_HICTK_CONVERT="${RUN_HICTK_CONVERT:-true}"
DNASE="${DNASE:-true}"

LABEL="${SAMPLE}_${TAG}"
RUN_DIR="${OUT_BASE%/}/${LABEL}"
WORK_DIR="${WORK_BASE%/}/${LABEL}/work"
RESULTS_DIR="${RUN_DIR}/results"
LOG_DIR="${RUN_DIR}/logs"
SAMPLESHEET="${RUN_DIR}/samplesheet_${LABEL}.csv"
SINGULARITY_CACHE="${SINGULARITY_CACHE:-${WORK_BASE%/}/singularity_cache}"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${RESULTS_DIR}" "${LOG_DIR}" "${SINGULARITY_CACHE}"

export NXF_HOME="${WORK_BASE%/}/${LABEL}/.nextflow"
export NXF_WORK="${WORK_DIR}"
export NXF_SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"
export NXF_ANSI_LOG=false

cat > "${SAMPLESHEET}" <<EOF
sample,fastq_1,fastq_2
${SAMPLE},${HIC_R1},${HIC_R2}
EOF

NFCORE_ARGS=(
  run nf-core/hic
  -r "${NXFCORE_HIC_REV}"
  -profile "${NXFCORE_PROFILE}"
  --input "${SAMPLESHEET}"
  --fasta "${FASTA}"
  --min_cis_dist "${MIN_CIS_DIST}"
  --outdir "${RESULTS_DIR}"
  -work-dir "${WORK_DIR}"
  -with-report "${LOG_DIR}/${LABEL}_report.html"
  -with-timeline "${LOG_DIR}/${LABEL}_timeline.html"
  -with-trace "${LOG_DIR}/${LABEL}_trace.txt"
)

if [[ "${DNASE}" == "true" ]]; then
  NFCORE_ARGS+=(--dnase)
fi

if [[ "${SPLIT_FASTQ}" == "true" ]]; then
  NFCORE_ARGS+=(--split_fastq)
fi

if [[ "${RESUME:-false}" == "true" ]]; then
  NFCORE_ARGS+=(-resume)
fi

echo "[INFO] Running nf-core/hic for ${LABEL}"
echo "[INFO] FASTA: ${FASTA}"
echo "[INFO] Results: ${RESULTS_DIR}"
echo "[INFO] Work: ${WORK_DIR}"

nextflow "${NFCORE_ARGS[@]}"

if [[ "${RUN_HICTK_CONVERT}" != "true" ]]; then
  echo "[INFO] hictk conversion disabled."
  echo "[DONE] ${LABEL}"
  exit 0
fi

COOL_DIR="${RESULTS_DIR}/contact_maps/cool"
MCOOL="${COOL_DIR}/${SAMPLE}.mcool"
HIC="${COOL_DIR}/${SAMPLE}.hic"
HICTK_SIF="${HICTK_SIF:-./containers/hictk_2.1.2.sif}"
HICTK_IMAGE="${HICTK_IMAGE:-docker://ghcr.io/paulsengroup/hictk:2.1.2}"
TMPDIR_HICTK="${TMPDIR_HICTK:-${WORK_BASE%/}/tmp_hictk_${LABEL}_${SLURM_JOB_ID:-manual}}"

if [[ ! -f "${MCOOL}" ]]; then
  echo "[WARN] Expected mcool not found: ${MCOOL}"
  echo "[WARN] Skipping hictk conversion."
  echo "[DONE] ${LABEL}"
  exit 0
fi

mkdir -p "$(dirname "${HICTK_SIF}")" "${TMPDIR_HICTK}"

if [[ ! -f "${HICTK_SIF}" ]]; then
  echo "[INFO] Pulling hictk container: ${HICTK_IMAGE}"
  singularity pull "${HICTK_SIF}" "${HICTK_IMAGE}"
fi

echo "[INFO] Converting ${MCOOL} to ${HIC}"
singularity run --cleanenv \
  -B "${COOL_DIR}:/data" \
  -B "${TMPDIR_HICTK}:/tmpdir" \
  "${HICTK_SIF}" \
  convert "/data/${SAMPLE}.mcool" "/data/${SAMPLE}.hic" \
  --threads "${SLURM_CPUS_PER_TASK:-32}" \
  --tmpdir /tmpdir \
  --force

echo "[DONE] ${LABEL}"
echo "Run dir : ${RUN_DIR}"
echo "Outdir  : ${RESULTS_DIR}"
echo "Hi-C    : ${HIC}"
