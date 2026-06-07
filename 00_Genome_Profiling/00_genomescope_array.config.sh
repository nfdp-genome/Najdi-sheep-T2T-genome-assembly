#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --job-name=genomescope
#SBATCH --array=1-2
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err

set -euo pipefail

# =========================
#   CONFIG
# =========================
# Config-driven copy of 00_genomescope_array.sh: all paths/params live in
# config.yaml next to this script, so nothing site-specific is baked in.
# Override with: CONFIG=/path/to/config.yaml sbatch 00_genomescope_array.config.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${CONFIG:-${SCRIPT_DIR}/config.yaml}"
[[ -f "$CONFIG" ]] || { echo "config not found: $CONFIG" >&2; exit 1; }

# Read a top-level scalar key from the YAML config.
yaml_scalar() {
  sed -nE "s/^[[:space:]]*$1:[[:space:]]*[\"']?([^\"'#]*[^\"'# ]).*/\1/p" "$CONFIG" | head -1
}

GENOMESCOPE_R="$(yaml_scalar genomescope_r)"
OUTPUT_BASE="$(yaml_scalar output_base)"
KMER="$(yaml_scalar kmer)"
PLOIDY="$(yaml_scalar ploidy)"
KMC_MEM_GB="$(yaml_scalar kmc_mem_gb)"
KMC_MAX_COUNT="$(yaml_scalar kmc_max_count)"

# Read the samples list ("<id> <hifi_path>" per entry).
mapfile -t SAMPLES < <(awk '
  /^samples:/ {inblk=1; next}
  inblk && /^[^[:space:]]/ {inblk=0}
  inblk && /^[[:space:]]*-[[:space:]]/ {
    sub(/^[[:space:]]*-[[:space:]]*/, "")
    gsub(/^["\x27]|["\x27][[:space:]]*$/, "")
    print
  }' "$CONFIG")

IDX=$((SLURM_ARRAY_TASK_ID - 1))
[[ $IDX -ge 0 && $IDX -lt ${#SAMPLES[@]} ]] || {
  echo "array index $SLURM_ARRAY_TASK_ID out of range (${#SAMPLES[@]} samples in config)" >&2
  exit 1
}
SAMPLE_ID="${SAMPLES[$IDX]%% *}"
SAMPLE_PATH="${SAMPLES[$IDX]#* }"

# =========================
#   MODULES
# =========================
module purge
module load kmc/3.1.2
module load R

# =========================
#   OUTPUT DIRECTORY
# =========================
BASE_OUT="${OUTPUT_BASE}/${SAMPLE_ID}"
mkdir -p "${BASE_OUT}/tmp"   # tmp is required by KMC
cd "${BASE_OUT}"

# =========================
#   CREATE hifi.fofn
# =========================
echo "${SAMPLE_PATH}" > hifi.fofn

# =========================
#   K-MER COUNTING
# =========================
kmc \
  -k"${KMER}" \
  -t"${SLURM_CPUS_PER_TASK}" \
  -m"${KMC_MEM_GB}" \
  -ci1 \
  -cs"${KMC_MAX_COUNT}" \
  @hifi.fofn \
  kmcdb tmp

# =========================
#   HISTOGRAM
# =========================
kmc_tools transform kmcdb histogram sample.histo -cx"${KMC_MAX_COUNT}"

# =========================
#   GENOMESCOPE
# =========================
Rscript "${GENOMESCOPE_R}" \
  -i sample.histo \
  -k "${KMER}" \
  -p "${PLOIDY}" \
  -o genomescope_out
