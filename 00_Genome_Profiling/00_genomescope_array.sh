#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --job-name=genomescope
#SBATCH --array=1-2
#SBATCH --output=/ibex/project/c2293/transit/mashael/t2t/logs/%x.%A_%a.out
#SBATCH --error=/ibex/project/c2293/transit/mashael/t2t/logs/%x.%A_%a.err

set -euo pipefail

# =========================
#   MODULES
# =========================
module purge
module load kmc/3.1.2
module load R

# =========================
#   INPUT FILES
# =========================
SAMPLES=(
"/ibex/project/c2293/transit/Najdi_T2T/input/1271/hifi/hifi.fastq.gz"
"/ibex/project/c2293/transit/Najdi_T2T/input/1268/hifi/hifi.fastq.gz"
)

SAMPLE_PATH=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
SAMPLE_ID=$(basename $(dirname $(dirname ${SAMPLE_PATH})))   # 1271 or 1268

# =========================
#   OUTPUT DIRECTORIES
# =========================
BASE_OUT="/ibex/project/c2293/transit/mashael/t2t/genomescope/${SAMPLE_ID}"
mkdir -p "${BASE_OUT}"
cd "${BASE_OUT}"

# REQUIRED for KMC
mkdir -p tmp

# =========================
#   CREATE hifi.fofn
# =========================
echo "${SAMPLE_PATH}" > hifi.fofn

# =========================
#   K-MER COUNTING
# =========================
kmc \
  -k21 \
  -t${SLURM_CPUS_PER_TASK} \
  -m64 \
  -ci1 \
  -cs10000000 \
  @hifi.fofn \
  kmcdb tmp

# =========================
#   HISTOGRAM
# =========================
kmc_tools transform kmcdb histogram sample.histo -cx10000000

# =========================
#   GENOMESCOPE
# =========================
Rscript \
/ibex/project/c2293/transit/mashael/t2t/tools/genomescope2.0/genomescope.R \
-i sample.histo \
-k 21 \
-p 2 \
-o genomescope_out
