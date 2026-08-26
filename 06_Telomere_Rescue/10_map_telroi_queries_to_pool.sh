#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 10_map_telroi_queries_to_pool.sh
#
# Wrapper: submits 11_mapping_assemblies.V1.12.sh (sibling script in this folder) as a SLURM
# array over all TELROI query FASTAs for each sample. Produces one PAF per TELROI query in the
# per-sample mapping output directory.
# ------------------------------------------------------------------------------

set -euo pipefail

WORK_DIR="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"

MAPPING_SCRIPT="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/code/scripts/mapping_assemblies.V1.12.sh"

cd "${WORK_DIR}"

if [[ ! -s "${MAPPING_SCRIPT}" ]]; then
    echo "[ERROR] Mapping script is missing or empty: ${MAPPING_SCRIPT}" >&2
    exit 1
fi

sbatch "${MAPPING_SCRIPT}" bam
