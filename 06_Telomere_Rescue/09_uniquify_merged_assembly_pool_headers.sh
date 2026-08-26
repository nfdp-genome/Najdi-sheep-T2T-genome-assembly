#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 09_uniquify_merged_assembly_pool_headers.sh
#
# Rewrite the merged per-sample assembly-pool FASTA deflines with a sequential ASMxxxx___
# prefix so all sequence IDs are unique.
# ------------------------------------------------------------------------------

# Add unique identifiers to all merged-assembly FASTA headers

set -euo pipefail

WORK_DIR="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"

IN_1268="${WORK_DIR}/1268_all_assemblies.merged.fasta"
IN_1271="${WORK_DIR}/1271_all_assemblies.merged.fasta"

OUT_1268="${WORK_DIR}/1268_all_assemblies.merged.unique.fasta"
OUT_1271="${WORK_DIR}/1271_all_assemblies.merged.unique.fasta"


# -------------------------------
# Add a sequential ASM prefix
# -------------------------------
add_unique_prefix() {
    local input_fasta="$1"
    local output_fasta="$2"

    if [[ ! -s "${input_fasta}" ]]; then
        echo "[ERROR] Input FASTA is missing or empty: ${input_fasta}" >&2
        exit 1
    fi

    awk '
        BEGIN {
            sequence_number = 0
        }

        /^>/ {
            sequence_number++
            asm_id = sprintf("ASM%03d", sequence_number)
            original_header = substr($0, 2)

            print ">" asm_id "___" original_header
            next
        }

        {
            print
        }
    ' "${input_fasta}" > "${output_fasta}"
}


add_unique_prefix "${IN_1268}" "${OUT_1268}"
add_unique_prefix "${IN_1271}" "${OUT_1271}"


echo "[OK] Unique merged-assembly FASTA files written:"
echo " - ${OUT_1268}"
echo " - ${OUT_1271}"
