#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 08_validate_extracted_telroi_fastas.sh
#
# QC of extracted TELROI FASTAs: file count, header format, absence of internal tags,
# DNA-only content.
# ------------------------------------------------------------------------------

# STEP 5: Validate extracted TELROI FASTA files

set -euo pipefail

BASE_DIR="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
ROI_TABLE="${BASE_DIR}/telomere_roi_table.tsv"
ROI_DIR="${BASE_DIR}/roi_sequences"

cd "${BASE_DIR}"

test -s "${ROI_TABLE}"
test -d "${ROI_DIR}"

echo "### 1) Expected TELROIs vs extracted FASTA files"

expected_count=$(
    tail -n +2 "${ROI_TABLE}" |
    grep -cve '^[[:space:]]*$'
)

observed_count=$(
    find "${ROI_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*_TELROI.fa' |
    wc -l
)

echo "Expected TELROIs: ${expected_count}"
echo "Extracted FASTAs: ${observed_count}"

if [[ "${expected_count}" -ne "${observed_count}" ]]; then
    echo "[FAIL] Extracted FASTA count does not match the TELROI table."
    exit 1
fi

echo "[PASS] Extracted FASTA count matches the TELROI table."


echo
echo "### 2) Example filenames"

find "${ROI_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*_TELROI.fa' \
    -printf '%f\n' |
sort |
head -n 10


echo
echo "### 3) Validate FASTA headers"

invalid_headers=0

while IFS= read -r fasta; do
    header=$(head -n 1 "${fasta}")

    if [[ ! "${header}" =~ ^\>@@TELROI\|.*@@$ ]]; then
        echo "[INVALID_HEADER] ${fasta}: ${header}"
        invalid_headers=$((invalid_headers + 1))
    fi
done < <(
    find "${ROI_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*_TELROI.fa' |
    sort
)

if [[ "${invalid_headers}" -ne 0 ]]; then
    echo "[FAIL] Invalid TELROI FASTA headers: ${invalid_headers}"
    exit 1
fi

echo "[PASS] All FASTA headers contain valid TELROI tags."


echo
echo "### 4) Confirm TELROI tags are absent from sequence lines"

internal_tag_count=$(
    grep -h '@@TELROI' "${ROI_DIR}"/*_TELROI.fa |
    grep -vc '^>' || true
)

if [[ "${internal_tag_count}" -ne 0 ]]; then
    echo "[FAIL] TELROI tags were found inside sequence lines."
    exit 1
fi

echo "[PASS] TELROI tags occur only in FASTA headers."


echo
echo "### 5) Validate sequence characters"

invalid_sequence_lines=$(
    grep -hEv '^>|^[ACGTNacgtn]+$|^[[:space:]]*$' \
        "${ROI_DIR}"/*_TELROI.fa |
    wc -l
)

if [[ "${invalid_sequence_lines}" -ne 0 ]]; then
    echo "[FAIL] Invalid sequence lines detected: ${invalid_sequence_lines}"
    exit 1
fi

echo "[PASS] All sequences contain only A, C, G, T or N."


echo
echo "### 6) Example 5′ TELROI"

example_5p=$(
    find "${ROI_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*_5p_TELROI.fa' |
    sort |
    head -n 1
)

if [[ -n "${example_5p}" ]]; then
    echo "File: ${example_5p}"
    head -n 2 "${example_5p}"
else
    echo "[INFO] No 5′ TELROI FASTA files found."
fi


echo
echo "### 7) Example 3′ TELROI"

example_3p=$(
    find "${ROI_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*_3p_TELROI.fa' |
    sort |
    head -n 1
)

if [[ -n "${example_3p}" ]]; then
    echo "File: ${example_3p}"
    head -n 1 "${example_3p}"
    tail -n 1 "${example_3p}"
else
    echo "[INFO] No 3′ TELROI FASTA files found."
fi


echo
echo "[PASS] TELROI FASTA sanity checks completed successfully."
