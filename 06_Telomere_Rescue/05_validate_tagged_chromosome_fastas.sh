#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 05_validate_tagged_chromosome_fastas.sh
#
# QC pass 1: expected-vs-observed placeholder lists (must be identical).
# ------------------------------------------------------------------------------

# Quality check for TELROI-tagged chromosome FASTA files

set -euo pipefail

BASE_DIR="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"

cd "${BASE_DIR}"

echo "### 1) Tagged FASTA count vs unique sample–chromosome count"

tagged_file_count=$(
    find tagged_chromosomes \
        -maxdepth 1 \
        -type f \
        -name '*_TELROI_TAGGED.fasta' |
    wc -l
)

expected_file_count=$(
    cut -f1,2 telomere_roi_table.tsv |
    tail -n +2 |
    sort -u |
    wc -l
)

echo "Tagged FASTA files: ${tagged_file_count}"
echo "Expected chromosomes: ${expected_file_count}"

if [[ "${tagged_file_count}" -ne "${expected_file_count}" ]]; then
    echo "[FAIL] Tagged FASTA count does not match the expected chromosome count."
    exit 1
fi

echo "[PASS] Tagged FASTA count matches the expected chromosome count."


echo
echo "### 2) Check FASTA headers"

grep -hE '^>' tagged_chromosomes/*_TELROI_TAGGED.fasta |
awk 'NR <= 10 {print}'

invalid_header_count=$(
    grep -hE '^>' tagged_chromosomes/*_TELROI_TAGGED.fasta |
    grep -vc 'TELROI_TAGGED' || true
)

if [[ "${invalid_header_count}" -ne 0 ]]; then
    echo "[FAIL] One or more FASTA headers do not contain TELROI_TAGGED."
    exit 1
fi

echo "[PASS] FASTA headers look valid."


echo
echo "### 3) Count expected placeholders"

test -s expected_placeholders.txt

expected_total=$(
    grep -cve '^[[:space:]]*$' expected_placeholders.txt
)

echo "Expected TELROI placeholders: ${expected_total}"


echo
echo "### 4) Extract observed TELROI placeholders"

grep -RhoE '@@TELROI\|[^@]*@@' tagged_chromosomes |
sort -u > observed_tags.txt

sort -u expected_placeholders.txt > expected_placeholders.sorted.txt

observed_total=$(wc -l < observed_tags.txt)
expected_unique_total=$(wc -l < expected_placeholders.sorted.txt)

echo "Expected unique placeholders: ${expected_unique_total}"
echo "Observed unique placeholders: ${observed_total}"


echo
echo "### 5) Compare expected and observed placeholders"

if diff -u expected_placeholders.sorted.txt observed_tags.txt; then
    echo "[PASS] All expected TELROI placeholders were found."
else
    echo "[FAIL] Expected and observed TELROI placeholders differ."
    exit 1
fi
