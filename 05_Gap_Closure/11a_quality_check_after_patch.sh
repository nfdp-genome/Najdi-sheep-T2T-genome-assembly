#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 11a_quality_check_after_patch.sh
#
# Quality check on the four patched chromosomes: total length, gap count,
# N's per 100 kbp.
# ------------------------------------------------------------------------------
# quality check:
cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq

#echo -n GCATGGGGACTGGCCTTTCCTGAGGCCACCAGAGCGGGTCCCTGAGGTCCCCGTCGTAAGTCGAGAGCACCTGCAGCATACTCGACTAGTGA | wc -c 
#grep -v "^>" 1268_chr02.gap_filled.fasta | grep -o 'N\{388\}' | wc -l
#grep -v "^>" 1268_chr02.fasta | grep -o 'N\{388\}' | wc -l
#stat -c %s 1268_chr02.fasta
#stat -c %s 1268_chr02.gap_filled.fasta
#echo $(( $(stat -c %s 1268_chr02.fasta) - $(stat -c %s 1268_chr02.gap_filled.fasta) ))
#grep -v "^>" 1268_chr02.fasta | grep -o 'N' | wc -l
#grep -v "^>" 1268_chr02.gap_filled.fasta| grep -o 'N' | wc -l
echo "Ns original oneline:"
grep -v "^>" 1268_chr02.oneline.fasta | tr -cd 'N' | wc -c

echo "Ns gapfilled oneline:"
grep -v "^>" 1268_chr02.gap_filled.oneline.fasta | tr -cd 'N' | wc -c
