#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 10_patch_1268_chr10.sh
#
# Patch the selected gap-filling contig into chromosome 10 of sample 1268
# (the multi-contig chromosome joined earlier by RagTag).
# ------------------------------------------------------------------------------
# 1268_chr10 gap 
## this is a true gap .. 
cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq

grep -v "^>" 1268_chr10_scaffolded_HuM.fasta | tr -cd 'N' | wc -c
sed 's/N\{100\}/caattaaata/g' 1268_chr10_scaffolded_HuM.fasta > 1268_chr10_scaffolded_HuM.gap_filled.fasta
grep -v "^>" 1268_chr10_scaffolded_HuM.gap_filled.fasta | tr -cd 'N' | wc -c
