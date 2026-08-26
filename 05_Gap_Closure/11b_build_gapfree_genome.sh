#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 11b_build_gapfree_genome.sh
#
# Assemble the final gap-free per-sample chromosome FASTAs into one whole-genome FASTA per
# sample (input to the telomere-rescue workflow in 06_Telomere_Rescue/).
# ------------------------------------------------------------------------------
# Building the gapfree genome 
cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq 

grep chr02 /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta

grep chr10 /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta


grep chr08 /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta

grep chr17 /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta

pwd 

ls 1268_chr02.gap_filled.fasta
ls 1268_chr10_scaffolded_HuM.gap_filled.fasta
ls 1271_chr08_scaffolded_HuF.gap_filled.fasta
ls 1271_chr17_scaffolded_HuF.no_gap.fasta


grep ">"  1268_chr02.gap_filled.fasta
grep ">"  1268_chr10_scaffolded_HuM.gap_filled.fasta
grep ">"  1271_chr08_scaffolded_HuF.gap_filled.fasta
grep ">"  1271_chr17_scaffolded_HuF.no_gap.fasta
