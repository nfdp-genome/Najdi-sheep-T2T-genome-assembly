#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 01_phase0_declare_inputs.py
#
# PHASE 0 - declare the input paths per sample (scaffolded chromosome FASTA + gap coordinates
# table). Reference cell before running any of the phases below.
# ------------------------------------------------------------------------------
# PHASE 0 — Inputs
    # genome with scaffolded contigs 
        # /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta
        # /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta

# working directory: /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling
# gaps informations 
| Sample | Chr | Final scaffold (bp) | Gap position (scaffold coordinates) |
| ------ | --- | ------------------: | ----------------------------------- |
| 1271   | 08  |          93,575,207 | 76,483,492–76,483,591               |
| 1271   | 17  |          75,039,272 | 39,008,832–39,008,931               |
| 1268   | 10  |          87,857,919 | 37,430,396–37,430,495               |

| Sample | Chr | Final scaffold (bp) | Gap position (bp)           |
|--------|-----|---------------------|-----------------------------|
| 1268   | 02  |        255,348,196    | 138,435,449–138,435,836     |
