#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# validate_v1_by_ONT_mapping.sh
#
# ONT-based validation of shipped v1.0 genomes: map v1.0 vs ONT reads with minimap2 -x map-ont,
# produce PAF, add a header row, and prepare the filtered manifest of terminal, extension-bearing
# reads (is_terminal == 1, tel_ext_bp >= 60, alnlen >= 120).
# ------------------------------------------------------------------------------
# mapping on ONT > output paf 
#cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/testing_telomeres_tools/hu_paper_approach/stage7d_telomerevalidation_minimap_on_v1.0
cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT

minimap2 -x map-ont -t 88 \
  /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/Najdi_T2T_NLFDP1268_v1.0/Najdi_T2T_NLFDP1268_v1.0_genome.fasta \
  /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/raw/1268/ont/ont.fastq.gz \
  > /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/Najdi_T2T_NLFDP1268_v1.0_genome_vs_ONT.paf \
  2> >(tee /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/minimap2.log >&2)

minimap2 -x map-ont -t 40 \
  /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/Najdi_T2T_NLFDP1271_v1.0/Najdi_T2T_NLFDP1271_v1.0_genome.fasta \
  /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/raw/1271/ont/ont.fastq.gz \
  > /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/Najdi_T2T_NLFDP1271_v1.0_genome_vs_ONT.paf \
  2> >(tee /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/minimap2.log >&2)

tail /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/minimap2.log
tail /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/minimap2.log

# addling header to paf: 
cut -f1-12 \
/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/Najdi_T2T_NLFDP1268_v1.0_genome_vs_ONT.paf \
> q

cut -f1-12 \
/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/Najdi_T2T_NLFDP1271_v1.0_genome_vs_ONT.paf \
> /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/Najdi_T2T_NLFDP1271_v1.0_genome_vs_ONT.clean.paf

# 1268
echo -e "qname\tqlen\tqstart\tqend\tstrand\ttname\ttlen\ttstart\ttend\tnmatch\talnlen\tmapq" > tmp && \
cat /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/Najdi_T2T_NLFDP1268_v1.0_genome_vs_ONT.clean.paf >> tmp && \
mv tmp /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1268_v1.0_genome/Najdi_T2T_NLFDP1268_v1.0_genome_vs_ONT.clean.paf


# 1271
echo -e "qname\tqlen\tqstart\tqend\tstrand\ttname\ttlen\ttstart\ttend\tnmatch\talnlen\tmapq" > tmp && \
cat /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/Najdi_T2T_NLFDP1271_v1.0_genome_vs_ONT.clean.paf >> tmp && \
mv tmp /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/ready_to_ship_sequences/map_on_ONT/Najdi_T2T_NLFDP1271_v1.0_genome/Najdi_T2T_NLFDP1271_v1.0_genome_vs_ONT.clean.paf

## calculating () 
  #qcov = (qend - qstart) / qlen 
  #is_5p = (tstart <= 2000) 
  #is_3p = ((tlen - tend) <= 2000) 
  #is_terminal = (is_5p OR is_3p) # Make a filtered manifest
Keep only:
- is_terminal == 1
- tel_ext_bp >= 60
- alnlen >= 120
