#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 09_patch_1268_chr02.sh
#
# Patch the selected gap-filling contig into chromosome 02 of sample 1268
# (338 bp true gap on the anchored chromosome).
# ------------------------------------------------------------------------------
# 1268_chr02 gap 
## this is a true gap .. 
cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq
grep -v "^>" 1268_chr02.fasta | tr -cd 'N' | wc -c
    #GCATGGGGACTGGCCTTTCCTGAGGCCACCAGAGCGGGTCCCTGAGGTCCCCGTCGTAAGTCGAGAGCACCTGCAGCATACTCGACTAGTGA
grep -v "^>" 1268_chr02.fasta | grep -o 'N\{388\}' | wc -l
grep -v "^>" 1268_chr02.fasta | grep -o 'N\{389,\}' | wc -l

# convert to oneline fasta:
#awk 'BEGIN{seq=""}
#     /^>/{hdr=$0; next}
#     {gsub(/[ \t\r\n]/,""); seq=seq $0}
#     END{print hdr; print seq}' 1268_chr02.fasta > 1268_chr02.oneline.fasta
grep -c '^>' 1268_chr02.oneline.fasta
wc -l 1268_chr02.oneline.fasta
grep -v "^>" 1268_chr02.oneline.fasta | tr -cd 'N' | wc -c
#grep -v "^>" 1268_chr02.oneline.fasta | grep -o 'N\{1,\}' | awk '{print length($0)}' | sort -nr | head

#sed 's/N\{388\}/GCATGGGGACTGGCCTTTCCTGAGGCCACCAGAGCGGGTCCCTGAGGTCCCCGTCGTAAGTCGAGAGCACCTGCAGCATACTCGACTAGTGA/' 1268_chr02.oneline.fasta > 1268_chr02.gap_filled.oneline.fasta
#grep -v "^>" 1268_chr02.gap_filled.oneline.fasta | tr -cd 'N' | wc -c

# convert back to multiline 
awk 'NR==1{print; next} {for(i=1;i<=length($0);i+=60) print substr($0,i,60)}' 1268_chr02.gap_filled.oneline.fasta > 1268_chr02.gap_filled.fasta
