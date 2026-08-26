#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 08_patch_1271_chr08_and_chr17.sh
#
# Patch selected gap-filling contigs into chromosomes 08 and 17 of sample 1271
# (the two multi-contig chromosomes that were joined earlier by RagTag).
# ------------------------------------------------------------------------------
#patch chromosome 08 and 17 to sample 1271:
set -euo pipefail
module load samtools

# -----------------------
# Inputs
# -----------------------
GEN="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta"

CHR08="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq/1271_chr08_scaffolded_HuF.gap_filled.fasta"
CHR17="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/gap_filled_seq/1271_chr17_scaffolded_HuF.no_gap.fasta"

OUT="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.CHR08_CHR17_REPLACED.fasta"

# -----------------------
# Ensure indexed + detect exact IDs
# -----------------------
[[ -f "${GEN}.fai" ]] || samtools faidx "$GEN"

ID08=$(cut -f1 "${GEN}.fai" | grep -E '^NLFDP1271_chr08($|_)' | head -n 1 || true)
ID17=$(cut -f1 "${GEN}.fai" | grep -E '^NLFDP1271_chr17($|_)' | head -n 1 || true)

echo "Picked ID08: ${ID08:-NONE}"
echo "Picked ID17: ${ID17:-NONE}"
[[ -n "${ID08:-}" ]] || { echo "[ERROR] chr08 ID not found in genome"; exit 1; }
[[ -n "${ID17:-}" ]] || { echo "[ERROR] chr17 ID not found in genome"; exit 1; }

# -----------------------
# Build new genome: all contigs except chr08/chr17
# -----------------------
samtools faidx "$GEN" $(cut -f1 "${GEN}.fai" | grep -vFx "$ID08" | grep -vFx "$ID17") > "$OUT"

# Append new chr08 + chr17, forcing headers to match genome IDs
awk -v id="$ID08" 'BEGIN{done=0} /^>/{if(!done){print ">"id; done=1; next}} {print}' "$CHR08" >> "$OUT"
awk -v id="$ID17" 'BEGIN{done=0} /^>/{if(!done){print ">"id; done=1; next}} {print}' "$CHR17" >> "$OUT"

# Index output
samtools faidx "$OUT"

echo "[DONE] Wrote: $OUT"
echo
echo "Now validating md5..."
echo "chr08:"
samtools faidx "$OUT" "$ID08" | grep -v '^>' | tr -d '\n' | md5sum
grep -v '^>' "$CHR08" | tr -d '\n' | md5sum
echo
echo "chr17:"
samtools faidx "$OUT" "$ID17" | grep -v '^>' | tr -d '\n' | md5sum
grep -v '^>' "$CHR17" | tr -d '\n' | md5sum
