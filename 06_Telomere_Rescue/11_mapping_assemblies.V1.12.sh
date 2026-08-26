#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# 11_mapping_assemblies.V1.12.sh
#
# SLURM array job invoked by script 10. Runs minimap2 -x asm5 for each TELROI query FASTA
# against the sample's merged assembly pool, producing a per-ROI PAF file (with optional sorted
# BAM for IGV inspection). Best-hit selection logic for GAP-type ROIs is embedded (besthit_gap
# awk block); non-GAP ROIs produce PAF only, to be processed by scripts 12-15.
# ------------------------------------------------------------------------------
#SBATCH --job-name=bam_map_ROI_to_mergedAssemblies
#SBATCH --output=/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/logs/bam_map_ROI_Unique_%A_%a.out
#SBATCH --error=/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/logs/bam_map_ROI_Unique_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --constraint=intel
#SBATCH --array=0-98%30

set -euo pipefail

# ------------------------------------------------------------
# Mode control: all | nobam | bam
# ------------------------------------------------------------
MODE="${1:-all}"

THREADS="${SLURM_CPUS_PER_TASK:-40}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"

module purge
module load minimap2/2.24
module load samtools/1.16.1

WORK="/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
RODIR="$WORK/roi_sequences"
LOGDIR="$WORK/logs"
mkdir -p "$LOGDIR"

REF_1268="$WORK/1268_all_assemblies.merged.unique.fasta"
REF_1271="$WORK/1271_all_assemblies.merged.unique.fasta"

LIST_1268="$RODIR/1268_roi_sequences_paths.txt"
LIST_1271="$RODIR/1271_roi_sequences_paths.txt"

[[ -s "$LIST_1268" ]] || { echo "[ERROR] Missing list: $LIST_1268"; exit 1; }
[[ -s "$LIST_1271" ]] || { echo "[ERROR] Missing list: $LIST_1271"; exit 1; }
[[ -s "$REF_1268" ]] || { echo "[ERROR] Missing reference: $REF_1268"; exit 1; }
[[ -s "$REF_1271" ]] || { echo "[ERROR] Missing reference: $REF_1271"; exit 1; }

# ------------------------------------------------------------
# Read ROI lists
# ------------------------------------------------------------
mapfile -t L1268 < <(grep -vE '^\s*($|#)' "$LIST_1268")
mapfile -t L1271 < <(grep -vE '^\s*($|#)' "$LIST_1271")

N1268="${#L1268[@]}"
N1271="${#L1271[@]}"
N=$((N1268 + N1271))

[[ "$N" -gt 0 ]] || { echo "[ERROR] No ROI paths found"; exit 1; }
[[ "$IDX" -lt "$N" ]] || { echo "[ERROR] Array index out of range"; exit 1; }

if [[ "$IDX" -lt "$N1268" ]]; then
  SAMPLE="1268"
  UNP="${L1268[$IDX]}"
  CHR="$REF_1268"
else
  SAMPLE="1271"
  UNP="${L1271[$((IDX - N1268))]}"
  CHR="$REF_1271"
fi

[[ -s "$UNP" ]] || { echo "[ERROR] Missing query fasta: $UNP"; exit 1; }

BASE="$(basename "$UNP")"
RID="${BASE%.fa}"
RID="${RID%.fasta}"

# merged refs are already indexed; this is a harmless check (fast if present)
samtools faidx "$CHR" || true
# query indexing is optional; harmless check
samtools faidx "$UNP" || true

OUTDIR="$WORK/mapping_to_mergedAssemblies_withBAM/$SAMPLE/$RID"
mkdir -p "$OUTDIR"

PAF="$OUTDIR/${RID}.vs_mergedAssemblies.paf"
BEST="$OUTDIR/${RID}.besthit.tsv"
BAM="$OUTDIR/${RID}.vs_mergedAssemblies.bam"

echo "======================================================"
echo "[ROI2mergedAssemblies] $(date)"
echo "  Mode     : $MODE"
echo "  Task     : $IDX / $((N-1))"
echo "  Sample   : $SAMPLE"
echo "  ROI      : $RID"
echo "======================================================"

# ------------------------------------------------------------
# Best-hit selectors
# ------------------------------------------------------------
besthit_gap() {
  awk -v LEFT_MAX=50000 -v RIGHT_ANCHOR_END=54000 'BEGIN{
    OFS="\t";
    print "roi_query","query_length_bp","query_start_0based","query_end_0based","strand",
          "target_contig_id","target_length_bp","target_start_0based","target_end_0based",
          "matching_bases","alignment_length_bp","mapping_quality","alignment_type",
          "minimizer_matches","chaining_score_best","chaining_score_second_best",
          "sequence_divergence","repetitive_seed_length_bp","identity","qcov";
  }
  $13=="tp:A:P" && $12>=20 && $4<=LEFT_MAX && $5>=RIGHT_ANCHOR_END && $11>=40000 {
    q=$1; qlen=$2;              # <-- FIX: define q (and qlen here)
    aln=$11; matches=$10;

    idn=(aln>0?matches/aln:0);
    qcov=(qlen>0?aln/qlen:0);

    if (!(q in best) || aln>best[q]) {
      best[q]=aln;
      line[q]=$0 OFS idn OFS qcov
    }
  }
  END{for (q in line) print line[q]}'
}

# ------------------------------------------------------------
# STEP 1 + STEP 2
# ------------------------------------------------------------
if [[ "$MODE" != "bam" ]]; then
  minimap2 -x asm5 -t "$THREADS" "$CHR" "$UNP" > "$PAF"

  if [[ "$RID" == *GAP* ]]; then
    besthit_gap < "$PAF" > "$BEST"
  else
    : > "$BEST"
  fi
fi

# ------------------------------------------------------------
# Extract best-hit FASTA
# ------------------------------------------------------------
extract_besthit_fasta() {
  [[ -s "$BEST" && $(wc -l < "$BEST") -ge 2 ]] || return 0
  CONTIG=$(awk 'NR==2{print $6}' "$BEST")
  samtools faidx "$CHR" "$CONTIG" > "$OUTDIR/${RID}.besthit.fasta"
}

# ------------------------------------------------------------
# BAM mapping (NO PIPES)
# ------------------------------------------------------------
run_bam_mapping() {

  local SAM="$OUTDIR/${RID}.vs_mergedAssemblies.sam"

  # 1) write SAM (no pipe)
  minimap2 --split-prefix "$OUTDIR/${RID}.mm2split" \
    -a -x asm5 -t "$THREADS" \
    "$CHR" "$UNP" > "$SAM"

  # optional hard fail early if header is broken
  # samtools view -H "$SAM" >/dev/null

  # 2) sort SAM -> BAM
  samtools sort -@ "$THREADS" -T "$OUTDIR/${RID}.tmp" -o "$BAM" "$SAM"

  # 3) index BAM
  samtools index -@ "$THREADS" "$BAM"

  # optional cleanup
  # rm -f "$SAM"
}

# ------------------------------------------------------------
# Mode dispatch
# ------------------------------------------------------------
case "$MODE" in
  all)   extract_besthit_fasta; run_bam_mapping ;;
  nobam) extract_besthit_fasta ;;
  bam)   run_bam_mapping ;;
esac

echo "[DONE]"
