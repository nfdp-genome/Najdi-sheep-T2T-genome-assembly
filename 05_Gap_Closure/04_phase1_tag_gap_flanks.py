#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 04_phase1_tag_gap_flanks.py
#
# PHASE 1A - insert @@LEFT@@ / @@RIGHT@@ boundary placeholders on a working copy of each
# chromosome FASTA so the gap ROI (flank + gap + flank) can be sliced unambiguously in the
# next step.
# ------------------------------------------------------------------------------
# PHASE 1A — TAG GAP REGIONS IN CHROMOSOME FASTAs (NO EXTRACTION YET) — TAG INSERTION (NO REPLACEMENT)
# What this does (precisely)
# Keeps the chromosome sequence unchanged
# The tags are inserted at the edges of the 50 kb flanks, so that everything between the tags is the sequence of interest that will later be replaced.

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

FLANK = 50_000

TAGGED_DIR = WORKDIR / "tagged_chromosomes"
TAGGED_DIR.mkdir(parents=True, exist_ok=True)

def load_fasta_dict(fasta):
    return {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}

def resolve_chr_id(records, sample, chr_num):
    prefix = f"NLFDP{sample}_chr{chr_num}"
    matches = [k for k in records if k.startswith(prefix)]
    assert len(matches) == 1, f"[ERROR] Chromosome ID not unique for {prefix}: {matches}"
    return matches[0]

def insert_tags_enclosing_roi(seq, gap_start, gap_end, flank, tag_left, tag_right):
    """
    Tags enclose ROI = [gap_start-flank, gap_end+flank] (1-based inclusive).
    Insert LEFT tag before ROI, RIGHT tag after ROI. ROI sequence unchanged.
    """
    L = len(seq)
    leftB_0  = max(0, (gap_start - flank) - 1)  # 0-based start
    rightB_0 = min(L, gap_end + flank)          # 0-based end (python slice end)

    return seq[:leftB_0] + tag_left + seq[leftB_0:rightB_0] + tag_right + seq[rightB_0:]


for _, row in gaps.iterrows():
    sample    = row["sample"]
    chr_num   = row["chr"]
    gap_start = int(row["gap_start"])
    gap_end   = int(row["gap_end"])

    fasta = FASTA_1268 if sample == "1268" else FASTA_1271
    records = load_fasta_dict(fasta)

    chr_id = resolve_chr_id(records, sample, chr_num)
    rec = records[chr_id]

    tag_left = (
        f"@@LEFT_BOUNDARY|ROI|{sample}|chr{chr_num}|"
        f"gap:{gap_start}-{gap_end}|flank:{FLANK}@@"
    )
    tag_right = (
        f"@@RIGHT_BOUNDARY|ROI|{sample}|chr{chr_num}|"
        f"gap:{gap_start}-{gap_end}|flank:{FLANK}@@"
    )

    tagged_seq = insert_tags_enclosing_roi(
        str(rec.seq),
        gap_start,
        gap_end,
        FLANK,
        tag_left,
        tag_right
    )

    out_fa = TAGGED_DIR / f"{sample}_chr{chr_num}_TAGGED.fasta"
    SeqIO.write(
        SeqRecord(tagged_seq, id=rec.id, description="ROI_TAGGED_FLANK50K"),
        out_fa,
        "fasta"
    )

    print(f"[TAGGED_OK] {sample} chr{chr_num} -> {out_fa}")
