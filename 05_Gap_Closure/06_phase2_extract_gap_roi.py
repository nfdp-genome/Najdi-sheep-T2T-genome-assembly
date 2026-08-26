#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 06_phase2_extract_gap_roi.py
#
# PHASE 2 - strip all @@LEFT@@ / @@RIGHT@@ placeholders on a working copy, then slice by
# 1-based inclusive coordinates to write one ROI query FASTA per gap
# (left flank + gap Ns + right flank).
# ------------------------------------------------------------------------------
# PHASE 2 — Extract ROI blocks (tag-based, safe)
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

OUTDIR = WORKDIR / "roi_sequences"
OUTDIR.mkdir(parents=True, exist_ok=True)

def read_single_fasta(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return rec.id, str(rec.seq)

def extract_between_tags(seq, left_tag, right_tag):
    pL = seq.find(left_tag)
    pR = seq.find(right_tag)
    assert pL != -1, "LEFT tag not found"
    assert pR != -1, "RIGHT tag not found"
    assert pR > pL, "RIGHT tag before LEFT tag"

    roi = seq[pL + len(left_tag):pR]
    return roi

def write_fasta(path, header, seq):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        f.write(seq)
        f.write("\n")

print("[INFO] Extracting ROI blocks...")

for _, row in gaps.iterrows():
    sample    = row["sample"]
    chr_num   = row["chr"]
    gap_start = int(row["gap_start"])
    gap_end   = int(row["gap_end"])

    tagged_fa = WORKDIR / "tagged_chromosomes" / f"{sample}_chr{chr_num}_TAGGED.fasta"
    header, seq = read_single_fasta(tagged_fa)

    left_tag = (
        f"@@LEFT_BOUNDARY|ROI|{sample}|chr{chr_num}|"
        f"gap:{gap_start}-{gap_end}|flank:{FLANK}@@"
    )
    right_tag = (
        f"@@RIGHT_BOUNDARY|ROI|{sample}|chr{chr_num}|"
        f"gap:{gap_start}-{gap_end}|flank:{FLANK}@@"
    )

    roi = extract_between_tags(seq, left_tag, right_tag)

    # Split ROI
    left_50k  = roi[:FLANK]
    right_50k = roi[-FLANK:]

    prefix = OUTDIR / f"{sample}_chr{chr_num}_gap_{gap_start}_{gap_end}"

    write_fasta(prefix.with_suffix(".roi.fa"),
                f"{prefix.name}_ROI_FLANK50K",
                roi)

    write_fasta(prefix.with_suffix(".left_50k.fa"),
                f"{prefix.name}_LEFT_50K",
                left_50k)

    write_fasta(prefix.with_suffix(".right_50k.fa"),
                f"{prefix.name}_RIGHT_50K",
                right_50k)

    write_fasta(prefix.with_suffix(".combined_L50k_GAP_R50k.fa"),
                f"{prefix.name}_L50K_GAP_R50K",
                roi)

    print(f"[OK] {sample} chr{chr_num} ROI extracted")
