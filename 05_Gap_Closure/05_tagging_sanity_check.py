#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 05_tagging_sanity_check.py
#
# Sanity check on the tagged chromosome FASTAs: expected @@LEFT@@ / @@RIGHT@@ pairs are
# present, in order, and bracket a stretch of N's of the expected length.
# ------------------------------------------------------------------------------
# tagging sanity check 
    # ROI = reigon of interest [ genome before ] @@LEFT@@  <—— ROI ——>  @@RIGHT@@ [ genome after ]

from Bio import SeqIO

def read_single_fasta(path):
    rec = next(SeqIO.parse(path, "fasta"))
    return rec.id, str(rec.seq)

def roi_len_between_tags(seq, left_tag, right_tag):
    pL = seq.find(left_tag)
    pR = seq.find(right_tag)
    assert pL != -1, "LEFT tag not found"
    assert pR != -1, "RIGHT tag not found"
    assert pR > pL, "RIGHT tag occurs before LEFT tag"
    roi = seq[pL+len(left_tag):pR]
    return len(roi)

print("sample\tchr\tgap_len\texpected_ROI_len\tobserved_ROI_len\tOK")

for _, row in gaps.iterrows():
    sample    = row["sample"]
    chr_num   = row["chr"]
    gap_start = int(row["gap_start"])
    gap_end   = int(row["gap_end"])
    gap_len   = gap_end - gap_start + 1

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

    expected = (2*FLANK) + gap_len
    observed = roi_len_between_tags(seq, left_tag, right_tag)

    ok = (observed == expected)
    print(f"{sample}\t{chr_num}\t{gap_len}\t{expected}\t{observed}\t{ok}")
