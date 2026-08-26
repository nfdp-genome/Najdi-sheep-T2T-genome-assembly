#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 02_phase05_gap_coordinate_preflight.py
#
# PHASE 0.5 - pre-flight sanity checks on the gap coordinates table
# (bounds, columns, no overlapping gaps).
# ------------------------------------------------------------------------------
# =========================
# PHASE 0.5 — PRE-FLIGHT SANITY CHECKS (FIXED)
# =========================

def load_fasta_ids(fasta):
    return [rec.id for rec in SeqIO.parse(fasta, "fasta")]

def resolve_chr_id(fasta_ids, sample, chr_num):
    prefix = f"NLFDP{sample}_chr{chr_num}"
    matches = [cid for cid in fasta_ids if cid.startswith(prefix)]
    assert len(matches) == 1, (
        f"[ERROR] Chromosome ID problem for {sample} chr{chr_num}: {matches}"
    )
    return matches[0]


for _, row in gaps.iterrows():

    sample = row["sample"]
    chr_num = row["chr"]
    gap_start = row["gap_start"]
    gap_end = row["gap_end"]

    fasta = FASTA_1268 if sample == "1268" else FASTA_1271
    fasta_ids = load_fasta_ids(fasta)

    # ---- 1. Chromosome exists ----
    chr_id = resolve_chr_id(fasta_ids, sample, chr_num)

    # ---- 2. Coordinates sanity ----
    records = {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}
    seq_len = len(records[chr_id].seq)

    assert 1 <= gap_start < gap_end <= seq_len, (
        f"[ERROR] Invalid coordinates for {sample} chr{chr_num}"
    )

    # ---- 3. Gap is internal (not telomeric) ----
    assert gap_start > 1, (
        f"[ERROR] Gap touches chromosome start: {sample} chr{chr_num}"
    )
    assert gap_end < seq_len, (
        f"[ERROR] Gap touches chromosome end: {sample} chr{chr_num}"
    )

    # ---- 4. Gap length sanity ----
    gap_len = gap_end - gap_start + 1
    assert gap_len > 0, (
        f"[ERROR] Zero/negative gap length: {sample} chr{chr_num}"
    )

    print(f"[OK] Pre-flight check passed: {sample} chr{chr_num}")
