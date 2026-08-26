#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 03_verify_gap_coordinates_on_fastas.py
#
# Coordinate-only verification on the original chromosome FASTAs (no modifications) - confirms
# every declared gap position is actually a run of N's of the expected length.
# ------------------------------------------------------------------------------
# coordinate verification only on the original FASTAs.
# Goal now: confirm for each gap that:
    # the chromosome ID is correct (including RagTag suffix)
    # the gap interval is an N-run
    # the gap boundaries match the true N-run (start/end)
    # optional: report the full N-run if your provided coords sit inside a longer N stretch

from pathlib import Path
from Bio import SeqIO
import pandas as pd

# ---- Explicit FASTA paths (authoritative) ----
FASTA_1268 = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta")
FASTA_1271 = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta")

def load_fasta_dict(fasta):
    return {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}

def resolve_chr_id(records, sample, chr_num):
    prefix = f"NLFDP{sample}_chr{chr_num}"
    matches = [k for k in records if k.startswith(prefix)]
    assert len(matches) == 1, f"[ERROR] Chromosome ID not unique for {prefix}: {matches}"
    return matches[0]

def is_all_N(s):
    s = str(s).upper()
    return len(s) > 0 and set(s) <= {"N"}

def expand_to_full_N_run(seq, start_1b, end_1b):
    """
    Given a candidate gap interval [start_1b, end_1b], expand left/right
    to cover the full contiguous N-run that contains it.
    Returns (full_start_1b, full_end_1b).
    """
    s = str(seq).upper()
    L = len(s)

    # convert to 0-based inclusive
    i = start_1b - 1
    j = end_1b - 1

    # if interval includes non-N, expansion isn't meaningful
    if not is_all_N(s[i:j+1]):
        return None

    # expand left
    a = i
    while a > 0 and s[a-1] == "N":
        a -= 1

    # expand right
    b = j
    while b < L-1 and s[b+1] == "N":
        b += 1

    return (a+1, b+1)  # back to 1-based inclusive

print("sample\tchr\tchr_id\tgiven_start\tgiven_end\tgiven_len\tgiven_all_N\tfullN_start\tfullN_end\tfullN_len")

for _, row in gaps.iterrows():
    sample = row["sample"]
    chr_num = row["chr"]
    start = int(row["gap_start"])
    end = int(row["gap_end"])

    fasta = FASTA_1268 if sample == "1268" else FASTA_1271
    records = load_fasta_dict(fasta)
    chr_id = resolve_chr_id(records, sample, chr_num)

    seq = records[chr_id].seq
    L = len(seq)

    assert 1 <= start < end <= L, f"[ERROR] Out-of-range coords: {sample} chr{chr_num} {start}-{end} (len={L})"

    gap_seq = seq[start-1:end]
    given_all_N = is_all_N(gap_seq)
    given_len = end - start + 1

    full = expand_to_full_N_run(seq, start, end)
    if full is None:
        full_start, full_end, full_len = "NA", "NA", "NA"
    else:
        full_start, full_end = full
        full_len = full_end - full_start + 1

    print(f"{sample}\t{chr_num}\t{chr_id}\t{start}\t{end}\t{given_len}\t{given_all_N}\t{full_start}\t{full_end}\t{full_len}")
