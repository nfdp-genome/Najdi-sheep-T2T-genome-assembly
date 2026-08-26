#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# telomere_scan_chromosomes.py
#
# Standalone chromosome-FASTA scan for TTAGGG/CCCTAA repeats
# (functionally equivalent to 01_telomere_anchor_scan.py).
# ------------------------------------------------------------------------------
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# -------------------------
# Parameters (locked)
# -------------------------
WINDOW = 1_000
SCAN_LEN = 50_000
ANCHOR_LEN = 20_000
MOTIF_THRESHOLD = 10

MOTIFS = ["TTAGGG", "CCCTAA"]  # reverse complements are symmetric here
MOTIFS = ["TTAGGG", "CCCTAA"]  # reverse complements are symmetric here

FASTA_1268 = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta")
FASTA_1271 = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/unplaced_contigs_vs_scaffolded_assembly/1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta")

#FASTA_1268 = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/1268/1268_NLFDP1268_chr04_5p_TELROI/1268_NLFDP1268_chr04_5p_TELROI.besthit.fasta")

OUT_TSV = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/telomere_anchor_scan_Test.tsv")

# -------------------------
# Helper functions
# -------------------------
def count_motifs(seq, motifs):
    s = seq.upper()
    return sum(s.count(m) for m in motifs)

def scan_end(seq, which_end):
    L = len(seq)
    rows = []

    if which_end == "5p":
        region = seq[:SCAN_LEN]
        offset = 0
    else:
        region = seq[-SCAN_LEN:]
        offset = L - SCAN_LEN

    for i in range(0, SCAN_LEN, WINDOW):
        win = region[i:i + WINDOW]
        rows.append({
            "window_index": i // WINDOW,
            "start_1b": offset + i + 1,
            "end_1b": offset + i + len(win),
            "motif_count": count_motifs(win, MOTIFS)
        })

    return rows

def first_telomeric_window(windows):
    for w in windows:
        if w["motif_count"] >= MOTIF_THRESHOLD:
            return w
    return None

def analyze_fasta(fasta, sample_label):
    records = []

    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq)
        chr_id = rec.id
        L = len(seq)

        for end in ["5p", "3p"]:
            windows = scan_end(seq, end)
            hit = first_telomeric_window(windows)

            tel_start = hit["start_1b"] if hit else None
            tel_motifs = hit["motif_count"] if hit else 0

            # Anchor definition (internal side of scan region)
            if end == "3p":
                anchor_end = L - SCAN_LEN + ANCHOR_LEN
                anchor_start = anchor_end - ANCHOR_LEN + 1
            else:
                anchor_start = SCAN_LEN - ANCHOR_LEN + 1
                anchor_end = SCAN_LEN

            records.append({
                "sample": sample_label,
                "chromosome": chr_id,
                "end": end,
                "chromosome_length": L,
                "scan_length_bp": SCAN_LEN,
                "window_size_bp": WINDOW,
                "motif_threshold": MOTIF_THRESHOLD,
                "first_telomeric_window_start_1b": tel_start,
                "motifs_per_kb_at_start": tel_motifs,
                "anchor_start_1b": anchor_start,
                "anchor_end_1b": anchor_end
            })

    return records

# -------------------------
# Run analysis
# -------------------------
rows = []
rows += analyze_fasta(FASTA_1268, "1268")
rows += analyze_fasta(FASTA_1271, "1271")

df = pd.DataFrame(rows)

# deterministic column order (important for downstream pipelines)
df = df[
    [
        "sample",
        "chromosome",
        "end",
        "chromosome_length",
        "scan_length_bp",
        "window_size_bp",
        "motif_threshold",
        "first_telomeric_window_start_1b",
        "motifs_per_kb_at_start",
        "anchor_start_1b",
        "anchor_end_1b",
    ]
]

# -------------------------
# Write TSV
# -------------------------
df.to_csv(OUT_TSV, sep="\t", index=False)

print(f"TSV written to: {OUT_TSV}")
print(f"Rows: {len(df)}")

df.head()
