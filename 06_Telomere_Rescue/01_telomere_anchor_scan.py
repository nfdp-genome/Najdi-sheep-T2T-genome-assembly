#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 01_telomere_anchor_scan.py
#
# Scan every chromosome end for TTAGGG/CCCTAA in 1 kb windows across the terminal 50 kb;
# record first inward window with >=10 motifs; define 20 kb internal anchor per end.
# ------------------------------------------------------------------------------
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# -------------------------
# Parameters
# -------------------------
WINDOW = 1_000
SCAN_LEN = 50_000
ANCHOR_LEN = 20_000
MOTIF_THRESHOLD = 10

MOTIFS = ["TTAGGG", "CCCTAA"]


# -------------------------
# Input and output paths
# -------------------------
INPUT_FASTA = Path(
    "/ibex/project/c2293/najdi_t2t_project/"
    "WP1_genome_assembly_qc/data/processed/hifiasm_assembly/"
    "chromosomal_reconstruction/anchor_based_gap_filling/"
    "mapping_to_mergedAssemblies/1268/"
    "1268_NLFDP1268_chr04_5p_TELROI/"
    "1268_NLFDP1268_chr04_5p_TELROI.besthit.fasta"
)

OUT_TSV = Path(
    "/ibex/project/c2293/najdi_t2t_project/"
    "WP1_genome_assembly_qc/data/processed/hifiasm_assembly/"
    "chromosomal_reconstruction/anchor_based_gap_filling/"
    "telomere_anchor_scan_Test.tsv"
)


# -------------------------
# Helper functions
# -------------------------
def count_motifs(seq, motifs):
    """Count occurrences of all specified telomeric motifs."""
    sequence = seq.upper()
    return sum(sequence.count(motif) for motif in motifs)


def scan_end(seq, which_end):
    """
    Scan a terminal sequence region in 1 kb windows,
    moving from the chromosome terminus inward.
    """
    sequence_length = len(seq)
    scan_length = min(SCAN_LEN, sequence_length)
    rows = []

    if which_end == "5p":
        region = seq[:scan_length]
        offset = 0
        positions = range(0, len(region), WINDOW)

    elif which_end == "3p":
        region = seq[-scan_length:]
        offset = sequence_length - scan_length
        positions = range(len(region), 0, -WINDOW)

    else:
        raise ValueError("which_end must be either '5p' or '3p'")

    for window_index, position in enumerate(positions):
        if which_end == "5p":
            window_start = position
            window_end = min(position + WINDOW, len(region))
        else:
            window_end = position
            window_start = max(0, position - WINDOW)

        window_sequence = region[window_start:window_end]

        rows.append(
            {
                "window_index": window_index,
                "start_1b": offset + window_start + 1,
                "end_1b": offset + window_end,
                "motif_count": count_motifs(window_sequence, MOTIFS),
            }
        )

    return rows


def first_telomeric_window(windows):
    """Return the first terminal window meeting the motif threshold."""
    for window in windows:
        if window["motif_count"] >= MOTIF_THRESHOLD:
            return window

    return None


def define_anchor(sequence_length, which_end):
    """Define the 20 kb anchor on the internal side of the 50 kb scan region."""
    if sequence_length < SCAN_LEN:
        raise ValueError(
            f"Sequence length ({sequence_length:,} bp) is shorter than "
            f"SCAN_LEN ({SCAN_LEN:,} bp); anchor coordinates cannot be defined."
        )

    if which_end == "5p":
        anchor_start = SCAN_LEN - ANCHOR_LEN + 1
        anchor_end = SCAN_LEN

    elif which_end == "3p":
        anchor_start = sequence_length - SCAN_LEN + 1
        anchor_end = anchor_start + ANCHOR_LEN - 1

    else:
        raise ValueError("which_end must be either '5p' or '3p'")

    return anchor_start, anchor_end


def analyze_fasta(fasta_path, sample_label):
    """Analyse both ends of every sequence in a FASTA file."""
    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence = str(record.seq)
        chromosome_id = record.id
        sequence_length = len(sequence)

        for chromosome_end in ["5p", "3p"]:
            windows = scan_end(sequence, chromosome_end)
            telomeric_hit = first_telomeric_window(windows)

            telomeric_window_start = (
                telomeric_hit["start_1b"] if telomeric_hit else None
            )
            telomeric_motif_count = (
                telomeric_hit["motif_count"] if telomeric_hit else 0
            )

            anchor_start, anchor_end = define_anchor(
                sequence_length,
                chromosome_end,
            )

            records.append(
                {
                    "sample": sample_label,
                    "chromosome": chromosome_id,
                    "end": chromosome_end,
                    "chromosome_length": sequence_length,
                    "scan_length_bp": SCAN_LEN,
                    "window_size_bp": WINDOW,
                    "motif_threshold": MOTIF_THRESHOLD,
                    "first_telomeric_window_start_1b": telomeric_window_start,
                    "motifs_per_kb_at_start": telomeric_motif_count,
                    "anchor_start_1b": anchor_start,
                    "anchor_end_1b": anchor_end,
                }
            )

    return records


# -------------------------
# Run analysis
# -------------------------
rows = analyze_fasta(INPUT_FASTA, "1268")

df = pd.DataFrame(rows)

# Deterministic column order for downstream processing
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
# Write output
# -------------------------
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_TSV, sep="\t", index=False)

print(f"TSV written to: {OUT_TSV}")
print(f"Rows: {len(df)}")

df.head()
