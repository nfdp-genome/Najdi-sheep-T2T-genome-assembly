#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 07_extract_telroi_query_sequences.py
#
# Strip all placeholders on a working copy, then slice by 1-based inclusive TELROI
# coordinates and write one query FASTA per end.
# ------------------------------------------------------------------------------
# STEP 4: Extract TELROI query sequences for mapping
# - Reads telomere_roi_table.tsv
# - Reads tagged chromosome FASTAs
# - Removes all TELROI placeholders
# - Extracts each ROI using 1-based inclusive coordinates
# - Writes one FASTA per TELROI

from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO


# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

ROI_TABLE = BASE_DIR / "telomere_roi_table.tsv"
TAGGED_DIR = BASE_DIR / "tagged_chromosomes"
OUT_DIR = BASE_DIR / "roi_sequences"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------
# Load TELROI table
# ----------------------------
df = pd.read_csv(
    ROI_TABLE,
    sep="\t",
    dtype={"sample": str},
)

required_columns = {
    "sample",
    "chromosome",
    "end",
    "tel_class",
    "actual_roi_len",
    "roi_start_1b",
    "roi_end_1b",
}

missing_columns = required_columns.difference(df.columns)

if missing_columns:
    raise ValueError(
        "TELROI table is missing required columns: "
        + ", ".join(sorted(missing_columns))
    )


# ----------------------------
# TELROI placeholder pattern
# ----------------------------
TELROI_PATTERN = re.compile(
    r"@@TELROI\|[^@]*@@"
)


def placeholder_from_row(row):
    """Build the expected TELROI FASTA header."""
    return (
        f"@@TELROI|{row['sample']}|{row['chromosome']}|"
        f"{row['end']}|{row['tel_class']}|"
        f"{int(row['actual_roi_len'])}|"
        f"roi:{int(row['roi_start_1b'])}-"
        f"{int(row['roi_end_1b'])}@@"
    )


def read_single_fasta_sequence(fasta_path):
    """Read a FASTA file containing exactly one sequence record."""
    records = list(SeqIO.parse(str(fasta_path), "fasta"))

    if len(records) != 1:
        raise ValueError(
            f"Expected exactly one FASTA record in {fasta_path}, "
            f"but found {len(records)}."
        )

    return records[0].id, str(records[0].seq)


# ----------------------------
# Extract TELROI sequences
# ----------------------------
written = 0

for _, row in df.iterrows():
    sample = str(row["sample"])
    chromosome = str(row["chromosome"])
    chromosome_end = row["end"]

    roi_start = int(row["roi_start_1b"])
    roi_end = int(row["roi_end_1b"])
    expected_roi_length = int(row["actual_roi_len"])

    tag = placeholder_from_row(row)

    tagged_fasta = (
        TAGGED_DIR
        / f"{sample}_{chromosome}_TELROI_TAGGED.fasta"
    )

    if not tagged_fasta.exists():
        raise FileNotFoundError(
            f"Missing tagged chromosome FASTA: {tagged_fasta}"
        )

    record_id, tagged_sequence = read_single_fasta_sequence(
        tagged_fasta
    )

    if record_id != chromosome:
        raise ValueError(
            f"FASTA header mismatch in {tagged_fasta.name}: "
            f"expected '{chromosome}', found '{record_id}'."
        )

    if tag not in tagged_sequence:
        raise ValueError(
            f"Expected TELROI tag was not found in "
            f"{tagged_fasta.name}:\n{tag}"
        )

    # Remove every TELROI placeholder to recover the original chromosome
    original_sequence = TELROI_PATTERN.sub("", tagged_sequence)

    chromosome_length = len(original_sequence)

    if (
        roi_start < 1
        or roi_end > chromosome_length
        or roi_start > roi_end
    ):
        raise ValueError(
            f"ROI coordinates are out of bounds for "
            f"{sample} {chromosome} {chromosome_end}: "
            f"{roi_start}-{roi_end} "
            f"(chromosome length: {chromosome_length})"
        )

    # Convert 1-based inclusive coordinates to Python slicing
    roi_sequence = original_sequence[roi_start - 1:roi_end]

    if len(roi_sequence) != expected_roi_length:
        raise ValueError(
            f"Extracted ROI length mismatch for "
            f"{sample} {chromosome} {chromosome_end}: "
            f"expected {expected_roi_length}, "
            f"observed {len(roi_sequence)}."
        )

    output_fasta = (
        OUT_DIR
        / f"{sample}_{chromosome}_{chromosome_end}_TELROI.fa"
    )

    with output_fasta.open("w") as handle:
        handle.write(f">{tag}\n")

        for position in range(0, len(roi_sequence), 80):
            handle.write(
                roi_sequence[position:position + 80] + "\n"
            )

    written += 1

    print(
        f"[EXTRACTED_OK] {sample} {chromosome} "
        f"{chromosome_end} -> {output_fasta}"
    )


print(f"[OK] Wrote {written} TELROI FASTA files to: {OUT_DIR}")

if not df.empty:
    example_path = (
        OUT_DIR
        / f"{df.iloc[0]['sample']}_"
        f"{df.iloc[0]['chromosome']}_"
        f"{df.iloc[0]['end']}_TELROI.fa"
    )

    print(f"[OK] Example: {example_path}")
