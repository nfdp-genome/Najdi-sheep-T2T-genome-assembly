#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 02_build_telroi_decision_table.py
#
# Classify each end (Class A >=150, Class B 10-149, Class C <10 motifs/kb); keep only B and C;
# compute 1-based ROI coordinates (50 kb for B, 200 kb for C), actual ROI length, and a unique
# TELROI_tag. No sequence is modified.
# ------------------------------------------------------------------------------
# STEP 1: Build the TELROI decision table only
# - Reads telomere_anchor_scan.tsv
# - Classifies chromosome ends as A, B, or C
# - Selects Class B and C ends
# - Computes ROI lengths and coordinates
# - Generates TELROI tags
# - Writes telomere_roi_table.tsv
# - Performs no sequence operations

from pathlib import Path

import pandas as pd


# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

IN_TSV = BASE_DIR / "telomere_anchor_scan.tsv"
OUT_TSV = BASE_DIR / "telomere_roi_table.tsv"


# ----------------------------
# Load input
# ----------------------------
df = pd.read_csv(IN_TSV, sep="\t")


# ----------------------------
# Telomere class definition
# ----------------------------
def classify_telomere(motif_count):
    if motif_count >= 150:
        return "A"
    elif motif_count >= 10:
        return "B"
    else:
        return "C"


df["tel_class"] = df["motifs_per_kb_at_start"].apply(classify_telomere)

# Keep only Class B and C chromosome ends
df = df[df["tel_class"].isin(["B", "C"])].copy()


# ----------------------------
# Requested ROI length
# ----------------------------
def choose_roi_length(row):
    if row["tel_class"] == "B":
        roi_length = 50_000
    else:  # Class C
        roi_length = 200_000

    return min(max(roi_length, 20_000), 200_000)


df["roi_len"] = df.apply(choose_roi_length, axis=1)


# ----------------------------
# ROI coordinates (1-based)
# ----------------------------
def compute_roi_coords(row):
    chromosome_length = int(row["chromosome_length"])
    roi_length = int(row["roi_len"])

    if row["end"] == "5p":
        roi_start = 1
        roi_end = min(chromosome_length, roi_length)

    elif row["end"] == "3p":
        roi_start = max(1, chromosome_length - roi_length + 1)
        roi_end = chromosome_length

    else:
        raise ValueError(
            f"Unexpected chromosome end: {row['end']}. "
            "Expected '5p' or '3p'."
        )

    return pd.Series(
        {
            "roi_start_1b": roi_start,
            "roi_end_1b": roi_end,
        }
    )


df[["roi_start_1b", "roi_end_1b"]] = df.apply(
    compute_roi_coords,
    axis=1,
)


# ----------------------------
# Actual ROI length
# ----------------------------
df["actual_roi_len"] = (
    df["roi_end_1b"] - df["roi_start_1b"] + 1
)


# ----------------------------
# TELROI tag
# ----------------------------
def make_telroi_tag(row):
    return (
        f"TELROI::{row['sample']}::{row['chromosome']}::"
        f"{row['end']}::{row['tel_class']}::{row['actual_roi_len']}"
    )


df["TELROI_tag"] = df.apply(make_telroi_tag, axis=1)


# ----------------------------
# Final table
# ----------------------------
final_cols = [
    "sample",
    "chromosome",
    "end",
    "tel_class",
    "motifs_per_kb_at_start",
    "chromosome_length",
    "roi_start_1b",
    "roi_end_1b",
    "roi_len",
    "actual_roi_len",
    "TELROI_tag",
]

df = (
    df[final_cols]
    .sort_values(["sample", "chromosome", "end"])
    .reset_index(drop=True)
)


# ----------------------------
# Write output
# ----------------------------
df.to_csv(OUT_TSV, sep="\t", index=False)

print(f"[OK] TELROI decision table written to:\n{OUT_TSV}")
print(f"[INFO] Total TELROIs (Class B + C): {len(df)}")

df.head()
