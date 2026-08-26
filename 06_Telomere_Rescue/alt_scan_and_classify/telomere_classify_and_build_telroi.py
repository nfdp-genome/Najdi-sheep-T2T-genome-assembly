#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# telomere_classify_and_build_telroi.py
#
# Standalone A/B/C classification and TELROI decision table
# (functionally equivalent to 02_build_telroi_decision_table.py).
# ------------------------------------------------------------------------------
# STEP 1 (fresh): Build TELROI decision table ONLY
# - Reads telomere_anchor_scan.tsv
# - Selects Class B & C
# - Computes ROI length and coordinates
# - Generates TELROI tags
# - Writes telomere_roi_table.tsv
# - NO sequence operations

import pandas as pd
from pathlib import Path

# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

IN_TSV  = BASE_DIR / "telomere_anchor_scan.tsv"
OUT_TSV = BASE_DIR / "telomere_roi_table.tsv"

# ----------------------------
# Load input
# ----------------------------
df = pd.read_csv(IN_TSV, sep="\t")

# ----------------------------
# Telomere class definition
# ----------------------------
def classify_telomere(m):
    if m >= 150:
        return "A"
    elif 10 <= m < 150:
        return "B"
    else:
        return "C"

df["tel_class"] = df["motifs_per_kb_at_start"].apply(classify_telomere)

# Keep only Class B & C
df = df[df["tel_class"].isin(["B", "C"])].copy()

# ----------------------------
# ROI length policy (Option C: hybrid)
# ----------------------------
def choose_roi_length(row):
    if row["tel_class"] == "B":
        L = 50_000
    else:  # Class C
        L = 200_000
    # caps
    return min(max(L, 20_000), 200_000)

df["roi_len"] = df.apply(choose_roi_length, axis=1)

# ----------------------------
# ROI coordinates (1-based)
# ----------------------------
def compute_roi_coords(row):
    chr_len = int(row["chromosome_length"])
    L = int(row["roi_len"])

    if row["end"] == "5p":
        start = 1
        end = min(chr_len, L)
    else:  # 3p
        start = max(1, chr_len - L + 1)
        end = chr_len

    return pd.Series({
        "roi_start_1b": start,
        "roi_end_1b": end
    })

df[["roi_start_1b", "roi_end_1b"]] = df.apply(compute_roi_coords, axis=1)

# ----------------------------
# Anchor length
# ----------------------------
df["anchor_len"] = df["roi_end_1b"] - df["roi_start_1b"] + 1

# ----------------------------
# TELROI tag
# ----------------------------
def make_telroi_tag(row):
    return (
        f"TELROI::{row['sample']}::{row['chromosome']}::"
        f"{row['end']}::{row['tel_class']}::{row['anchor_len']}"
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
    "anchor_len",
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
