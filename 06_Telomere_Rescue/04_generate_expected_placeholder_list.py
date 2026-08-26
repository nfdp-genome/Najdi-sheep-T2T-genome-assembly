#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 04_generate_expected_placeholder_list.py
#
# QC helper: enumerate the placeholders that must appear given the decision table.
# ------------------------------------------------------------------------------
# Quality check: generate the expected TELROI placeholders

from pathlib import Path

import pandas as pd


# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

ROI_TABLE = BASE_DIR / "telomere_roi_table.tsv"
OUT_TXT = BASE_DIR / "expected_placeholders.txt"


# ----------------------------
# Load TELROI table
# ----------------------------
df = pd.read_csv(ROI_TABLE, sep="\t")


# ----------------------------
# Build expected placeholders
# ----------------------------
df["expected_placeholder"] = df.apply(
    lambda row: (
        f"@@TELROI|{row['sample']}|{row['chromosome']}|"
        f"{row['end']}|{row['tel_class']}|"
        f"{int(row['actual_roi_len'])}|"
        f"roi:{int(row['roi_start_1b'])}-{int(row['roi_end_1b'])}@@"
    ),
    axis=1,
)


# ----------------------------
# Write output
# ----------------------------
df["expected_placeholder"].to_csv(
    OUT_TXT,
    index=False,
    header=False,
)

print(f"[OK] Expected placeholders written to:\n{OUT_TXT}")
print(f"[INFO] Total expected placeholders: {len(df)}")
