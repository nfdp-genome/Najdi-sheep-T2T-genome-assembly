#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 06_validate_telroi_tag_placement.py
#
# QC pass 2: row-level sanity table checking each placeholder's position, coordinates,
# ROI length, chromosome-ID match, and chromosome bounds.
# ------------------------------------------------------------------------------
# STEP 3: Validate TELROI placeholder placement and metadata
# - Reads the TELROI decision table
# - Parses TELROI placeholders from tagged chromosome FASTAs
# - Checks tag position, coordinates, ROI length, chromosome ID and table matching
# - Writes detailed and summary sanity tables

from pathlib import Path
import re

import pandas as pd


# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

TAGGED_DIR = BASE_DIR / "tagged_chromosomes"
ROI_TABLE = BASE_DIR / "telomere_roi_table.tsv"

OUT_TSV = BASE_DIR / "telroi_tagging_sanity_table.tsv"
OUT_SUMMARY_TSV = BASE_DIR / "telroi_tagging_sanity_summary.tsv"


if not TAGGED_DIR.exists():
    raise FileNotFoundError(f"Missing tagged FASTA directory: {TAGGED_DIR}")

if not ROI_TABLE.exists():
    raise FileNotFoundError(f"Missing TELROI decision table: {ROI_TABLE}")


# ----------------------------
# Load expected TELROI table
# ----------------------------
roi = pd.read_csv(
    ROI_TABLE,
    sep="\t",
    dtype=str,
)

roi.columns = [column.strip() for column in roi.columns]

# Use a consistent internal chromosome-column name
if "chromosome" in roi.columns and "chrom" not in roi.columns:
    roi = roi.rename(columns={"chromosome": "chrom"})

required_columns = [
    "sample",
    "chrom",
    "end",
    "tel_class",
    "actual_roi_len",
    "roi_start_1b",
    "roi_end_1b",
    "chromosome_length",
]

missing_columns = [
    column
    for column in required_columns
    if column not in roi.columns
]

if missing_columns:
    raise ValueError(
        "TELROI decision table is missing columns: "
        + ", ".join(missing_columns)
    )

numeric_columns = [
    "actual_roi_len",
    "roi_start_1b",
    "roi_end_1b",
    "chromosome_length",
]

for column in numeric_columns:
    roi[column] = roi[column].astype(int)


# ----------------------------
# TELROI placeholder pattern
# ----------------------------
# Format:
# @@TELROI|sample|chrom|end|class|actual_roi_len|roi:start-end@@

TELROI_PATTERN = re.compile(
    r"@@TELROI\|"
    r"(?P<sample>[^|]+)\|"
    r"(?P<chrom>[^|]+)\|"
    r"(?P<end>5p|3p)\|"
    r"(?P<tel_class>[A-Z])\|"
    r"(?P<actual_roi_len>\d+)\|"
    r"roi:"
    r"(?P<roi_start_1b>\d+)-"
    r"(?P<roi_end_1b>\d+)"
    r"@@"
)


# ----------------------------
# Read a single-record FASTA
# ----------------------------
def read_single_record_fasta(fasta_path):
    with fasta_path.open("r") as handle:
        header = handle.readline().strip()

        if not header.startswith(">"):
            raise ValueError(
                f"Invalid FASTA header in: {fasta_path}"
            )

        header_id = header[1:].split()[0]

        sequence = "".join(
            line.strip()
            for line in handle
            if not line.startswith(">")
        )

    return header_id, sequence


# ----------------------------
# Parse observed TELROI tags
# ----------------------------
observed_rows = []

fasta_files = sorted(
    TAGGED_DIR.glob("*_TELROI_TAGGED.fasta")
)

if not fasta_files:
    raise FileNotFoundError(
        f"No tagged chromosome FASTA files found in: {TAGGED_DIR}"
    )

for fasta_path in fasta_files:
    header_id, tagged_sequence = read_single_record_fasta(fasta_path)
    tagged_sequence_length = len(tagged_sequence)

    for match in TELROI_PATTERN.finditer(tagged_sequence):
        tag_data = match.groupdict()

        tag_data["actual_roi_len"] = int(
            tag_data["actual_roi_len"]
        )
        tag_data["roi_start_1b"] = int(
            tag_data["roi_start_1b"]
        )
        tag_data["roi_end_1b"] = int(
            tag_data["roi_end_1b"]
        )

        tag_start_0b = match.start()
        tag_end_0b = match.end()

        chromosome_end = tag_data["end"]

        if chromosome_end == "5p":
            tag_at_expected_edge = tag_start_0b == 0
        else:
            tag_at_expected_edge = (
                tag_end_0b == tagged_sequence_length
            )

        coordinate_span_length = (
            tag_data["roi_end_1b"]
            - tag_data["roi_start_1b"]
            + 1
        )

        span_length_matches_roi = (
            coordinate_span_length
            == tag_data["actual_roi_len"]
        )

        observed_rows.append(
            {
                "file": fasta_path.name,
                "chrom_header_id": header_id,
                "tag_pos_start_1b_in_tagged_seq": (
                    tag_start_0b + 1
                ),
                "tag_pos_end_1b_in_tagged_seq": tag_end_0b,
                "tagged_seq_len_characters": tagged_sequence_length,
                "tag_at_expected_edge": tag_at_expected_edge,
                "span_len_bp": coordinate_span_length,
                "span_len_matches_actual_roi": (
                    span_length_matches_roi
                ),
                "chrom_id_matches_header": (
                    header_id == tag_data["chrom"]
                ),
                **tag_data,
            }
        )


tags = pd.DataFrame(observed_rows)

if tags.empty:
    raise ValueError(
        "No TELROI placeholders were detected in the tagged FASTAs. "
        "Check the placeholder format and tagged FASTA files."
    )


# ----------------------------
# Validate observed columns
# ----------------------------
observed_required_columns = [
    "sample",
    "chrom",
    "end",
    "tel_class",
    "actual_roi_len",
    "roi_start_1b",
    "roi_end_1b",
]

missing_observed_columns = [
    column
    for column in observed_required_columns
    if column not in tags.columns
]

if missing_observed_columns:
    raise ValueError(
        "Observed TELROI tags are missing fields: "
        + ", ".join(missing_observed_columns)
    )


# ----------------------------
# Compare observed tags with table
# ----------------------------
merge_columns = [
    "sample",
    "chrom",
    "end",
    "tel_class",
    "actual_roi_len",
    "roi_start_1b",
    "roi_end_1b",
]

merged = tags.merge(
    roi,
    on=merge_columns,
    how="left",
    indicator=True,
    suffixes=("", "_expected"),
)

merged["tag_matches_table"] = merged["_merge"].eq("both")

merged["roi_within_bounds_ok"] = (
    merged["roi_start_1b"].ge(1)
    & merged["roi_end_1b"].le(merged["chromosome_length"])
    & merged["roi_start_1b"].le(merged["roi_end_1b"])
)


# ----------------------------
# Detailed sanity table
# ----------------------------
output_columns = [
    "file",
    "chrom_header_id",
    "chrom",
    "sample",
    "end",
    "tel_class",
    "roi_start_1b",
    "roi_end_1b",
    "actual_roi_len",
    "span_len_bp",
    "span_len_matches_actual_roi",
    "tag_at_expected_edge",
    "roi_within_bounds_ok",
    "chrom_id_matches_header",
    "tag_matches_table",
]

sanity_table = (
    merged[output_columns]
    .sort_values(["sample", "chrom", "end"])
    .reset_index(drop=True)
)

display(sanity_table)


# ----------------------------
# Summary table
# ----------------------------
summary = pd.DataFrame(
    [
        {
            "n_tagged_fasta_files": len(fasta_files),
            "n_observed_tags": len(tags),
            "n_expected_tags_in_table": len(roi),
            "n_tags_matched_to_table": int(
                merged["tag_matches_table"].sum()
            ),
            "n_tags_unmatched": int(
                (~merged["tag_matches_table"]).sum()
            ),
            "n_wrong_edge": int(
                (~merged["tag_at_expected_edge"]).sum()
            ),
            "n_roi_length_mismatch": int(
                (~merged["span_len_matches_actual_roi"]).sum()
            ),
            "n_bounds_fail": int(
                (~merged["roi_within_bounds_ok"]).sum()
            ),
            "n_header_mismatch": int(
                (~merged["chrom_id_matches_header"]).sum()
            ),
        }
    ]
)

display(summary)


# ----------------------------
# Display failed rows
# ----------------------------
failed_rows = sanity_table[
    (~sanity_table["tag_matches_table"])
    | (~sanity_table["tag_at_expected_edge"])
    | (~sanity_table["span_len_matches_actual_roi"])
    | (~sanity_table["roi_within_bounds_ok"])
    | (~sanity_table["chrom_id_matches_header"])
]

print(f"\nFAIL rows: {len(failed_rows)}")
display(failed_rows)


# ----------------------------
# Save outputs
# ----------------------------
sanity_table.to_csv(
    OUT_TSV,
    sep="\t",
    index=False,
)

summary.to_csv(
    OUT_SUMMARY_TSV,
    sep="\t",
    index=False,
)

print(f"[OK] Sanity table written to: {OUT_TSV}")
print(f"[OK] Summary written to: {OUT_SUMMARY_TSV}")
