#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 03_tag_chromosome_fastas_with_telroi.py
#
# Insert @@TELROI|sample|chrom|end|class|actual_roi_len|roi:start-end@@ placeholders on a
# working copy of each selected chromosome FASTA (5' before the sequence, 3' after).
# ------------------------------------------------------------------------------
# STEP 2: Add TELROI placeholders to chromosome FASTA files
# - Reads telomere_roi_table.tsv
# - Loads the chromosome assemblies for samples 1268 and 1271
# - Adds placeholders at the selected 5p and/or 3p chromosome ends
# - Writes one tagged FASTA file per chromosome
# - Performs no sequence extraction or replacement

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ----------------------------
# Paths
# ----------------------------
BASE_DIR = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)

ROI_TABLE = BASE_DIR / "telomere_roi_table.tsv"
OUT_DIR = BASE_DIR / "tagged_chromosomes"

OUT_DIR.mkdir(parents=True, exist_ok=True)


FASTA_BY_SAMPLE = {
    "1268": Path(
        "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
        "hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/"
        "unplaced_contigs_vs_scaffolded_assembly/"
        "1268.Hifiasm_HiC_salsa_racon_cleaned_clean_ChrDefline_scaffolded.fasta"
    ),
    "1271": Path(
        "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
        "hifiasm_assembly/chromosomal_reconstruction/unplaced_contigs_investigation/"
        "unplaced_contigs_vs_scaffolded_assembly/"
        "1271.Hifiasm_HiC_yahs_racon_tgsgapcloser_cleaned_ChrDefline_scaffolded.fasta"
    ),
}


# ----------------------------
# Load TELROI decision table
# ----------------------------
df = pd.read_csv(ROI_TABLE, sep="\t")

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
        "Missing required columns from TELROI table: "
        + ", ".join(sorted(missing_columns))
    )


# ----------------------------
# Build TELROI placeholder
# ----------------------------
def make_telroi_placeholder(row):
    """
    Build a TELROI placeholder without FASTA header characters.
    """
    return (
        f"@@TELROI|{row['sample']}|{row['chromosome']}|"
        f"{row['end']}|{row['tel_class']}|"
        f"{int(row['actual_roi_len'])}|"
        f"roi:{int(row['roi_start_1b'])}-{int(row['roi_end_1b'])}@@"
    )


# ----------------------------
# Tag chromosomes
# ----------------------------
for sample, fasta_path in FASTA_BY_SAMPLE.items():
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"Assembly FASTA not found for sample {sample}: {fasta_path}"
        )

    # Load all chromosome records once for this sample
    records = {
        record.id: record
        for record in SeqIO.parse(str(fasta_path), "fasta")
    }

    sample_rows = df[
        df["sample"].astype(str) == str(sample)
    ].copy()

    for chromosome, group in sample_rows.groupby("chromosome"):
        chromosome = str(chromosome)

        if chromosome not in records:
            raise KeyError(
                f"Chromosome '{chromosome}' was not found in: {fasta_path}"
            )

        sequence = str(records[chromosome].seq)

        # Process the 5′ placeholder before the 3′ placeholder
        end_order = {"5p": 0, "3p": 1}

        group = group.sort_values(
            by="end",
            key=lambda values: values.map(end_order),
        )

        for _, row in group.iterrows():
            chromosome_end = row["end"]

            if chromosome_end not in end_order:
                raise ValueError(
                    f"Unexpected chromosome end '{chromosome_end}' for "
                    f"{sample} {chromosome}. Expected '5p' or '3p'."
                )

            placeholder = make_telroi_placeholder(row)

            if chromosome_end == "5p":
                sequence = placeholder + sequence
            else:
                sequence = sequence + placeholder

        # Remove embedded newline characters before writing
        sequence = sequence.replace("\n", "").replace("\r", "")

        output_fasta = (
            OUT_DIR / f"{sample}_{chromosome}_TELROI_TAGGED.fasta"
        )

        output_record = SeqRecord(
            Seq(sequence),
            id=chromosome,
            description="TELROI_TAGGED",
        )

        SeqIO.write(
            output_record,
            str(output_fasta),
            "fasta",
        )

        print(
            f"[TAGGED_OK] {sample} {chromosome} -> {output_fasta}"
        )
