#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 07_filter_paf_and_rank_gap_candidates.py
#
# TELROI + GAP ROI PAF filtering & candidate ranking. Reads the aggregated PAF from mapping
# gap ROIs against the merged assembly pool, ranks candidates that span the full gap plus both
# flanking regions, and writes the per-gap best-hit table used by the patching steps below.
# ------------------------------------------------------------------------------
# ============================================================
# TELROI + GAP ROI PAF FILTERING & SANITY PREP (ROBUST BLOCK)
# ============================================================

import pandas as pd
import numpy as np
import re

# ----------------------------
# USER INPUT
# ----------------------------
PAF_PATH = "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/mapping_to_mergedAssemblies/1268/chr02_gap_138435449_138435836/1268.chr02_gap_138435449_138435836.ROIcombined_vs_mergedAssemblies.paf"
TELROI_TABLE = "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling/telomere_roi_table.tsv"

MAPQ_GAP_MIN = 20
MAPQ_TEL_MIN = 5
MIN_QCOV = 0.40
TOP_N_HITS = 10

# ----------------------------
# LOAD PAF (first 12 cols only)
# ----------------------------
cols = [
    "qname","qlen","qstart","qend","strand",
    "tname","tlen","tstart","tend",
    "aln_bp","block_len","mapq"
]

# Read full line, but only keep first 12 cols; ignore minimap tags
paf = pd.read_csv(
    PAF_PATH,
    sep="\t",
    header=None,
    usecols=list(range(12)),
    names=cols,
    dtype=str,            # read as str first (safe), then coerce numeric
    engine="python"
)

# force numeric for the integer columns
num_cols = ["qlen","qstart","qend","tlen","tstart","tend","aln_bp","block_len","mapq"]
for c in num_cols:
    paf[c] = pd.to_numeric(paf[c], errors="coerce")

# drop any rows that failed coercion (corrupt lines)
bad = paf[num_cols].isna().any(axis=1)
if bad.any():
    print(f"[WARN] Dropping {bad.sum()} malformed PAF rows (non-numeric in required fields).")
    paf = paf.loc[~bad].copy()

# metrics
#paf["qcov"] = paf["aln_bp"] / paf["qlen"]
paf["qcov"] = paf["block_len"] / paf["qlen"]   # ✅ use aligned block length, not matches
paf["is_telroi"] = paf["qname"].str.startswith("@@TELROI")

# ----------------------------
# TELROI PARSER
# ----------------------------
def parse_telroi_qname(q):
    """
    @@TELROI|sample|chrom|end|class|anchor_len|roi:start-end@@
    """
    m = re.match(
        r"@@TELROI\|(\d+)\|([^|]+)\|(5p|3p)\|([A-Z])\|(\d+)\|roi:(\d+)-(\d+)@@",
        q
    )
    if not m:
        return None
    return {
        "sample": m.group(1),
        "chrom": m.group(2),
        "end": m.group(3),
        "tel_class": m.group(4),
        "anchor_len_in_tag": int(m.group(5)),
        "roi_start_in_tag_1b": int(m.group(6)),
        "roi_end_in_tag_1b": int(m.group(7)),
    }

# ----------------------------
# PARSE TELROI QNAMES SAFELY
# ----------------------------
tel_parsed = paf.loc[paf["is_telroi"], ["qname"]].drop_duplicates().copy()
tel_parsed["parsed"] = tel_parsed["qname"].apply(parse_telroi_qname)
tel_parsed = tel_parsed[tel_parsed["parsed"].notna()]

if tel_parsed.empty:
    tel_meta = pd.DataFrame(columns=[
        "qname","sample","chrom","end","tel_class",
        "anchor_len_in_tag","roi_start_in_tag_1b","roi_end_in_tag_1b"
    ])
else:
    tel_meta = pd.concat(
        [
            tel_parsed[["qname"]].reset_index(drop=True),
            tel_parsed["parsed"].apply(pd.Series).reset_index(drop=True)
        ],
        axis=1
    )

# merge TELROI metadata (left join, only telroi rows get filled)
paf = paf.merge(
    tel_meta[
        ["qname","sample","chrom","end","tel_class",
         "anchor_len_in_tag","roi_start_in_tag_1b","roi_end_in_tag_1b"]
    ],
    on="qname",
    how="left"
)

# ----------------------------
# GAP ROI SAMPLE DETECTION
# ----------------------------
paf.loc[~paf["is_telroi"], "sample"] = paf.loc[
    ~paf["is_telroi"], "qname"
].str.extract(r"^(\d+)")

# ----------------------------
# LOAD TELROI TABLE (robust whitespace-separated header)
# ----------------------------
tel_tbl = pd.read_csv(TELROI_TABLE, sep=r"\s+", engine="python")
tel_tbl = tel_tbl[["sample","chromosome","chromosome_length","anchor_len"]].drop_duplicates()

# make sure sample is string to match parsed 'sample'
tel_tbl["sample"] = tel_tbl["sample"].astype(str)

# attach chr length + anchor_len (table) for telroi positional sanity
paf = paf.merge(
    tel_tbl,
    left_on=["sample","tname"],
    right_on=["sample","chromosome"],
    how="left"
)

# ----------------------------
# TELROI POSITIONAL CORRECTNESS
# Use anchor_len from telomere_roi_table.tsv (your rule)
# ----------------------------
paf["dist_to_start"] = paf["tstart"]
paf["dist_to_end"] = paf["chromosome_length"] - paf["tend"]

paf["tel_position_ok"] = np.where(
    paf["is_telroi"],
    (
        ((paf["end"] == "5p") & (paf["dist_to_start"] <= paf["anchor_len"])) |
        ((paf["end"] == "3p") & (paf["dist_to_end"]   <= paf["anchor_len"]))
    ),
    pd.NA
)

# ----------------------------
# FILTERING LOGIC (your agreed rules)
# ----------------------------

# GAP ROI filtering
gap_filt = paf[
    (~paf["is_telroi"]) &
    (paf["mapq"] >= MAPQ_GAP_MIN) &
    (paf["qcov"] >= MIN_QCOV)
]

# TELROI filtering
tel_filt = paf[
    (paf["is_telroi"]) &
    (paf["mapq"] >= MAPQ_TEL_MIN) &
    (paf["tel_position_ok"] == True)
]

# combine
paf_filt = pd.concat([gap_filt, tel_filt], ignore_index=True)

# ----------------------------
# BEST 10 HITS PER QUERY
# rank by aln_bp desc, tie-break mapq desc (your rule)
# ----------------------------
paf_filt = (
    paf_filt
    .sort_values(["qname","aln_bp","mapq"], ascending=[True, False, False])
    .groupby("qname", group_keys=False)
    .head(TOP_N_HITS)
    .reset_index(drop=True)
)

# ----------------------------
# QUICK SUMMARY + PREVIEW
# ----------------------------
print("=== SUMMARY ===")
print(f"Total PAF rows        : {len(paf)}")
print(f"TELROI rows           : {int(paf['is_telroi'].sum())}")
print(f"GAP ROI rows          : {int((~paf['is_telroi']).sum())}")
print(f"Filtered rows kept    : {len(paf_filt)}")
print(f"Filtered TELROI hits  : {int(paf_filt['is_telroi'].sum())}")
print(f"Filtered GAP hits     : {int((~paf_filt['is_telroi']).sum())}")

display(
    paf_filt[
        ["qname","tname","tstart","tend","aln_bp","qcov","mapq",
         "is_telroi","tel_position_ok"]
    ].head(12)
)
