#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 13_filter_best_hits_by_rescue_list.py
#
# Restricts the top-hits TSV from script 12 to only those (sample, chromosome, end) rows
# that appear in the rescue list, i.e. the chromosome ends selected for rescue after manual
# review. Everything else is dropped.
# ------------------------------------------------------------------------------
"""
Filter a best-hits TSV to only rows whose (sample, chromosome, end) are in a rescue list.
Does not modify the selector script; use after select_two_best_hits_per_roi*.py.
Rescue list TSV: SampleID, Chromosome, end (e.g. 1268  chr01  3p)
Best-hits TSV: must have sample (or SampleID), chrom (or chromosome), end.
Usage:
  python filter_best_hits_by_rescue_list.py \\
    --input telomeres_two_best_hits.tsv \\
    --rescue-list telomeres_rescue_list.tsv \\
    --output telomeres_two_best_hits_filtered.tsv
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Optional
ROOT = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling")
def extract_chr_token(s: str) -> Optional[str]:
    """Normalize to chrNN or chrX/chrY from strings like 'chr01', '1268_chr01_3p_TELROI'."""
    if not s:
        return None
    m = re.search(r"(?:^|[_\W])chr(\d+)(?:[_\W]|$)", s)
    if m:
        return f"chr{int(m.group(1)):02d}"
    if re.search(r"(?:^|[_\W])chr([XY])(?:[_\W]|$)", s, re.IGNORECASE):
        return re.search(r"chr([XY])", s, re.IGNORECASE).group(0)
    return None
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--input", type=Path, required=True, help="Best-hits TSV (from selector)")
    ap.add_argument("--rescue-list", type=Path, default=ROOT / "telomeres_rescue_list.tsv", help="Rescue list TSV (SampleID, Chromosome, end)")
    ap.add_argument("--output", type=Path, required=True, help="Filtered TSV output")
    args = ap.parse_args()
    if not args.rescue_list.exists():
        print(f"[ERROR] Rescue list not found: {args.rescue_list}", file=sys.stderr)
        sys.exit(1)
    if not args.input.exists():
        print(f"[ERROR] Input TSV not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    # Load rescue set: (sample, chr_tok, end)
    rescue_set: set[tuple[str, str, str]] = set()
    with open(args.rescue_list, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            sample = (row.get("SampleID") or "").strip()
            chrom_raw = (row.get("Chromosome") or "").strip()
            end = (row.get("end") or "").strip()
            if not sample or end not in {"5p", "3p"}:
                continue
            tok = extract_chr_token(chrom_raw) or chrom_raw
            rescue_set.add((sample, tok, end))
    print(f"[INFO] Rescue list: {len(rescue_set)} (sample, chr, end) entries from {args.rescue_list.name}")
    # Filter best-hits TSV
    kept = 0
    skipped = 0
    out_fieldnames = None
    out_rows = []
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        out_fieldnames = list(reader.fieldnames or [])
        for row in reader:
            sample = (row.get("sample") or row.get("SampleID") or "").strip()
            chrom_raw = (row.get("chrom") or row.get("chromosome") or "").strip()
            end = (row.get("end") or "").strip()
            tok = extract_chr_token(chrom_raw) or extract_chr_token(row.get("paf_file") or "")
            if not tok:
                skipped += 1
                continue
            if (sample, tok, end) in rescue_set:
                out_rows.append(row)
                kept += 1
            else:
                skipped += 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    print(f"[INFO] Wrote {kept} rows (skipped {skipped}) to {args.output}")
if __name__ == "__main__":
    main()
