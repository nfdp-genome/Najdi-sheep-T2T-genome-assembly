#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 12_select_two_best_hits_per_roi.py
#
# Reads the aggregated all-hits TSV compiled from the per-ROI PAFs (script 11) and selects
# the top 1 or 2 candidate hits per TELROI. Ranks by: distance from contig end (asc), alignment
# length (desc), sequence identity (desc). Requires each candidate to reach the telomere-facing
# end of the query (terminal_bp), pass a minimum alignment length, and (optionally) a min/max
# overhang length to avoid internal telomere-like repeats.
# ------------------------------------------------------------------------------
"""
Select two best hits per ROI from the all-hits TSV. V1.1.

Criteria:
- Target terminality: overhang at end of contig. 5p → tstart ≤ terminal_bp; 3p → (tlen - tend) ≤ terminal_bp.
- Min extension: overhang at contig end must be ≥ min_extension (5p: tstart; 3p: tlen-tend).
- Max extension (optional): reject if overhang > max_extension to avoid internal telomere-like repeats.
- Allow tp=P and tp=S (no filter on tp).
- qcov not used for ranking (repeats). Rank by: distance from contig end (asc), then alnlen (desc), then identity (desc).
- Pick top 1 or 2 hits per ROI (--top). If 2, prefer second on a different tname than the first.
- Optional --ensure-one-per-roi: for any ROI with no hit after filters, add the best available hit (relaxed) and set rescue=1.

Output: same columns as input plus hit_rank (1 or 2). If --ensure-one-per-roi used, plus rescue (0 or 1).

Usage:
  python select_two_best_hits_per_roi_V1.1.py --input telroi_allhits_unfiltered.tsv --output telomeres_two_best_hits.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

ROOT = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument(
        "--input",
        type=Path,
        default=ROOT / "telroi_allhits_unfiltered.tsv",
        help="All-hits TSV (from PAF)",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=ROOT / "telomeres_two_best_hits.tsv",
        help="Output TSV with hit_rank column",
    )
    ap.add_argument(
        "--terminal-bp",
        type=int,
        default=500,
        help="Max distance from contig end for target terminality (5p: tstart; 3p: tlen-tend). Default 500.",
    )
    ap.add_argument(
        "--min-aln",
        type=int,
        default=20_000,
        help="Min alignment length (bp). Default 20000.",
    )
    ap.add_argument(
        "--min-qcov",
        type=float,
        default=0.0,
        help="Min query coverage (0 = no filter). Default 0.",
    )
    ap.add_argument(
        "--min-extension",
        type=int,
        default=0,
        help="Min overhang (bp) at contig end: 5p = tstart, 3p = tlen-tend. Default 0.",
    )
    ap.add_argument(
        "--max-extension",
        type=int,
        default=0,
        metavar="BP",
        help="Max overhang (bp); reject if larger (avoids internal repeats). 0 = no cap. Default 0.",
    )
    ap.add_argument(
        "--top",
        type=int,
        default=2,
        choices=(1, 2),
        help="Number of best hits per ROI to keep (1 or 2). Default 2.",
    )
    ap.add_argument(
        "--ensure-one-per-roi",
        action="store_true",
        help="If set, every ROI gets at least one hit; missing ROIs get best available (rescue=1).",
    )
    args = ap.parse_args()

    if not args.input.exists():
        print(f"[ERROR] Input not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    rows_by_roi: dict[str, list[dict[str, Any]]] = {}
    all_rows_by_roi: dict[str, list[dict[str, Any]]] = {}  # for rescue when --ensure-one-per-roi

    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            paf_file = (row.get("paf_file") or "").strip()
            end = (row.get("end") or "").strip()
            if not paf_file or end not in {"5p", "3p"}:
                continue

            try:
                tstart = int(row.get("tstart", 0))
                tend = int(row.get("tend", 0))
                tlen = int(row.get("tlen", 0))
                alnlen = int(row.get("alnlen", 0))
            except (ValueError, TypeError):
                continue
            try:
                qcov = float(row.get("qcov", 0) or 0)
            except (ValueError, TypeError):
                qcov = 0.0

            # Compute dist_end for ranking (and for all_rows when rescue)
            if end == "5p":
                dist_end = tstart
            else:  # 3p
                dist_end = tlen - tend
                if dist_end < 0:
                    continue

            try:
                identity = float(row.get("identity", 0) or 0)
            except (ValueError, TypeError):
                identity = 0.0

            row_copy = dict(row)
            row_copy["_dist_end"] = dist_end
            row_copy["_alnlen"] = alnlen
            row_copy["_identity"] = identity

            if args.ensure_one_per_roi:
                all_rows_by_roi.setdefault(paf_file, []).append(row_copy)

            # Strict filters for main selection
            if end == "5p":
                if tstart > args.terminal_bp:
                    continue
            else:
                if dist_end > args.terminal_bp:
                    continue
            if dist_end < args.min_extension:
                continue
            if args.max_extension > 0 and dist_end > args.max_extension:
                continue
            if alnlen < args.min_aln:
                continue
            if qcov < args.min_qcov:
                continue

            row["_dist_end"] = dist_end
            row["_alnlen"] = alnlen
            row["_identity"] = identity
            rows_by_roi.setdefault(paf_file, []).append(row)

    # Sort each group: dist_end asc, alnlen desc, identity desc
    for roi in rows_by_roi:
        rows_by_roi[roi].sort(
            key=lambda r: (r["_dist_end"], -r["_alnlen"], -r["_identity"])
        )

    # Pick top 1 or 2 per ROI (rank 2 on different tname when possible)
    out_rows: list[dict[str, Any]] = []
    for roi, group in sorted(rows_by_roi.items()):
        if not group:
            continue
        first = group[0]
        first["hit_rank"] = 1
        # Remove internal keys before output
        for k in ("_dist_end", "_alnlen", "_identity"):
            first.pop(k, None)
        out_rows.append(first)

        if args.top < 2 or len(group) < 2:
            continue
        first_tname = (first.get("tname") or "").strip()
        second = None
        for r in group[1:]:
            if (r.get("tname") or "").strip() != first_tname:
                second = r
                break
        if second is None:
            second = group[1]
        second["hit_rank"] = 2
        for k in ("_dist_end", "_alnlen", "_identity"):
            second.pop(k, None)
        out_rows.append(second)

    # Ensure every ROI has at least one hit (rescue with best available)
    if args.ensure_one_per_roi and all_rows_by_roi:
        out_rois = {r.get("paf_file", "").strip() for r in out_rows}
        for roi in sorted(all_rows_by_roi.keys()):
            if roi in out_rois:
                continue
            group = all_rows_by_roi[roi]
            group.sort(key=lambda r: (r["_dist_end"], -r["_alnlen"], -r["_identity"]))
            first = dict(group[0])
            first["hit_rank"] = 1
            first["rescue"] = 1
            for k in ("_dist_end", "_alnlen", "_identity"):
                first.pop(k, None)
            out_rows.append(first)
            out_rois.add(roi)
            if args.top >= 2 and len(group) >= 2:
                first_tname = (first.get("tname") or "").strip()
                second = None
                for r in group[1:]:
                    if (r.get("tname") or "").strip() != first_tname:
                        second = r
                        break
                if second is None:
                    second = group[1]
                second = dict(second)
                second["hit_rank"] = 2
                second["rescue"] = 1
                for k in ("_dist_end", "_alnlen", "_identity"):
                    second.pop(k, None)
                out_rows.append(second)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_fieldnames = fieldnames + ["hit_rank"]
    if args.ensure_one_per_roi:
        out_fieldnames = fieldnames + ["hit_rank", "rescue"]
        for r in out_rows:
            if "rescue" not in r:
                r["rescue"] = 0

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    n_roi = len(rows_by_roi)
    n_rescue = len(all_rows_by_roi) - n_roi if args.ensure_one_per_roi else 0
    n_two = sum(1 for r in out_rows if r.get("hit_rank") == 2)
    n_rescue_rows = sum(1 for r in out_rows if r.get("rescue") == 1)
    msg = f"Wrote {len(out_rows)} rows ({n_roi} ROIs, top={args.top}, {n_two} with rank 2)"
    if args.ensure_one_per_roi:
        msg += f"; {n_rescue} ROIs rescued ({n_rescue_rows} rescue rows)"
    msg += f" to {out_path}"
    print(msg)


if __name__ == "__main__":
    main()
