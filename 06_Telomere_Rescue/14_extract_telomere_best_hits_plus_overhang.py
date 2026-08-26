#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 14_extract_telomere_best_hits_plus_overhang.py
#
# For each surviving best-hit from script 13, extracts four FASTA products from the
# corresponding candidate contig: (1) full contig sequence, (2) aligned block only, (3) aligned
# block + telomere-side overhang (legacy), and (4) the extension-only segment (the telomeric
# overhang beyond the aligned block). The extension-only FASTA (product 4) is the input to
# script 15.
# ------------------------------------------------------------------------------
"""
Extract telomere best-hit regions plus telomeric overhang AND (optionally) the
extension-only segment (for Option A: append-only telomeric extension).

Inputs:
- Best-hit TSV (default: telomeres_best_hit_batch5.tsv)
- Optional rescue-list TSV to restrict which chromosome ends are processed:
    columns: SampleID, Chromosome, end
    e.g. 1268  chr01  3p

Outputs:
1) telomeres_best_hits_fullSeq/        -> full contig sequence (entire tname)
2) telomeres_best_hits_matchingCoord/  -> aligned block only [tstart, tend]
3) telomeres_best_hits_extension_only/ -> extension-only segment (overhang beyond aligned block)
4) telomeres_best_hits/                -> aligned+overhang (legacy: matchingCoord + extension)

Notes (Option A planning):
- extension-only is the segment outside the aligned block, toward the contig end:
    * end=5p: [tstart - tel_ext_bp, tstart)
    * end=3p: [tend, tend + tel_ext_bp)
- reverse-complements extension (and full extract) when strand is '-'.

Usage:
  python extract_telomere_best_hits_V1.2.py \
      --tsv /path/telomeres_best_hit_batch5.tsv \
      --rescue-list /path/telomeres_rescue_list.tsv

Dependencies: pyfaidx, biopython
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from Bio.Seq import Seq
from pyfaidx import Fasta

ROOT = Path("/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling")
ASSEMBLIES = {
    "1268": ROOT / "1268_all_assemblies.merged.unique.fasta",
    "1271": ROOT / "1271_all_assemblies.merged.unique.fasta",
}

# How close must alignment touch contig boundary to be "terminal" (sanity flag)
TERMINAL_BP_TOL = 10


def roi_id_from_paf_file(paf_file: str) -> str:
    """Strip .vs_mergedAssemblies.paf to get ROI ID."""
    if paf_file.endswith(".vs_mergedAssemblies.paf"):
        return paf_file[:-len(".vs_mergedAssemblies.paf")]
    return paf_file.replace(".paf", "")


def extract_chr_token(s: str) -> Optional[str]:
    """
    Extract chr token like 'chr01' from strings such as:
      - 'NLFDP1271_chr25'
      - 'chr25'
      - '1271_chr25_5p'
      - '1268_NLFDP1268_chr01_3p_TELROI'
    Returns normalized 'chrNN' with zero padding for numeric chroms; passes through chrX/chrY.
    """
    if not s:
        return None
    # Numeric: allow chr preceded by start or underscore (ROI IDs use _ before chr)
    m = re.search(r"(?:^|[_\W])chr(\d+)(?:[_\W]|$)", s)
    if m:
        n = int(m.group(1))
        return f"chr{n:02d}"
    # Non-numeric (e.g. chrX, chrY): require whole-token match
    if re.search(r"(?:^|[_\W])chr([XY])(?:[_\W]|$)", s, re.IGNORECASE):
        return re.search(r"chr([XY])", s, re.IGNORECASE).group(0)  # e.g. chrX
    return None


def load_rescue_list(path: Path) -> Set[Tuple[str, str, str]]:
    """
    Load rescue list rows as (sample, chrNN, end).
    Expected columns: SampleID, Chromosome, end
    """
    wanted: Set[Tuple[str, str, str]] = set()
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in ("SampleID", "Chromosome", "end") if c not in reader.fieldnames]
        if missing:
            raise ValueError(f"Rescue-list TSV missing columns: {missing}. Found: {reader.fieldnames}")
        for row in reader:
            sample = (row.get("SampleID") or "").strip()
            chrom_raw = (row.get("Chromosome") or "").strip()
            end = (row.get("end") or "").strip()
            if not sample or not chrom_raw or end not in {"5p", "3p"}:
                continue
            chr_tok = extract_chr_token(chrom_raw) or chrom_raw
            wanted.add((sample, chr_tok, end))
    return wanted


def is_terminal_on_contig(end: str, tstart: int, tend: int, tlen: int) -> bool:
    """
    Sanity: does the aligned block touch the contig boundary on the side you care about?
    This is not a pass/fail for extraction, just a flag in defline.
    - For end=5p: aligned start near 0 suggests extension is truly terminal.
    - For end=3p: aligned end near tlen suggests extension is truly terminal.
    """
    if end == "5p":
        return tstart <= TERMINAL_BP_TOL
    if end == "3p":
        return (tlen - tend) <= TERMINAL_BP_TOL
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--tsv", default=ROOT / "telomeres_best_hit_batch5.tsv", type=Path, help="Best-hit TSV")
    ap.add_argument("--out-dir", default=ROOT / "telomeres_best_hits", type=Path, help="Output dir for aligned+overhang (legacy)")
    ap.add_argument(
        "--out-dir-full-seq",
        default=ROOT / "telomeres_best_hits_fullSeq",
        type=Path,
        help="Output dir for full contig sequence per hit",
    )
    ap.add_argument(
        "--out-dir-matching-coord",
        default=ROOT / "telomeres_best_hits_matchingCoord",
        type=Path,
        help="Output dir for aligned-block-only (matching coordinates) FASTAs",
    )
    ap.add_argument(
        "--out-dir-ext-only",
        default=ROOT / "telomeres_best_hits_extension_only",
        type=Path,
        help="Output dir for extension-only FASTAs",
    )
    ap.add_argument(
        "--rescue-list",
        default=ROOT / "telomeres_rescue_list.tsv",
        type=Path,
        help="Optional: TSV restricting which SampleID/Chromosome/end to process. If missing, no filtering.",
    )
    ap.add_argument(
        "--disable-rescue-filter",
        action="store_true",
        help="If set, ignores --rescue-list even if present.",
    )
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir_full_seq.mkdir(parents=True, exist_ok=True)
    args.out_dir_matching_coord.mkdir(parents=True, exist_ok=True)
    args.out_dir_ext_only.mkdir(parents=True, exist_ok=True)

    if not args.tsv.exists():
        print(f"[ERROR] TSV not found: {args.tsv}", file=sys.stderr)
        sys.exit(1)

    rescue_set: Optional[Set[Tuple[str, str, str]]] = None
    if not args.disable_rescue_filter and args.rescue_list.exists():
        rescue_set = load_rescue_list(args.rescue_list)
        print(f"[INFO] Rescue filter enabled: {len(rescue_set)} entries from {args.rescue_list.name}")
    elif not args.disable_rescue_filter and not args.rescue_list.exists():
        print(f"[INFO] Rescue list not found at {args.rescue_list} -> no filtering applied")
    else:
        print("[INFO] Rescue filter disabled by --disable-rescue-filter")

    n_ok = 0
    n_full_ok = 0
    n_match_ok = 0
    n_ext_ok = 0
    n_skip = 0
    n_fail = 0
    fa_cache: Dict[str, Fasta] = {}

    with open(args.tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            paf_file = (row.get("paf_file") or "").strip()
            sample = (row.get("sample") or "").strip()
            tname = (row.get("tname") or "").strip()
            end = (row.get("end") or "").strip()
            strand = (row.get("strand") or "+").strip()

            # Derive ROI ID and chromosome token for filtering/output clarity
            roi_id = roi_id_from_paf_file(paf_file)
            chrom_raw = (row.get("chrom", "") or row.get("chromosome", "") or "").strip()
            chr_tok = (
                extract_chr_token(chrom_raw)
                or extract_chr_token(roi_id or "")
                or extract_chr_token(tname or "")
            )

            # Apply rescue filter if available
            if rescue_set is not None:
                if not chr_tok:
                    print(f"[SKIP] {roi_id}: could not derive chromosome token for rescue filtering", file=sys.stderr)
                    n_skip += 1
                    continue
                if (sample, chr_tok, end) not in rescue_set:
                    n_skip += 1
                    continue

            # Parse ints
            try:
                tstart = int(row.get("tstart", 0))
                tend = int(row.get("tend", 0))
                tlen = int(row.get("tlen", 0))
            except ValueError:
                print(f"[SKIP] {roi_id}: invalid coordinates", file=sys.stderr)
                n_skip += 1
                continue

            try:
                tel_ext_bp = int(row.get("tel_ext_bp", 0) or 0)
            except ValueError:
                print(f"[SKIP] {roi_id}: invalid tel_ext_bp", file=sys.stderr)
                n_skip += 1
                continue

            if not roi_id:
                print(f"[SKIP] {paf_file}: could not derive ROI ID", file=sys.stderr)
                n_skip += 1
                continue

            ref_path = ASSEMBLIES.get(sample)
            if not ref_path or not ref_path.exists():
                print(f"[SKIP] {roi_id}: assembly not found for sample {sample}", file=sys.stderr)
                n_skip += 1
                continue

            if end not in {"5p", "3p"}:
                print(f"[SKIP] {roi_id}: invalid end={end}", file=sys.stderr)
                n_skip += 1
                continue

            if tstart >= tend or tstart < 0 or tend < 0:
                print(f"[SKIP] {roi_id}: invalid range tstart={tstart} tend={tend}", file=sys.stderr)
                n_skip += 1
                continue
            if tend > tlen:
                print(f"[SKIP] {roi_id}: tend={tend} > tlen={tlen}", file=sys.stderr)
                n_skip += 1
                continue

            # --- 1) aligned + overhang (legacy) ---
            ext_start = tstart
            ext_end = tend
            if tel_ext_bp > 0:
                if end == "5p":
                    ext_start = max(0, tstart - tel_ext_bp)
                else:  # 3p
                    ext_end = min(tlen, tend + tel_ext_bp)

            if ext_start >= ext_end:
                print(f"[SKIP] {roi_id}: invalid extended range ext_start={ext_start} ext_end={ext_end}", file=sys.stderr)
                n_skip += 1
                continue

            # --- 2) extension-only (Option A) ---
            # This is ONLY the extra piece beyond the aligned block.
            ext_only_start = None
            ext_only_end = None
            if tel_ext_bp > 0:
                if end == "5p":
                    ext_only_start = max(0, tstart - tel_ext_bp)
                    ext_only_end = tstart
                else:  # 3p
                    ext_only_start = tend
                    ext_only_end = min(tlen, tend + tel_ext_bp)

                if ext_only_start >= ext_only_end:
                    ext_only_start = None
                    ext_only_end = None

            terminal_flag = is_terminal_on_contig(end=end, tstart=tstart, tend=tend, tlen=tlen)

            try:
                if sample not in fa_cache:
                    fa_cache[sample] = Fasta(str(ref_path), one_based_attributes=False)
                fa = fa_cache[sample]

                if tname not in fa:
                    print(f"[SKIP] {roi_id}: contig {tname} not in assembly", file=sys.stderr)
                    n_skip += 1
                    continue

                # --- 1) Full contig sequence ---
                seq_contig = str(fa[tname][0:tlen])
                if strand == "-":
                    seq_contig = str(Seq(seq_contig).reverse_complement())
                out_path_full = args.out_dir_full_seq / f"{roi_id}.fullSeq.fa"
                defline_full_seq = (
                    f">{roi_id}|{tname}"
                    f"|chr:{chr_tok or 'NA'}"
                    f"|full_contig:0-{tlen}"
                    f"|len:{tlen}"
                    f"|end:{end}"
                    f"|strand:{strand}"
                    f"|terminal_ok:{'YES' if terminal_flag else 'NO'}"
                )
                with open(out_path_full, "w") as out:
                    out.write(defline_full_seq + "\n")
                    for i in range(0, len(seq_contig), 80):
                        out.write(seq_contig[i:i + 80] + "\n")
                n_full_ok += 1

                # --- 2) Matching coordinates only (aligned block) ---
                seq_match = str(fa[tname][tstart:tend])
                if strand == "-":
                    seq_match = str(Seq(seq_match).reverse_complement())
                out_path_match = args.out_dir_matching_coord / f"{roi_id}.matchingCoord.fa"
                defline_match = (
                    f">{roi_id}|{tname}"
                    f"|chr:{chr_tok or 'NA'}"
                    f"|align:{tstart}-{tend}"
                    f"|len:{tend - tstart}"
                    f"|end:{end}"
                    f"|strand:{strand}"
                    f"|terminal_ok:{'YES' if terminal_flag else 'NO'}"
                )
                with open(out_path_match, "w") as out:
                    out.write(defline_match + "\n")
                    for i in range(0, len(seq_match), 80):
                        out.write(seq_match[i:i + 80] + "\n")
                n_match_ok += 1

                # --- 3) Aligned+overhang (legacy) ---
                seq_full = str(fa[tname][ext_start:ext_end])
                if strand == "-":
                    seq_full = str(Seq(seq_full).reverse_complement())
                out_path = args.out_dir / f"{roi_id}.besthit.fa"
                defline_full = (
                    f">{roi_id}|{tname}"
                    f"|chr:{chr_tok or 'NA'}"
                    f"|align:{tstart}-{tend}"
                    f"|ext:{ext_start}-{ext_end}"
                    f"|end:{end}"
                    f"|tel_ext_bp:{tel_ext_bp}"
                    f"|strand:{strand}"
                    f"|terminal_ok:{'YES' if terminal_flag else 'NO'}"
                )
                with open(out_path, "w") as out:
                    out.write(defline_full + "\n")
                    for i in range(0, len(seq_full), 80):
                        out.write(seq_full[i:i + 80] + "\n")
                n_ok += 1

                # --- 4) Extension-only FASTA (if available) ---
                if ext_only_start is not None and ext_only_end is not None:
                    seq_ext = str(fa[tname][ext_only_start:ext_only_end])
                    if strand == "-":
                        seq_ext = str(Seq(seq_ext).reverse_complement())

                    out_path_ext = args.out_dir_ext_only / f"{roi_id}.extension_only.fa"
                    defline_ext = (
                        f">{roi_id}|{tname}"
                        f"|chr:{chr_tok or 'NA'}"
                        f"|end:{end}"
                        f"|extension_only:{ext_only_start}-{ext_only_end}"
                        f"|len:{len(seq_ext)}"
                        f"|tel_ext_bp:{tel_ext_bp}"
                        f"|strand:{strand}"
                        f"|terminal_ok:{'YES' if terminal_flag else 'NO'}"
                    )
                    with open(out_path_ext, "w") as out:
                        out.write(defline_ext + "\n")
                        for i in range(0, len(seq_ext), 80):
                            out.write(seq_ext[i:i + 80] + "\n")

                    n_ext_ok += 1

                print(f"[OK] {roi_id} -> fullSeq + matchingCoord + besthit | ext_only={'YES' if ext_only_start is not None else 'NO'}")
            except Exception as e:
                print(f"[FAIL] {roi_id}: {e}", file=sys.stderr)
                n_fail += 1

    print(f"\n[DONE] fullSeq={n_full_ok} matchingCoord={n_match_ok} besthit={n_ok} ext_only={n_ext_ok} skip={n_skip} fail={n_fail}")


if __name__ == "__main__":
    main()
