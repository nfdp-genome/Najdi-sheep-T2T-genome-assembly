#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# 15_attach_telomere_extensions.py
#
# Final step: for each tagged chromosome FASTA (produced by script 3), strips the
# @@TELROI|...@@ placeholders to recover the clean chromosome sequence, then prepends the 5'
# extension-only FASTA (if that end was on the rescue list) and appends the 3' extension-only
# FASTA (if that end was on the rescue list). Writes one *_TELEXTENDED.fasta per chromosome
# plus an attach_summary.tsv logging what was attached and the new chromosome lengths.
# ------------------------------------------------------------------------------
"""
attach_telomere_extensions_v1.0.py

Attach rescued telomere extension-only sequences to tagged chromosomes.

Inputs (defaults are absolute paths for the Najdi T2T project):
- Rescue list TSV: telomeres_rescue_list.tsv
    columns: SampleID, Chromosome, end
- Tagged chromosomes: one FASTA per chromosome with TELROI tags embedded
    e.g. 1268_NLFDP1268_chr01_TELROI_TAGGED.fasta
- Extension-only FASTAs: telomeres_best_hits_extension_only/
    e.g. 1268_NLFDP1268_chr01_3p_TELROI.extension_only.fa

For each tagged chromosome:
- Strip all @@TELROI|...@@ tags from the sequence to get a clean base.
- If (sample, chr, 5p) is in the rescue list and the corresponding
  {prefix}_5p_TELROI.extension_only.fa exists, prepend that extension.
- If (sample, chr, 3p) is in the rescue list and the corresponding
  {prefix}_3p_TELROI.extension_only.fa exists, append that extension.

Outputs:
- <base>/telomeres_extended/fasta/<prefix>_TELEXTENDED.fasta
    where prefix is, e.g., 1268_NLFDP1268_chr01
- <base>/telomeres_extended/logs/attach_summary.tsv
    one row per chromosome with flags and lengths.

Notes:
- Every tagged chromosome is written once. Chromosomes with no rescued
  ends are written unchanged (tags removed, no extension).
- Chromosomes with only one rescued end get an extension only on that end.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ---------- ABSOLUTE DEFAULT BASE (edit here only if project moves) ----------
DEFAULT_BASE = Path(
    "/ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/data/processed/"
    "hifiasm_assembly/chromosomal_reconstruction/anchor_based_gap_filling"
)
# ---------------------------------------------------------------------------

TAG_START = b"@@TELROI|"
TAG_END = b"@@"
DEFAULT_WRAP = 60


def read_fasta_header_and_seq_bytes(fa_path: Path) -> Tuple[str, bytes]:
    """
    Read FASTA: return (header_without_>, sequence_bytes_without_newlines).
    Reads only the first record. Works for wrapped or one-line FASTA.
    """
    header = None
    seq_chunks: List[bytes] = []
    with fa_path.open("rb") as f:
        for line in f:
            if line.startswith(b">"):
                if header is None:
                    header = line[1:].strip().decode("utf-8", errors="replace")
                else:
                    break
            else:
                if header is None:
                    continue
                line = line.strip()
                if line:
                    seq_chunks.append(line)
    if header is None:
        raise ValueError(f"No FASTA header found in {fa_path}")
    return header, b"".join(seq_chunks)


def read_single_fasta_sequence_bytes(fa_path: Path) -> bytes:
    """Read first sequence record from a FASTA; return sequence bytes without newlines."""
    _, seq = read_fasta_header_and_seq_bytes(fa_path)
    return seq


def find_all_telroi_tags(seq: bytes) -> List[Tuple[int, int, bytes]]:
    """
    Find all @@TELROI|...@@ tags in sequence (byte positions).
    Returns list of tuples: (start, end_excl, tag_bytes)
    """
    hits: List[Tuple[int, int, bytes]] = []
    i = 0
    while True:
        j = seq.find(TAG_START, i)
        if j == -1:
            break
        k = seq.find(TAG_END, j + len(TAG_START))
        if k == -1:
            raise ValueError("Found @@TELROI| but no closing @@")
        k_excl = k + len(TAG_END)
        tag_bytes = seq[j:k_excl]
        hits.append((j, k_excl, tag_bytes))
        i = k_excl
    return hits


def strip_telroi_tags(seq: bytes) -> bytes:
    """Remove all @@TELROI|...@@ tags from sequence bytes."""
    hits = find_all_telroi_tags(seq)
    if not hits:
        return seq
    chunks: List[bytes] = []
    prev = 0
    for start, end_excl, _ in hits:
        if start > prev:
            chunks.append(seq[prev:start])
        prev = end_excl
    if prev < len(seq):
        chunks.append(seq[prev:])
    return b"".join(chunks)


def write_wrapped_fasta(out_path: Path, header: str, seq: bytes, wrap: int) -> int:
    """Write a single FASTA record with wrapping; return sequence length."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_len = len(seq)
    with out_path.open("wb") as out:
        out.write(b">" + header.encode("utf-8") + b"\n")
        line_len = 0
        mv = memoryview(seq)
        idx = 0
        while idx < len(mv):
            space = wrap - line_len
            chunk = mv[idx : idx + space]
            out.write(chunk)
            idx += len(chunk)
            line_len += len(chunk)
            if line_len == wrap:
                out.write(b"\n")
                line_len = 0
        if line_len != 0:
            out.write(b"\n")
    return total_len


def infer_prefix_from_tagged_filename(tagged_path: Path) -> str:
    """
    From: 1268_NLFDP1268_chr01_TELROI_TAGGED.fasta -> 1268_NLFDP1268_chr01
    Supports: 1271_NLFDP1271_chr08_Hu-F_RagTag_TELROI_TAGGED.fasta
    """
    suffix = "_TELROI_TAGGED.fasta"
    name = tagged_path.name
    if not name.endswith(suffix):
        raise ValueError(f"Unexpected tagged filename format: {name}")
    return name[: -len(suffix)]


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
    m = re.search(r"(?:^|[_\W])chr(\d+)(?:[_\W]|$)", s)
    if m:
        n = int(m.group(1))
        return f"chr{n:02d}"
    if re.search(r"(?:^|[_\W])chr([XY])(?:[_\W]|$)", s, re.IGNORECASE):
        return re.search(r"chr([XY])", s, re.IGNORECASE).group(0)
    return None


def load_rescue_list(path: Path) -> Set[Tuple[str, str, str]]:
    """
    Load rescue list rows as (sample, chrNN, end).
    Expected columns: SampleID, Chromosome, end
    """
    wanted: Set[Tuple[str, str, str]] = set()
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Rescue-list TSV has no header")
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


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Attach telomere extension-only sequences to tagged chromosomes."
    )

    ap.add_argument("--base", default=str(DEFAULT_BASE), help="Project base directory (absolute default set).")
    ap.add_argument("--rescue-list", default=None, help="Default: <base>/telomeres_rescue_list.tsv")
    ap.add_argument("--tagged-dir", default=None, help="Default: <base>/tagged_chromosomes")
    ap.add_argument(
        "--extension-dir",
        default=None,
        help="Default: <base>/telomeres_best_hits_extension_only",
    )
    ap.add_argument(
        "--out-dir",
        default=None,
        help="Default: <base>/telomeres_extended",
    )
    ap.add_argument(
        "--wrap",
        type=int,
        default=DEFAULT_WRAP,
        help=f"FASTA wrap length (default {DEFAULT_WRAP})",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and log lengths only; do not write FASTA files.",
    )

    args = ap.parse_args()

    base = Path(args.base).resolve()
    rescue_list_path = (
        Path(args.rescue_list).resolve()
        if args.rescue_list
        else (base / "telomeres_rescue_list.tsv").resolve()
    )
    tagged_dir = (
        Path(args.tagged_dir).resolve()
        if args.tagged_dir
        else (base / "tagged_chromosomes").resolve()
    )
    extension_dir = (
        Path(args.extension_dir).resolve()
        if args.extension_dir
        else (base / "telomeres_best_hits_extension_only").resolve()
    )
    out_dir = (
        Path(args.out_dir).resolve()
        if args.out_dir
        else (base / "telomeres_extended").resolve()
    )

    out_fasta_dir = out_dir / "fasta"
    out_logs_dir = out_dir / "logs"
    out_logs_dir.mkdir(parents=True, exist_ok=True)

    attach_summary = out_logs_dir / "attach_summary.tsv"

    rescue_set = load_rescue_list(rescue_list_path)

    tagged_files = sorted(tagged_dir.glob("*_TELROI_TAGGED.fasta"))
    if not tagged_files:
        raise SystemExit(f"No tagged chromosomes found in: {tagged_dir}")

    with attach_summary.open("w", encoding="utf-8") as s:
        s.write(
            "prefix\tsample\tchr\tattached_5p\tattached_3p\t"
            "len_5p_ext\tlen_3p_ext\tbase_len\tnew_len\tout_fasta\n"
        )

    for tagged_path in tagged_files:
        prefix = infer_prefix_from_tagged_filename(tagged_path)
        sample = prefix.split("_", 1)[0]
        chr_tok = extract_chr_token(prefix)

        if not chr_tok:
            # Cannot map to rescue list; treat as no rescue.
            has_5p = False
            has_3p = False
        else:
            has_5p = (sample, chr_tok, "5p") in rescue_set
            has_3p = (sample, chr_tok, "3p") in rescue_set

        ext5_path = extension_dir / f"{prefix}_5p_TELROI.extension_only.fa"
        ext3_path = extension_dir / f"{prefix}_3p_TELROI.extension_only.fa"

        use_5p = has_5p and ext5_path.exists()
        use_3p = has_3p and ext3_path.exists()

        ext5_seq = b""
        ext3_seq = b""
        len_5p_ext = 0
        len_3p_ext = 0

        if use_5p:
            ext5_seq = read_single_fasta_sequence_bytes(ext5_path)
            len_5p_ext = len(ext5_seq)
        if use_3p:
            ext3_seq = read_single_fasta_sequence_bytes(ext3_path)
            len_3p_ext = len(ext3_seq)

        header, seq = read_fasta_header_and_seq_bytes(tagged_path)
        base_seq = strip_telroi_tags(seq)
        base_len = len(base_seq)

        final_seq = ext5_seq + base_seq + ext3_seq
        new_len = len(final_seq)

        attached_5p = "Y" if use_5p else "N"
        attached_3p = "Y" if use_3p else "N"

        if not args.dry_run:
            out_fa = out_fasta_dir / f"{prefix}_TELEXTENDED.fasta"
            new_header = f"{header} | TELEXTENDED | 5p_ext:{len_5p_ext} | 3p_ext:{len_3p_ext}"
            write_wrapped_fasta(out_fa, new_header, final_seq, args.wrap)
            out_path_str = str(out_fa)
        else:
            out_path_str = "DRY_RUN"

        with attach_summary.open("a", encoding="utf-8") as s:
            s.write(
                f"{prefix}\t{sample}\t{chr_tok or 'NA'}\t"
                f"{attached_5p}\t{attached_3p}\t"
                f"{len_5p_ext}\t{len_3p_ext}\t"
                f"{base_len}\t{new_len}\t{out_path_str}\n"
            )


if __name__ == "__main__":
    main()
