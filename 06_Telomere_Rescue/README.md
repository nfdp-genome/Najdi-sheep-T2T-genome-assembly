# 06_Telomere_Rescue

Custom anchor-based workflow used to detect, classify and (where possible)
extend telomeric repeat arrays at chromosome termini of the Najdi male
(NLFDP1268) and female (NLFDP1271) assemblies.

The workflow is applied **after** Stage 3 (chromosome reconstruction and
internal-gap resolution) and produces the final telomere-rescued
chromosome FASTAs from which the shipped v1.0 assemblies are derived.

## Workflow overview

The rescue is anchor-based: at each chromosome end where the assembled
sequence does not itself carry a canonical telomeric array, the terminal
50 kb (or 200 kb, for ends with no detectable telomere) is used as a
query to search the per-sample pool of prior intermediate assemblies and
unplaced contigs for a candidate that (i) aligns to the terminal region
and (ii) extends outward past it with a telomere-bearing overhang. Only
the outward telomeric overhang from a validated candidate is attached to
the anchored chromosome end.

## Main workflow — run in numerical order

### Stage 1 — build and validate the TELROI queries

| # | Script | What it does |
|---|--------|--------------|
| 01 | `01_telomere_anchor_scan.py` | Scan every chromosome end for TTAGGG/CCCTAA in 1 kb windows across the terminal 50 kb; record first inward window with ≥10 motifs; define 20 kb internal anchor per end. |
| 02 | `02_build_telroi_decision_table.py` | Classify each end (Class A ≥150, Class B 10–149, Class C <10 motifs/kb); keep only B and C; compute 1-based ROI coordinates (50 kb for B, 200 kb for C), actual ROI length, and a unique `TELROI_tag`. No sequence is modified. |
| 03 | `03_tag_chromosome_fastas_with_telroi.py` | Insert `@@TELROI\|sample\|chrom\|end\|class\|actual_roi_len\|roi:start-end@@` placeholders on a working copy of each selected chromosome FASTA (5' before the sequence, 3' after). |
| 04 | `04_generate_expected_placeholder_list.py` | QC helper: enumerate the placeholders that must appear given the decision table. |
| 05 | `05_validate_tagged_chromosome_fastas.sh` | QC pass 1: expected-vs-observed placeholder lists (must be identical). |
| 06 | `06_validate_telroi_tag_placement.py` | QC pass 2: row-level sanity table checking each placeholder's position, coordinates, ROI length, chromosome-ID match, and chromosome bounds. |
| 07 | `07_extract_telroi_query_sequences.py` | Strip **all** placeholders on a working copy, then slice by 1-based inclusive TELROI coordinates and write one query FASTA per end. |
| 08 | `08_validate_extracted_telroi_fastas.sh` | QC of extracted TELROI FASTAs: file count, header format, absence of internal tags, DNA-only content. |

### Stage 2 — search the assembly pool and rank candidates

| # | Script | What it does |
|---|--------|--------------|
| 09 | `09_uniquify_merged_assembly_pool_headers.sh` | Rewrite the merged per-sample assembly-pool FASTA deflines with a sequential `ASMxxxx___` prefix so all sequence IDs are unique. |
| 10 | `10_map_telroi_queries_to_pool.sh` | Wrapper: submits `11_mapping_assemblies.V1.12.sh` as a SLURM array over all TELROI queries per sample. |
| 11 | `11_mapping_assemblies.V1.12.sh` | The actual per-query SLURM job: runs `minimap2 -x asm5` for each TELROI query vs its sample's merged assembly pool, producing a PAF (and optionally a sorted BAM for IGV inspection). |
| 12 | `12_select_two_best_hits_per_roi.py` | Reads the aggregated all-hits TSV and picks the top 1 or 2 candidates per TELROI. Ranks by: distance from contig end (asc), alignment length (desc), sequence identity (desc). Applies terminality, minimum alignment length, and min/max overhang filters. |
| 13 | `13_filter_best_hits_by_rescue_list.py` | Filters the ranked hits down to only the (sample, chromosome, end) rows on the rescue list, i.e. the ends actually selected for rescue after manual review. |

### Stage 3 — extract extensions and attach to chromosomes

| # | Script | What it does |
|---|--------|--------------|
| 14 | `14_extract_telomere_best_hits_plus_overhang.py` | Extracts four FASTA products per surviving best-hit from the candidate contig: (a) full contig, (b) aligned block only, (c) aligned block + telomere-side overhang (legacy), and (d) the extension-only segment — the telomeric overhang beyond the aligned block. Product (d) is the input to script 15. |
| 15 | `15_attach_telomere_extensions.py` | For each tagged chromosome FASTA from script 3: strips the placeholders to recover the clean chromosome, then prepends the 5' extension-only FASTA (if that end is on the rescue list) and appends the 3' extension-only FASTA (if that end is on the rescue list). Writes one `*_TELEXTENDED.fasta` per chromosome plus an `attach_summary.tsv` logging what was attached and the new chromosome lengths. |

## Companion scripts (not part of the main run order)

### `alt_scan_and_classify/`

Standalone equivalents of scripts 01 and 02. Kept for reference; not part of the main run order.

| Script | What it does |
|--------|--------------|
| `telomere_scan_chromosomes.py` | Standalone chromosome-FASTA scan for TTAGGG/CCCTAA repeats (functionally equivalent to `01`). |
| `telomere_classify_and_build_telroi.py` | Standalone A/B/C classification and TELROI decision table (functionally equivalent to `02`). |

### `validation/`

Independent ONT-based validation of the shipped v1.0 genomes.

| Script | What it does |
|--------|--------------|
| `validate_v1_by_ONT_mapping.sh` | Maps v1.0 vs ONT reads with `minimap2 -x map-ont`, produces PAF, adds a header row, and prepares the filtered manifest of terminal, extension-bearing reads. |
