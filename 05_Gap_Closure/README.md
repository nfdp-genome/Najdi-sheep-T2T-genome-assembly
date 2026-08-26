# 05_Gap_Closure

Custom anchor-based workflow used to close internal chromosome gaps in the
Najdi male (NLFDP1268) and female (NLFDP1271) assemblies after chromosome
reconstruction (Stage 3 in the manuscript). The output of this folder is a
gap-free per-sample whole-genome FASTA that feeds into the telomere-rescue
workflow in `06_Telomere_Rescue/`.

## Workflow overview

The gap closure is anchor-based: for each declared internal gap on a
reconstructed chromosome, a 50 kb region on either side of the gap
(the "flank + gap + flank" ROI) is used as a query to search the
per-sample pool of prior intermediate assemblies and unplaced contigs
for a candidate contig that spans the gap plus both flanks in a
non-repetitive context. The selected candidate segment is patched into
the gap and the reconstructed chromosome is spliced back into the
whole-genome FASTA.

Four gap sites were resolved this way for the shipped v1.0 genomes:
- 1268 chr02 — a 338 bp true gap
- 1268 chr10 — a multi-contig chromosome joined earlier by RagTag
- 1271 chr08 — a multi-contig chromosome joined earlier by RagTag
- 1271 chr17 — a multi-contig chromosome joined earlier by RagTag

## Run order

### Stage 1 — declare inputs and verify gap coordinates

| # | Script | What it does |
|---|--------|--------------|
| 01 | `01_phase0_declare_inputs.py` | Declare the per-sample input paths (scaffolded chromosome FASTA + gap coordinates table). |
| 02 | `02_phase05_gap_coordinate_preflight.py` | Pre-flight sanity checks on the gap coordinates table (bounds, columns, no overlapping gaps). |
| 03 | `03_verify_gap_coordinates_on_fastas.py` | Coordinate-only verification on the original chromosome FASTAs (no modification) — confirms every declared gap position is actually a run of N's of the expected length. |

### Stage 2 — build gap ROI queries

| # | Script | What it does |
|---|--------|--------------|
| 04 | `04_phase1_tag_gap_flanks.py` | Insert `@@LEFT@@` / `@@RIGHT@@` boundary placeholders on a working copy of each chromosome FASTA so the gap ROI (flank + gap + flank) can be sliced unambiguously. |
| 05 | `05_tagging_sanity_check.py` | Sanity check on the tagged FASTAs: expected `@@LEFT@@` / `@@RIGHT@@` pairs are present, in order, and bracket a stretch of N's of the expected length. |
| 06 | `06_phase2_extract_gap_roi.py` | Strip all placeholders on a working copy, then slice by 1-based inclusive coordinates to write one ROI query FASTA per gap. |

### Stage 3 — search the assembly pool and rank candidates

The mapping step itself reuses the same script and per-sample assembly
pool as the telomere-rescue workflow — see
`06_Telomere_Rescue/11_mapping_assemblies.V1.12.sh`. This folder picks
up from the aggregated PAF output.

| # | Script | What it does |
|---|--------|--------------|
| 07 | `07_filter_paf_and_rank_gap_candidates.py` | Reads the aggregated PAF from mapping gap ROIs against the merged assembly pool, ranks candidates that span the full gap plus both flanking regions, and writes the per-gap best-hit table used by the patching steps below. |

### Stage 4 — patch selected candidates into chromosomes

| # | Script | What it does |
|---|--------|--------------|
| 08 | `08_patch_1271_chr08_and_chr17.sh` | Patch the selected gap-filling contigs into chromosomes 08 and 17 of sample 1271. |
| 09 | `09_patch_1268_chr02.sh` | Patch the selected gap-filling contig into chromosome 02 of sample 1268 (338 bp true gap). |
| 10 | `10_patch_1268_chr10.sh` | Patch the selected gap-filling contig into chromosome 10 of sample 1268. |

### Stage 5 — assemble the gap-free per-sample genome

| # | Script | What it does |
|---|--------|--------------|
| 11a | `11a_quality_check_after_patch.sh` | Quality check on the four patched chromosomes: total length, gap count, N's per 100 kbp. |
| 11b | `11b_build_gapfree_genome.sh` | Assemble the final gap-free per-sample chromosome FASTAs into one whole-genome FASTA per sample. Output feeds into `06_Telomere_Rescue/`. |
