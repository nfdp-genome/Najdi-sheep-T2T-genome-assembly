# Command inventory (provenance for `pipeline/steps/*.sh`)

Frozen tool invocations are copied from legacy `[scripts/](../scripts)`. Path arguments are supplied via environment variables set by the driver from YAML.

## hifiasm (telomere + ONT ultra-long)

**Source:** `[scripts/hifiasm_assembly_telomere_flag_noDual-scaf_v1.0.sh](../scripts/hifiasm_assembly_telomere_flag_noDual-scaf_v1.0.sh)`

```bash
hifiasm \
  -o "${PREFIX}.asm" \
  -t "${SLURM_CPUS_PER_TASK}" \
  --telo-m TTAGGG \
  --ul "${ONT}" \
  "${HIFI}" \
  2> "${PREFIX}.hifiasm.log"
```

**Env mapping:** `HIFIASM_PREFIX` → `-o ${HIFIASM_PREFIX}.asm`, `THREADS`, `HIFIASM_ONT`, `HIFIASM_HIFI`, `HIFIASM_LOG`.

## GFA → FASTA (primary contigs)

**Source:** `[scripts/hifiasm_assembly_quality_v1.3.sh](../scripts/hifiasm_assembly_quality_v1.3.sh)` — `gfa_to_fa` uses `awk` + `seqkit replace`; pipeline step uses **awk only** for portability (documented deviation: no `seqkit replace`).

```bash
awk '$1=="S"{print ">"$2"\n"$3}' "$in_gfa" > "$out_fa"
```

## YaHS

**Source:** `[scripts/hiC_Yahs_scaffolding_v1.0.sh](../scripts/hiC_Yahs_scaffolding_v1.0.sh)`

```bash
yahs "${ASM_BASENAME}" "${BAM}" -o "${PREFIX}"
# then copy ${PREFIX}_scaffolds_final.fa → CLEAN_NAME
```

## SALSA (`run_pipeline.py`)

**Source:** `[scripts/hiC_Salsa_scaffolding_v1.1.sh](../scripts/hiC_Salsa_scaffolding_v1.1.sh)`

```bash
"${SALSA_RUN}" \
    -a "${ASM_FASTA_BASENAME}" \
    -l "${LENGTHS_FILE}" \
    -b "${BED_SORTED}" \
    -e DNASE \
    -i 3 \
    -m yes \
    -o "${OUTDIR}"
```

## Racon (ONT polish)

**Source:** `[scripts/polish_racon_hifiasm_hic_with_ont_V1.0.sh](../scripts/polish_racon_hifiasm_hic_with_ont_V1.0.sh)`

```bash
ROUNDS=2  # fixed in legacy
minimap2 -x map-ont -t "${THREADS}" "${CURRENT_ASM}" "ont.fastq.gz" > "${PAF_OUT}"
racon -t "${THREADS}" "ont.fastq.gz" "${PAF_OUT}" "${CURRENT_ASM}" > "${ASM_OUT}"
```

## Medaka

**Source:** `[scripts/polish_medaka_hifiasm_hic_with_ont_V1.0.sh](../scripts/polish_medaka_hifiasm_hic_with_ont_V1.0.sh)`

```bash
minimap2 -ax map-ont -t "${THREADS}" "${REF_FA}" "ont.fastq.gz" \
    | samtools sort -@ "${THREADS}" -o align.bam
samtools index align.bam
medaka_consensus \
    -i "ont.fastq.gz" \
    -d "${REF_FA}" \
    -o "${OUT_DIR}" \
    -t "${THREADS}" \
    -m "${MEDAKA_MODEL}"
```

## ntLink

**Source:** `[scripts/gapfilling_hifiasm_hic_ont_ntlink_v1.0.sh](../scripts/gapfilling_hifiasm_hic_ont_ntlink_v1.0.sh)`

```bash
"${NTLINK_BIN}" scaffold gap_fill \
    target="${TARGET_FA}" \
    reads="${READS_FAQ}" \
    prefix="${PREFIX}" \
    t="${THREADS}" \
    k=32 \
    w=250 \
    > "${PREFIX}.ntlink.log" 2> "${PREFIX}.ntlink.err"
```

## TGS-GapCloser

**Source:** `[scripts/gapfilling_hifiasm_hic_ont_tgsgapcloser_v1.1.sh](../scripts/gapfilling_hifiasm_hic_ont_tgsgapcloser_v1.1.sh)`

```bash
"${TGSGC_BIN}" \
  --scaff  "${TAG}.fasta" \
  --reads  "${ONT_FASTA}" \
  --output "${OUT_PREFIX}" \
  --racon  "${RACON_BIN}" \
  --thread "${THREADS}" \
  > "${OUT_PREFIX}.log" 2> "${OUT_PREFIX}.err"
```

## funannotate clean (V3)

**Source:** `[scripts/funannotate_clean_V3.sh](../scripts/funannotate_clean_V3.sh)`

```bash
micromamba run -n funannotate_clean_env \
  funannotate clean \
  -i "${GENOME}" \
  -o "${OUT_FASTA}" \
  --minlen 15000000 \
  --pident 95 \
  --cov 95 \
  --exhaustive
```

**Env:** `FUNANNOTATE_RUNNER` default `micromamba run -n funannotate_clean_env` (override in config via wrapper env).

## BUSCO + QUAST (QC wave)

**Source:** `[scripts/run_quast_busco_v1.1.sh](../scripts/run_quast_busco_v1.1.sh)`

```bash
busco \
  -i "${FA}" \
  -o "${OUTNAME1}" \
  -m genome \
  -l "${L1_PATH}" \
  -c "${CPUS}" \
  --out_path "${BUSCO_OUT_ROOT}"
```

```bash
quast --large -t "${CPUS}" -o "$qdir" "$fa" &> "$qlog"
```

## Dotplot (nucmer + DotPrep)

**Source:** `[scripts/dotplot_references_LR.T1.5.sh](../scripts/dotplot_references_LR.T1.5.sh)`

```bash
nucmer \
  --mum \
  -l 100 \
  -c 500 \
  -t "$SLURM_CPUS_PER_TASK" \
  "$reference" \
  "$sample_file" \
  -p "$outprefix"

python3 "$DOTPREP" \
  --delta "${outprefix}.delta" \
  --out   "$outprefix"
```

## Multi-sample QUAST

**Source:** `[scripts/quast_all_in_one.sh](../scripts/quast_all_in_one.sh)`

```bash
quast.py \
  "${GENOMES[@]}" \
  --large \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --output-dir "${OUT_DIR}"
```

Driver passes explicit file list + optional `--labels` built from assembly IDs (same argv pattern as `quast_all_in_one`; labels added only as extra flags mirroring QUAST CLI, file list from driver).