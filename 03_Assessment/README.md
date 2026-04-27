# 03_Assessment

Standalone assembly quality assessment steps (BUSCO/QUAST/dotplot).

## Files in this folder

- `step_qc_wave.sh` + `step_qc_wave.config.yaml`
- `step_quast_multi.sh` + `step_quast_multi.config.yaml`
- `step_dotplot.sh` + `step_dotplot.config.yaml`
- `step_qc_bundle.sh` + `step_qc_bundle.config.yaml`

## What each step does

- `step_qc_wave.sh`: BUSCO (2 lineages) + QUAST for one assembly.
- `step_quast_multi.sh`: QUAST on multiple assemblies (list-driven).
- `step_dotplot.sh`: nucmer + DotPrep dotplot generation.
- `step_qc_bundle.sh`: convenience wrapper (`step_qc_wave` + optional dotplots).

## How to run

Edit the matching YAML first, then submit with `sbatch`.

```bash
sbatch step_qc_wave.sh --config step_qc_wave.config.yaml
sbatch step_quast_multi.sh --config step_quast_multi.config.yaml
sbatch step_dotplot.sh --config step_dotplot.config.yaml
sbatch step_qc_bundle.sh --config step_qc_bundle.config.yaml
```

## Notes

- `step_qc_bundle.sh` is optional but recommended for one-command QC.
- BUSCO lineage paths and dotplot reference defaults are prefilled in templates.
- `LABEL` is an identifier, not a filesystem path.
