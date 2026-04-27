# 01_Assembly

Standalone assembly steps for generating draft assemblies.

## Files in this directory

- `step_hifiasm.sh` + `step_hifiasm.config.yaml`
- `step_gfa_primary.sh` + `step_gfa_primary.config.yaml`
- `step_hic_map.sh` + `step_hic_map.config.yaml`
- `step_salsa.sh` + `step_salsa.config.yaml`
- `step_yahs.sh` + `step_yahs.config.yaml`

## What each step does

- `step_hifiasm.sh`: HiFi/ONT de novo assembly.
- `step_gfa_primary.sh`: converts primary GFA to FASTA.
- `step_hic_map.sh`: maps Hi-C reads to assembly and produces sorted/indexed BAM.
- `step_salsa.sh`: SALSA scaffolding using Hi-C BAM.
- `step_yahs.sh`: YaHS scaffolding using Hi-C BAM.

## How to run

Edit the matching YAML first, then submit with `sbatch`.

```bash
sbatch step_hifiasm.sh --config step_hifiasm.config.yaml
sbatch step_gfa_primary.sh --config step_gfa_primary.config.yaml
sbatch step_hic_map.sh --config step_hic_map.config.yaml
sbatch step_salsa.sh --config step_salsa.config.yaml
sbatch step_yahs.sh --config step_yahs.config.yaml
```

## Notes

- Each script expects `--config <file>.config.yaml`.
- Main step log path comes from `log_file` in YAML.
- Slurm stdout/stderr are written by each script's `#SBATCH` settings.
- `TAG` values in configs are identifiers (not filesystem paths).
