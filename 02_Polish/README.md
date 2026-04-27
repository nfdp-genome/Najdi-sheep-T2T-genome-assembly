# 02_Polish

Standalone polishing and gap-improvement steps after initial assembly/scaffolding.

## Files in this folder

- `step_polish_racon.sh` + `step_polish_racon.config.yaml`
- `step_polish_medaka.sh` + `step_polish_medaka.config.yaml`
- `step_gap_tgsgapcloser.sh` + `step_gap_tgsgapcloser.config.yaml`
- `step_gap_ntlink.sh` + `step_gap_ntlink.config.yaml`
- `step_funannotate_clean.sh` + `step_funannotate_clean.config.yaml`

## What each step does

- `step_polish_racon.sh`: ONT-based iterative racon polishing.
- `step_polish_medaka.sh`: ONT-based medaka consensus polishing.
- `step_gap_tgsgapcloser.sh`: gap closing with TGS-GapCloser.
- `step_gap_ntlink.sh`: ntLink gap fill/scaffold improvement.
- `step_funannotate_clean.sh`: contig cleanup/dedup style filtering.

## How to run

Edit the matching YAML first, then submit with `sbatch`.

```bash
sbatch step_polish_racon.sh --config step_polish_racon.config.yaml
sbatch step_polish_medaka.sh --config step_polish_medaka.config.yaml
sbatch step_gap_tgsgapcloser.sh --config step_gap_tgsgapcloser.config.yaml
sbatch step_gap_ntlink.sh --config step_gap_ntlink.config.yaml
sbatch step_funannotate_clean.sh --config step_funannotate_clean.config.yaml
```

## Notes

- `MEDAKA_MODEL` is configured in YAML (`step_polish_medaka.config.yaml`).
- `POLISHER` and `TAG` are names/identifiers, not paths.
- `NTLINK_BIN` and `FUNANNOTATE_RUNNER` are prefilled in templates; adjust only if environment changes.
