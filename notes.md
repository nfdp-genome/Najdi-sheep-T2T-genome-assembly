<<<<<<< HEAD

# Najdi T2T — config-driven assembly pipeline

YAML-driven DAG: HiFiasm → GFA→FASTA → Hi-C scaffold (YaHS / SALSA) → polish (Racon / Medaka) → gap fill (ntLink / TGS-GapCloser) → optional funannotate clean → QUAST. Supports **per-phase** continuation from HiFiasm primary / hap1 / hap2 via `hifiasm.contig_phases`.

This folder is a **self-contained export**: copy or `git init` here and push to GitHub as its own repository (no `tools/`, `results/`, or conda installs).

## Quick start

```bash
cd najdi_t2t_pipeline
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
cp config.example.yaml my_run.yaml
# Edit my_run.yaml: project.workdir, inputs.*, execution, stages

./.venv/bin/python driver.py --config my_run.yaml --dry-run
./.venv/bin/python driver.py --config my_run.yaml
```

- **SLURM**: `execution.mode: slurm` — tune `execution.slurm.modules_by_step` for your cluster.
- **Local**: `execution.mode: local` — optional Docker; see [docker/README.md](docker/README.md).

## Layout


| Path                        | Purpose                                            |
| --------------------------- | -------------------------------------------------- |
| `config.example.yaml`       | Template — copy to `my_run.yaml` (gitignored)      |
| `driver.py`                 | Config parser, job DAG, SLURM/local runner, resume |
| `steps/step_*.sh`           | One wrapper per tool                               |
| `docs/command_inventory.md` | Command provenance                                 |


## Push this folder to GitHub

From **inside** this directory:

```bash
cd najdi_t2t_pipeline
git init
git add -A
git status   # confirm no secrets / no my_run.yaml
git commit -m "Initial import: Najdi sheep T2T assembly pipeline driver and step scripts"
git branch -M main
git remote add origin https://github.com/nfdp-genome/Najdi-sheep-T2T-genome-assembly.git
git push -u origin main
```

If `git remote add` fails because `origin` already exists: `git remote set-url origin https://github.com/nfdp-genome/Najdi-sheep-T2T-genome-assembly.git`

Keep private paths and sample IDs only in **local** YAML files that you do not commit (`my_run.yaml` is listed in `.gitignore`).

## Observability

# Per-job logs, driver log, `failures.tsv`, `manifest.jsonl`, and `run_summary.json` under `workdir/logs/` (see `config.example.yaml` → `logging`).

# Najdi-sheep-T2T-genome-assembly

Public repository

> > > > > > > 6b2ea2ffcbae0f5e2f887e5fcfbe9b22c23c2dc6

