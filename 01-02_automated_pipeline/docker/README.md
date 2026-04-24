# Local Docker (optional)

Docker is **only** intended when `execution.mode: local` and `execution.local.use_docker: true`. On IBEX, use `mode: slurm` and cluster modules instead.

1. Build or pull images that contain the same tools as your `pipeline/steps` scripts (hifiasm, quast, busco, yahs, minimap2/racon, medaka, ntLink, tgsgapcloser, mummer, funannotate, etc.). The example image keys in `config.example.yaml` are placeholders — replace with images you maintain.
2. Set `execution.local.bind_mounts` so the pipeline `workdir` and any read-only reference paths are visible inside the container at the **same absolute paths** as on the host (the driver passes host paths in environment variables).
3. Run: `python pipeline/driver.py --config my_run.yaml` with `use_docker: true`.

There is no requirement that local Docker versions match IBEX modules; use Docker for laptop/small-server smoke tests only if you accept tool version drift.