# 00_Genome_Profiling

Pre-assembly k-mer genome profiling (genome size, heterozygosity, ploidy) with KMC + GenomeScope2.

- `00_genomescope_array.sh` — original Slurm array job: KMC (k=21) counts per HiFi sample, builds the k-mer histogram, then runs GenomeScope2.
- `00_genomescope_array.config.sh` — same workflow, but path-free: all inputs and parameters are read from `config.yaml`, so nothing site-specific is baked into the script.
- `config.yaml` — sample list (`<id> <hifi_path>`) plus GenomeScope2 path and KMC/GenomeScope parameters. Fill in the placeholder paths and keep `#SBATCH --array=1-N` in sync with the number of samples.

Run: `sbatch 00_genomescope_array.config.sh` (or `CONFIG=/path/to/config.yaml sbatch 00_genomescope_array.config.sh`).
