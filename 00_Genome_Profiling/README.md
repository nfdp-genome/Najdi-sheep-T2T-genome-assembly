# 00_Genome_Profiling

Pre-assembly k-mer genome profiling (genome size, heterozygosity, ploidy) with KMC + GenomeScope2.

- `00_genomescope_array.sh` — Slurm array job: KMC k=21 counts per HiFi sample (1271, 1268), builds the k-mer histogram, then runs GenomeScope2.
