#!/bin/bash
#SBATCH --job-name=najdi_driver
#SBATCH --output=/ibex/user/x_altuway/najdi_pipeline/1268/logs/driver.out
#SBATCH --error=/ibex/user/x_altuway/najdi_pipeline/1268/logs/driver.err
#SBATCH --time=336:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --partition=batch

cd /ibex/project/c2293/najdi_t2t_project/WP1_genome_assembly_qc/pipeline/Najdi-sheep-T2T-genome-assembly/01-02_automated_pipeline

python driver.py --config my_run_1268.yaml
