#!/bin/bash
#
# SLURM directives:
#SBATCH --job-name=extract_eqtlgen_mucin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

cd /data/gen1/reference/eqtlgen
> /scratch/gen1/yz735/updated_coloc/eqtlgen_mucin.txt
zcat 2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | awk '$3 == "11" && $4 >= 1012240 && $4 <= 1295601'  >> /scratch/gen1/yz735/updated_coloc/eqtlgen_mucin.txt

hostname
date
