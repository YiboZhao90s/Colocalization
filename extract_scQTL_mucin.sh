#!/bin/bash
#
# SLURM directives:
#SBATCH --job-name=extract_scQTL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

cd /data/gen1/yz735/scQTL_limix

> /scratch/gen1/yz735/updated_coloc/scQTL_limix_mucin.txt

for file in *.txt; do
    echo "Processing $file..." >&2
    zcat "$file" | awk -v fname="$file" '$15 == "11" && $2 >= 1012240 && $2 <= 1295601 { print $0 "\t" fname }' >> /scratch/gen1/yz735/updated_coloc/scQTL_limix.txt
done

hostname
date
