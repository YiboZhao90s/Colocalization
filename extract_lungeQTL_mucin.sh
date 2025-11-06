#!/bin/bash
#
# SLURM directives:
#SBATCH --job-name=extract_lungeQTL_mucin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

> /scratch/gen1/yz735/updated_coloc/lungeQTL_mucin.txt

for file in *.txt.gz; do
    zcat /data/gen1/reference/lung_eQTL/METAANALYSIS_Laval_UBC_Groningen_chr11_formatted.txt.gz | awk -F '$3 == 11 && $4 >= 1012240  && $4 <= 1274371' >> /scratch/gen1/yz735/updated_coloc/lungeQTL_mucin.txt
done

hostname
date
