#!/bin/bash
#
# SLURM directives:
#SBATCH --job-name=extract_decode_pQTL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

cd /data/gen1/reference/decode_pQTL

> /scratch/gen1/yz735/updated_coloc/decode_pQTL_mucin.txt

for file in *b37.txt.gz; do
    echo "Processing $file..." >&2
   
for file in *b37.txt.gz; do
    echo "Processing $file..." >&2
    zcat "$file" | awk -v fname="$file" '$1 == "chr11" && $2 >= 1012240 && $2 <= 1274371 { print $0 "\t" fname }' >> /scratch/gen1/yz735/updated_coloc/decode_pQTL_mucin.txt
done

done

hostname
date
