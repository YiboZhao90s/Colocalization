## extract significant SNPs within the mucin region 
## GWAS-GBMI (b38)
zcat /data/gen1/yz735/updated_coloc/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | awk '$1==11 && $2>1012241 && $2<1274371' > /scratch/gen1/yz735/Colocalization/GWAS2022_mucin.txt

## GWAS-2019 (b37)
zcat /data/gen1/yz735/updated_coloc/Shrine_30552067_moderate-severe_asthma.txt.gz | awk '$2==11 && $3>1012241 && $3<1295601' > /scratch/gen1/yz735/Colocalization/GWAS2019_mucin.txt
### 51 SNPs

## liftover b37 list to b38 and merge them (see R_liftover_b37tob38.R and R_merge_sig_GWAS_b38.R)
