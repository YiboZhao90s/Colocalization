# Colocalization
Colocalization Analysis on Moderate-Severe Asthma GWAS signals and xQTLs. This is part of WP1 of MRC grant (the MuCoSA project) at University of Leicester

# Data Description (extract_XXX.sh)
## QTL
- GTExb38Lung_mucin N = 99268 (b38)
- GTExb38Blood_mucin N = 109537 (b38)
- GTExb37Lung_mucin N = 83598 (b37)
- GTExb37Blood_mucin N = 76182 (b37)
- ukb_pQTL_mucin N = 1436909 (b37)
- decode_pQTL_mucin N = 4752089 (b38 & b37)
- sc-eQTL for 38 cell types (b38)
- lungeQTL_mucin N = 165826 (unknown source and reference genome, removed)
- eqtlgen N = 87377 (b37)

## GWAS 
- GWAS2019 N = 4237 (b37)
- GWAS2022 N = 3223 (b38)

# Main Workflow
1. Generate a "meta-QTL" including all QTL files (R_meta_QTL.R)
2. Add b38 or b37-based position into QTL files and reformat for automatic pipeline (R_lifeover_QTL.R, R_reformat.R)
3. Colocalization between QTLs and 2 GWAS separately (R_coloc_GTEx.R, R_coloc_pQTL.R, R_coloc_scQTL.R)
