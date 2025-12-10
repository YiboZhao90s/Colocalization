# reformat for colocalization

# Load GWAS
GWAS2019 <- read.delim("/scratch/gen1/yz735/Colocalization/GWAS2019_mucin.txt", header = F)
colnames(GWAS2019) <- c("SNP", "CHR", "POS_b37", "Coded", "nonCoded", "Coded_freq", "INFO", "beta", "SE_GC","OR" , "OR_L95", "OR_U95", "P_GC") 									
GWAS2019$SNP37 <- paste(GWAS2019$CHR, GWAS2019$POS_b37, GWAS2019$nonCoded, GWAS2019$Coded, sep = "_")
colnames(GWAS2019)[6] <- "MAF"
colnames(GWAS2019)[9] <- "VARBETA"
colnames(GWAS2019)[13] <- "PVAL"
colnames(GWAS2019)[8] <- "BETA"


GWAS2022 <- read.delim("/scratch/gen1/yz735/Colocalization/GWAS2022_mucin.txt", header = F)
colnames(GWAS2022) <- c("CHR","POS","REF","ALT","rsid","all_meta_AF","inv_var_meta_beta",	"inv_var_meta_sebeta",	"inv_var_meta_p","inv_var_het_p","direction","N_caseN_ctrl","n_dataset","n_bbk", "n_bbk1","is_strand_flip","is_diff_AF_gnomAD")
GWAS2022$SNP38 <- paste(GWAS2022$CHR, GWAS2022$POS, GWAS2022$REF, GWAS2022$ALT, sep = "_")
colnames(GWAS2022)[6] <- "MAF"
colnames(GWAS2022)[7] <- "BETA"
colnames(GWAS2022)[8] <- "VARBETA"

# Load QTL
meta_QTLs <- readRDS("/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")

# prepare matrix for colocalization
summary(meta_QTLs)

## GTEx
head(meta_QTLs$GTExb38Lung_mucin)
colnames(meta_QTLs$GTExb38Lung_mucin)[7] <-"MAF"
colnames(meta_QTLs$GTExb38Lung_mucin)[10] <-"PVAL"
colnames(meta_QTLs$GTExb38Lung_mucin)[11] <-"BETA"
colnames(meta_QTLs$GTExb38Lung_mucin)[12] <-"VARBETA"
colnames(meta_QTLs$GTExb38Lung_mucin)[14] <-"SNP38"

head(meta_QTLs$GTExb38Blood_mucin)
colnames(meta_QTLs$GTExb38Blood_mucin)[7] <-"MAF"
colnames(meta_QTLs$GTExb38Blood_mucin)[10] <-"PVAL"
colnames(meta_QTLs$GTExb38Blood_mucin)[11] <-"BETA"
colnames(meta_QTLs$GTExb38Blood_mucin)[12] <-"VARBETA"
colnames(meta_QTLs$GTExb38Blood_mucin)[14] <-"SNP38"

head(meta_QTLs$GTExb37Lung_mucin)
colnames(meta_QTLs$GTExb37Lung_mucin)[7] <-"MAF"
colnames(meta_QTLs$GTExb37Lung_mucin)[10] <-"PVAL"
colnames(meta_QTLs$GTExb37Lung_mucin)[11] <-"BETA"
colnames(meta_QTLs$GTExb37Lung_mucin)[12] <-"VARBETA"
colnames(meta_QTLs$GTExb37Lung_mucin)[4] <-"SNP37"

head(meta_QTLs$GTExb37Blood_mucin)
colnames(meta_QTLs$GTExb37Blood_mucin)[7] <-"MAF"
colnames(meta_QTLs$GTExb37Blood_mucin)[10] <-"PVAL"
colnames(meta_QTLs$GTExb37Blood_mucin)[11] <-"BETA"
colnames(meta_QTLs$GTExb37Blood_mucin)[12] <-"VARBETA"
colnames(meta_QTLs$GTExb37Blood_mucin)[4] <-"SNP37"

# other QTLs
head(meta_QTLs$ukb_pQTL_mucin)
meta_QTLs$ukb_pQTL_mucin$SNP37 <- paste(meta_QTLs$ukb_pQTL_mucin$CHROM, 
                                        meta_QTLs$ukb_pQTL_mucin$GENPOS, 
                                        meta_QTLs$ukb_pQTL_mucin$ALLELE0,
                                        meta_QTLs$ukb_pQTL_mucin$ALLELE1,
                                        sep = "_")
meta_QTLs$ukb_pQTL_mucin$SNP38 <- paste(meta_QTLs$ukb_pQTL_mucin$CHROM, 
                                        meta_QTLs$ukb_pQTL_mucin$GENPOSb38, 
                                        meta_QTLs$ukb_pQTL_mucin$ALLELE0,
                                        meta_QTLs$ukb_pQTL_mucin$ALLELE1,
                                        sep = "_")
colnames(meta_QTLs$ukb_pQTL_mucin)[12] <- "VARBETA"
colnames(meta_QTLs$ukb_pQTL_mucin)[8] <- "MAF"
meta_QTLs$ukb_pQTL_mucin$PVAL <- -10^meta_QTLs$ukb_pQTL_mucin$LOG10P

head(meta_QTLs$decode_pQTL_mucin)
meta_QTLs$decode_pQTL_mucin$SNP37 <- paste(substr(meta_QTLs$decode_pQTL_mucin$CHR, 4,5), 
                                        meta_QTLs$decode_pQTL_mucin$Pos_b37, 
                                        meta_QTLs$decode_pQTL_mucin$otherAllele,
                                        meta_QTLs$decode_pQTL_mucin$effectAllele,
                                        sep = "_")
meta_QTLs$decode_pQTL_mucin$SNP38 <- paste(substr(meta_QTLs$decode_pQTL_mucin$CHR, 4,5), 
                                           meta_QTLs$decode_pQTL_mucin$Pos_b38, 
                                           meta_QTLs$decode_pQTL_mucin$otherAllele,
                                           meta_QTLs$decode_pQTL_mucin$effectAllele,
                                           sep = "_")
colnames(meta_QTLs$decode_pQTL_mucin)[7] <- "BETA"
colnames(meta_QTLs$decode_pQTL_mucin)[10] <- "VARBETA"
colnames(meta_QTLs$decode_pQTL_mucin)[8] <- "PVAL"
colnames(meta_QTLs$decode_pQTL_mucin)[12] <- "MAF"

head(meta_QTLs$scQTL_limix_mucin)
## in scQTL, no alt allele was provided. SNP is will be CHR_POS_REF only
meta_QTLs$scQTL_limix_mucin$SNP37 <- paste(meta_QTLs$scQTL_limix_mucin$snp_chromosome, 
                                           meta_QTLs$scQTL_limix_mucin$snp_positionb37, 
                                           meta_QTLs$scQTL_limix_mucin$assessed_allele,
                                           sep = "_")
meta_QTLs$scQTL_limix_mucin$SNP38 <- paste(meta_QTLs$scQTL_limix_mucin$snp_chromosome, 
                                           meta_QTLs$scQTL_limix_mucin$snp_position, 
                                           meta_QTLs$scQTL_limix_mucin$assessed_allele,
                                           sep = "_")
colnames(meta_QTLs$scQTL_limix_mucin)[4] <- "BETA"
colnames(meta_QTLs$scQTL_limix_mucin)[5] <- "VARBETA"
colnames(meta_QTLs$scQTL_limix_mucin)[6] <- "PVAL"
colnames(meta_QTLs$scQTL_limix_mucin)[19] <- "MAF"

head(meta_QTLs$eqtlgen_mucin)

# eqtlgen_mucin doesn't have beta, se and maf wait and see whether it has overlapped SNPs with GWAS then extract
meta_QTLs$eqtlgen_mucin$SNP37 <- paste(meta_QTLs$eqtlgen_mucin$SNPChr, 
                                           meta_QTLs$eqtlgen_mucin$SNPPos, 
                                           meta_QTLs$eqtlgen_mucin$OtherAllele,
                                           meta_QTLs$eqtlgen_mucin$AssessedAllele,
                                           sep = "_")
meta_QTLs$eqtlgen_mucin$SNP38 <- paste(meta_QTLs$eqtlgen_mucin$SNPChr, 
                                       meta_QTLs$eqtlgen_mucin$SNPPosb38, 
                                       meta_QTLs$eqtlgen_mucin$OtherAllele,
                                       meta_QTLs$eqtlgen_mucin$AssessedAllele,
                                       sep = "_")
overlap1 <- intersect(meta_QTLs$eqtlgen_mucin$SNP38, GWAS2022$SNP38) 
overlap2 <- intersect(meta_QTLs$eqtlgen_mucin$SNP37, GWAS2019$SNP37)
library(biomaRt)
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
rsid <- meta_QTLs$eqtlgen_mucin$SNP
snp_info <- getBM(
                attributes = c("refsnp_id", "allele", "minor_allele_freq", "minor_allele"),
                filters = "snp_filter",
                values = rsid,
                mart = snp_mart)
MAF <- integer()
for (i in 1:nrow(meta_QTLs$eqtlgen_mucin)){
       posi <- match(meta_QTLs$eqtlgen_mucin$SNP[i], snp_info$refsnp_id)
       MAF[i] <- snp_info$minor_allele_freq[posi]
       
}
meta_QTLs$eqtlgen_mucin$MAF <- MAF
meta_QTLs$eqtlgen_mucin$VarG <- 2 * meta_QTLs$eqtlgen_mucin$MAF * (1 - meta_QTLs$eqtlgen_mucin$MAF)
meta_QTLs$eqtlgen_mucin$SE <- 1 / sqrt(meta_QTLs$eqtlgen_mucin$NrSamples * meta_QTLs$eqtlgen_mucin$VarG)
meta_QTLs$eqtlgen_mucin$BETA <- meta_QTLs$eqtlgen_mucin$Zscore * meta_QTLs$eqtlgen_mucin$SE
colnames(meta_QTLs$eqtlgen_mucin)[21] <- "VARBETA"

saveRDS(meta_QTLs, file = "/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")   
GWAS <- list("GWAS2019" = GWAS2019, "GWAS2022" = GWAS2022)
saveRDS(GWAS, file = "/scratch/gen1/yz735/Colocalization/GWAS.rds")   
