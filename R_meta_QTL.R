# Create meta-QTL file
# It includes: summary statistics, SNP position, sample size, related gene/protein, QTL type, need position coloumn with both b38 and b37 to overlap with GWAS2019 and GWAS2020
library(data.table)
library(tidyr)

## GTEx eQTLs: Lung and Whole-blood
path_GTEx = "/data/gen1/Kayesha/resources/eQTL/GTeX"

file_GTExb38Lung = "Lung.v8.EUR.allpairs.chr11.hg38.txt.gz"
file_GTExb38Blood = "Whole_Blood.v8.EUR.allpairs.chr11.hg38.txt.gz"
file_GTExb37Lung = "Lung.v8.EUR.allpairs.chr11.hg19.txt.gz"
file_GTExb37Blood = "Whole_Blood.v8.EUR.allpairs.chr11.hg19.txt.gz"

GTExb38Lung <- fread(paste(path_GTEx, file_GTExb38Lung, sep = "/"))
GTExb38Blood <- fread(paste(path_GTEx, file_GTExb38Blood, sep = "/"))
GTExb38Lung_mucin <- subset(GTExb38Lung, POS>1012240& POS<1295601)
GTExb38Blood_mucin <- subset(GTExb38Blood, POS>1012240&POS<1295601)
rm(GTExb38Blood)
rm(GTExb38Lung)

GTExb37Lung <- fread(paste(path_GTEx, file_GTExb37Lung, sep = "/"))
GTExb37Blood <- fread(paste(path_GTEx, file_GTExb37Blood, sep = "/"))
GTExb37Lung_mucin <- subset(GTExb37Lung, POS>1012240& POS<1274370)
GTExb37Blood_mucin <- subset(GTExb37Blood, POS>1012240&POS<1274370)
rm(GTExb37Blood)
rm(GTExb37Lung)

#write.table(GTExb37Blood_mucin, "/scratch/gen1/yz735/Colocalization/GTExb37Blood_mucin.txt", sep = "\t", quote = F, row.names = F)
#write.table(GTExb37Lung_mucin, "/scratch/gen1/yz735/Colocalization/GTExb37Lung_mucin.txt", sep = "\t", quote = F, row.names = F)

ukb_pQTL_mucin <- read.delim("/scratch/gen1/yz735/Colocalization/ukb_pQTL_mucin.txt")

# lungeQTL_mucin <- read.delim("/scratch/gen1/yz735/Colocalization/lungeQTL_mucin.txt", header = F)
# colnames(lungeQTL_mucin) <- c("#Probe",	"MarkerName","CHR","BP","Allele1","Allele2","Freq1","Effect","StdErr","P","Direction") # position doesn't match with b38 or b37.

#rm(lungeQTL_mucin)
decode_pQTL_mucin <- read.delim("/scratch/gen1/yz735/Colocalization/decode_pQTL_mucin.txt", header = F)
colnames(decode_pQTL_mucin) <-  c("CHR", "Pos_b38", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "minus_long10_pval", "SE", "N", "ImpMAF", "Pos_b37", "source_file") 

##sc QTL (limix)
scQTL_limix_mucin <- read.delim("/scratch/gen1/yz735/Colocalization/scQTL_limix_mucin.txt", header = F)
colnames(scQTL_limix_mucin) <- c("feature_id", "snp_id", "p_value", "beta", "beta_se", "empirical_feature_p_value", "feature_chromosome", "feature_start", "feature_end", "n_samples", "n_e_samples", "alpha_param", "beta_param", "rho", "snp_chromosome", "snp_position", "assessed_allele", "call_rate", "maf", "hwe_p", "source_file")

## eqtlgen
eqtlgen_mucin <- read.delim("/scratch/gen1/yz735/Colocalization/eqtlgen_mucin.txt", header = F)
colnames(eqtlgen_mucin) <- c("Pvalue", "SNP" , "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene", "GeneSymbol", "GeneChr", "GenePos", "NrCohorts", "NrSamples", "FDR", "BonferroniP")

mucin_QTL <- list("GTExb38Lung_mucin" = GTExb38Lung_mucin, 
                  "GTExb38Blood_mucin" = GTExb38Blood_mucin, 
                  "GTExb37Lung_mucin" = GTExb37Lung_mucin, 
                  "GTExb37Blood_mucin" = GTExb37Blood_mucin, 
                  "ukb_pQTL_mucin" = ukb_pQTL_mucin, 
                  "decode_pQTL_mucin" = decode_pQTL_mucin, 
                  "scQTL_limix_mucin" = scQTL_limix_mucin, 
                  "eqtlgen_mucin" = eqtlgen_mucin)
saveRDS(mucin_QTL, file = "/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")                          
