## Colocalization with GTEx
meta_QTLs <- readRDS("/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")
GWAS <- readRDS("/scratch/gen1/yz735/Colocalization/GWAS.rds")
library(devtools)
library(data.table)
library(coloc)
library(mirrorplot)
library(dplyr)
library(tidyr)
library(arrow)
#library(qqman)
library(locuszoomr)

## Here we need to run colocalization for each gene/protein separately

## GTEx Lung 1-b38; 2-b37
pheno_list1 <- unique(meta_QTLs$GTExb38Lung_mucin$phenotype_id)
pheno_list2 <- unique(meta_QTLs$GTExb37Lung_mucin$phenotype_id)

GWAS1 <- GWAS$GWAS2022
GWAS2 <- GWAS$GWAS2019
library(dplyr)
fit_coloc <- function(tissue = "Lung", build = "b38"){
        if (build == "b38"){
                df1 <- GWAS1
                s = 153763/(153763+1647022)
        } else if (build == "b37"){
                df1 <- GWAS2
                s = (5135+5414)/(5135+5414+25675+21471)
        }
        file_qtl <- paste0("GTEx", build, tissue, "_mucin")
        df2 <- meta_QTLs[[file_qtl]]
        pheno_list <- unique(df2$phenotype_id)
        n = 1
        PPH4 <- integer()
        PPH3 <- integer()
        PPH2 <- integer()
        PPH1 <- integer()
        for (gene in pheno_list){
                temp_df2 <- subset(df2, phenotype_id==gene)
                if (build == "b38"){
                        overlap <- intersect(temp_df2$SNP38, df1$SNP38)
                        df1_trimmed <- df1[match(overlap, df1$SNP38),]
                        temp_df2_trimmed <- temp_df2[match(overlap, temp_df2$SNP38),]
                } else if (build == "b37"){
                        overlap <- intersect(temp_df2$SNP37, df1$SNP37)
                        df1_trimmed <- df1[match(overlap, df1$SNP37),]
                        temp_df2_trimmed <- temp_df2[match(overlap, temp_df2$SNP37),]
                }
                
                gwas <- list(snp = overlap,
                             beta = df1_trimmed$BETA,
                             varbeta = df1_trimmed$VARBETA^2,
                             MAF = df1_trimmed$MAF,
                             type = "cc",
                             s = s)
                qtl <- list(snp = overlap,
                            beta = temp_df2_trimmed$BETA,
                            varbeta = temp_df2_trimmed$VARBETA^2,
                            MAF = temp_df2_trimmed$MAF,
                            type = "quant",
                            N = mean(temp_df2_trimmed$ma_samples))
                res <- coloc.abf(gwas, qtl)
                PPH1[n] <- res$summary[3]
                PPH2[n] <- res$summary[4]
                PPH3[n] <- res$summary[5]
                PPH4[n] <- res$summary[6]
                if (res$summary[6]>=0.5){
                        output <- res$results
                        outfile <- paste(file_qtl, gene,"coloc.txt", sep = "_")
                        write.table(output, paste("/scratch/gen1/yz735/Colocalization/colocres", outfile, sep = "/"), sep = "\t", quote = F, row.names = F)
                        print(paste("colocalization found for", gene, sep = " "))
                } else {
                        print("no colocalization")
                }
                n = n+1
        }
        PPtable <- cbind(PPH1,PPH2, PPH3, PPH4)
        rownames(PPtable) <- pheno_list
        return(PPtable)
}

eres1 <- fit_coloc("Blood", "b38")
eres2 <- fit_coloc("Blood", "b37")
eres3 <- fit_coloc("Lung", "b38")
eres4 <- fit_coloc("Lung", "b37")

write.table(eres4, "GTEx_Lung_b37.txt", sep = "\t", quote = F, row.names = T)
