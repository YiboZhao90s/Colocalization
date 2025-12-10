## Colocalization: pQTLs
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

scQTL <- meta_QTLs$scQTL_limix_mucin
scQTL <- separate(scQTL, col = "source_file", into = c("celltype"), sep = "_")
celltype <- unique(scQTL$celltype)
library(dplyr)
fit_coloc <- function(GWAS, build, cell){
        df1 <- GWAS
        if (build == "b38"){
                s = 153763/(153763+1647022)
        } else if (build == "b37"){
                s = (5135+5414)/(25675+21471+5135+5414)
        }
        df2 <- subset(scQTL, celltype==cell)
        pheno_list <- unique(df2$feature_id)
        n = 1
        PPH4 <- integer()
        PPH3 <- integer()
        PPH2 <- integer()
        PPH1 <- integer()
        for (protein in pheno_list){
                temp_df2 <- subset(df2, feature_id==protein)
                if (build == "b38"){
                        overlap <- intersect(temp_df2$SNP38, substr(df1$SNP38, 1, 12))
                        df1_trimmed <- df1[match(overlap, substr(df1$SNP38, 1, 12)),]
                        temp_df2_trimmed <- temp_df2[match(overlap, temp_df2$SNP38),]
                } else if (build == "b37"){
                        overlap <- intersect(temp_df2$SNP37, substr(df1$SNP37, 1, 12))
                        df1_trimmed <- df1[match(overlap, substr(df1$SNP37, 1, 12)),]
                        temp_df2_trimmed <- temp_df2[match(overlap, temp_df2$SNP37),]
                }
                
                gwas <- list(snp = overlap,
                             beta = df1_trimmed$BETA,
                             varbeta = df1_trimmed$VARBETA^2,
                             MAF = ifelse(df1_trimmed$MAF==1, 0.999, df1_trimmed$MAF),
                             type = "cc",
                             s = s)
                qtl <- list(snp = overlap,
                            beta = temp_df2_trimmed$BETA,
                            varbeta = temp_df2_trimmed$VARBETA^2,
                            MAF = temp_df2_trimmed$MAF,
                            type = "quant",
                            N = mean(temp_df2_trimmed$n_samples))
                res <- coloc.abf(gwas, qtl)
                PPH1[n] <- res$summary[3]
                PPH2[n] <- res$summary[4]
                PPH3[n] <- res$summary[5]
                PPH4[n] <- res$summary[6]
                if (res$summary[6]>=0.5){
                        output <- res$results
                        outfile <- paste(cell, protein,"coloc.txt", sep = "_")
                        write.table(output, paste("/scratch/gen1/yz735/Colocalization/colocres", outfile, sep = "/"), sep = "\t", quote = F, row.names = F)
                        print(paste("colocalization found for", protein, sep = " "))
                } else {
                        print("no colocalization")
                }
                n = n+1
        }
        PPtable <- cbind(PPH1,PPH2, PPH3, PPH4)
        rownames(PPtable) <- pheno_list
        return(PPtable)
}

"endothelial" "epithelial"  "immune"      "mesenchymal"
scres1 <- fit_coloc(GWAS2022, "b38", "endothelial")
scres2 <- fit_coloc(GWAS2019, "b37", "endothelial")

scres3 <- fit_coloc(GWAS2022, "b38", "epithelial")
scres4 <- fit_coloc(GWAS2019, "b37", "epithelial")

scres5 <- fit_coloc(GWAS2022, "b38", "immune")
scres6 <- fit_coloc(GWAS2019, "b37", "immune")

scres7 <- fit_coloc(GWAS2022, "b38", "mesenchymal")
scres8 <- fit_coloc(GWAS2019, "b37", "mesenchymal")

write.table(scres1, "scQTL_endothelial_b38.txt", sep = "\t", quote = F, row.names = T)
write.table(scres2, "scQTL_endothelial_b37.txt", sep = "\t", quote = F, row.names = T)
write.table(scres3, "scQTL_epithelial_b38.txt", sep = "\t", quote = F, row.names = T)
write.table(scres4, "scQTL_epithelial_b37.txt", sep = "\t", quote = F, row.names = T)
write.table(scres5, "scQTL_immune_b38.txt", sep = "\t", quote = F, row.names = T)
write.table(scres6, "scQTL_immune_b37.txt", sep = "\t", quote = F, row.names = T)
write.table(scres7, "scQTL_mesenchymal_b38.txt", sep = "\t", quote = F, row.names = T)
write.table(scres8, "scQTL_mesenchymal_b37.txt", sep = "\t", quote = F, row.names = T)

