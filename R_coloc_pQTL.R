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

pQTL1 <- meta_QTLs$ukb_pQTL_mucin
pQTL2 <- meta_QTLs$decode_pQTL_mucin
pQTL2 <- separate(pQTL2, col = "source_file", into = c("header1", "header2", "PROTEIN", "PROTEIN2"),sep = "_")
library(dplyr)
fit_coloc <- function(pQTL, GWAS, build, out){
        df1 <- GWAS
        if (build == "b38"){
                s = 153763/(153763+1647022)
        } else if (build == "b37"){
                s = (5135+5414)/(25675+21471+5135+5414)
        }
        df2 <- pQTL
        pheno_list <- unique(pQTL$PROTEIN)
        n = 1
        PPH4 <- integer()
        PPH3 <- integer()
        PPH2 <- integer()
        PPH1 <- integer()
        for (protein in pheno_list){
                temp_df2 <- subset(df2, PROTEIN==protein)
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
                             MAF = ifelse(df1_trimmed$MAF==1, 0.999, df1_trimmed$MAF),
                             type = "cc",
                             s = s)
                qtl <- list(snp = overlap,
                            beta = temp_df2_trimmed$BETA,
                            varbeta = temp_df2_trimmed$VARBETA^2,
                            MAF = temp_df2_trimmed$MAF,
                            type = "quant",
                            N = mean(temp_df2_trimmed$N))
                res <- coloc.abf(gwas, qtl)
                PPH1[n] <- res$summary[3]
                PPH2[n] <- res$summary[4]
                PPH3[n] <- res$summary[5]
                PPH4[n] <- res$summary[6]
                if (res$summary[6]>=0.5){
                        output <- res$results
                        outfile <- paste(out, protein,"coloc.txt", sep = "_")
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
pres1 <- fit_coloc(pQTL1, GWAS2022, "b38", "ukb.b38")
pres2 <- fit_coloc(pQTL1, GWAS2019, "b37", "ukb.b37")
pres3 <- fit_coloc(pQTL2, GWAS2022, "b38", "decode.b38")
pres4 <- fit_coloc(pQTL2, GWAS2019, "b37", "decode.b37")
write.table(pres1, "ukb_pQTL_b38.txt")
write.table(pres2, "ukb_pQTL_b37.txt")
write.table(pres3, "decode_pQTL_b38.txt")
write.table(pres2, "decode_pQTL_b37.txt")
