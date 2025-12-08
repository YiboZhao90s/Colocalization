# to automise all the analyses, need to add b37 and b38 for all QTL sets
meta_QTLs <- readRDS("/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")

# apart from decode_pQTL_mucin, all the other files need to add b37/b38 positions
#download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
#download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz","hg38ToHg19.over.chain.gz")

library(GenomicRanges)

# b37 -> b38: ukb_pQTL_mucin, eqtlgen

# snps_b37 <- data.frame("chr" = paste0("chr",meta_QTLs$ukb_pQTL_mucin$CHROM),
#                        "pos" = meta_QTLs$ukb_pQTL_mucin$GENPOS)
# gr_b37 <- GRanges(
#         seqnames = snps_b37$chr,
#         ranges = IRanges(start = snps_b37$pos, end = snps_b37$pos)
# )
# #chain_file <- "hg19ToHg38.over.chain.gz"
# #R.utils::gunzip(chain_file, overwrite = TRUE)
# chain <- import("hg19ToHg38.over.chain")
# mapped <- liftOver(gr_b37, chain)
# meta_QTLs$ukb_pQTL_mucin$GENPOSb38 <- sapply(mapped, function(x) if (length(x) > 0) start(x)[1] else NA)

snps_b37 <- data.frame("chr" = paste0("chr",meta_QTLs$eqtlgen_mucin$SNPChr),
                       "pos" = meta_QTLs$eqtlgen_mucin$SNPPos)
gr_b37 <- GRanges(
        seqnames = snps_b37$chr,
        ranges = IRanges(start = snps_b37$pos, end = snps_b37$pos)
)
#chain_file <- "hg19ToHg38.over.chain.gz"
#R.utils::gunzip(chain_file, overwrite = TRUE)
chain <- import("hg19ToHg38.over.chain")
mapped <- liftOver(gr_b37, chain)
meta_QTLs$eqtlgen_mucin$SNPPosb38 <- sapply(mapped, function(x) if (length(x) > 0) start(x)[1] else NA)

# b38 -> b37: sc-eQTL
snps_b38 <- data.frame("chr" = paste0("chr",meta_QTLs$scQTL_limix_mucin$snp_chromosome),
                       "pos" = meta_QTLs$scQTL_limix_mucin$snp_position)
gr_b38 <- GRanges(
        seqnames = snps_b38$chr,
        ranges = IRanges(start = snps_b38$pos, end = snps_b38$pos)
)
# chain_file <- "hg38ToHg19.over.chain.gz"
# R.utils::gunzip(chain_file, overwrite = TRUE)
chain <- import("hg19ToHg38.over.chain")
mapped <- liftOver(gr_b38, chain)
meta_QTLs$scQTL_limix_mucin$snp_positionb37 <- sapply(mapped, function(x) if (length(x) > 0) start(x)[1] else NA)

saveRDS(meta_QTLs, file = "/scratch/gen1/yz735/Colocalization/mucin_QTL.rds")     
