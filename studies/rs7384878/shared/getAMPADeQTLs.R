library(EndophenotypeExplorer)
targetGene <- "PILRA"
etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=FALSE)
# 7:100334426 (GRCh38)
# 7:99932049 (GRCh37)

chrom <- "chr7"
tag.snp.hg19 <- 99932049
shoulder <- 1000000
start <- tag.snp.hg19 - shoulder
end <- tag.snp.hg19 + shoulder

tbl.eqtls <- etx$getEQTLsInRegion(chrom, start, end)
dim(tbl.eqtls)
save(tbl.eqtls, file="~/github/tms-makeAndBreak/studies/rs7384878/shared/ampad.eqtls.pilra.plusMinus.1M.RData")

