library(ebi.eqtls)
ee <- ebi.eqtls$new()
tag.snp <- "rs143080277"
tag.snp.chrom <- "chr2"
tag.snp.hg38 <- 105749599
shoulder <- 500000
toi <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
         "GTEx_V8.Brain_Cerebellum",
         "GTEx_V8.Brain_Cortex",
         "GTEx_V8.Brain_Frontal_Cortex_BA9",
         "GTEx_V8.Brain_Hippocampus",
         "GTEx_V8.Brain_Hypothalamus")

fetch <- function(tissue){
    ee$fetch.eqtls.in.chunks(chrom=tag.snp.chrom,
                             start=tag.snp.hg38 - shoulder,
                             end=tag.snp.hg38 + shoulder,
                             study=tissue,
                             simplify=TRUE,
                             chunk.size=5000)
    }
tbls <- lapply(toi, fetch)

tbl.eqtl.gtex <- do.call(rbind, tbls)
rownames(tbl.eqtl.gtex) <- NULL
dim(tbl.eqtl.gtex)
head(tbl.eqtl.gtex)
tbl.eqtl.gtex$score <-with(tbl.eqtl.gtex, -log10(pvalue) * beta)

as.data.frame(sort(table(subset(tbl.eqtl.gtex, pvalue < 1e-5)$gene)))
save(tbl.eqtl.gtex, file="gtex-eqtls-tbl.RData")

