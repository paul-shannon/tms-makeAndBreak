library(ebi.eqtls)
ee <- ebi.eqtls$new()
rsid.loc <- "chr7:100,334,425"
tbl.cat <- ee$getCatalog()
gtex.brain.tissues <- unique(grep("GTEx_V8.Brain", tbl.cat$unique_id, v=TRUE))
toi <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
         "GTEx_V8.Brain_Cerebellum",
         "GTEx_V8.Brain_Cortex",
         "GTEx_V8.Brain_Frontal_Cortex_BA9",
         "GTEx_V8.Brain_Hippocampus",
         "GTEx_V8.Brain_Hypothalamus")

fetch <- function(tissue){
    ee$fetch.eqtls.in.chunks(chrom="chr7",
                             start=99334425,
                             end=101334426,
                             study=tissue,
                             simplify=TRUE,
                             chunk.size=5000)
    }
tbls <- lapply(toi, fetch)

tbl.eqtl.gtex <- do.call(rbind, tbls)
rownames(tbl.eqtl.gtex) <- NULL
dim(tbl.eqtl.gtex)
head(tbl.eqtl.gtex)
as.data.frame(sort(table(subset(tbl.eqtl.gtex, pvalue < 1e-20)$gene)))
save(tbl.eqtl.gtex, file="gtex-eqtls-tbl.RData")

