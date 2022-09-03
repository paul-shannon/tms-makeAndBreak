library(trena)
library(ghdb)
library(RUnit)

targetGene <- "NDUFS2"

    #-------------------------------
    # initialize with fimo
    #-------------------------------
ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg38")
tbl.fimo <- get(load("../shared/tbl.fimo.NDUFS2.RData"))
p.value.col <- grep("^p\\.value$", colnames(tbl.fimo))
if(length(p.value.col) == 1)
    colnames(tbl.fimo)[p.value.col] <- "fimo.pval"

ft$setFundamentalRegions(tbl.fimo)

    #-----------------------------------------------------
    #  add rosmap eqtls
    #  TODO: these eqtls do not have betas.  we need them
    #-----------------------------------------------------

library(EndophenotypeExplorer)
etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
tbl.eqtl <- etx$get.ampad.EQTLsForGene()
dim(tbl.eqtl)
tbl.eqtl.sub <- subset(tbl.eqtl, pvalue < 0.1 & study=="ampad-rosmap")
dim(tbl.eqtl.sub)

feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue")
checkTrue(all(as.character(feature.guide) %in% colnames(tbl.eqtl.sub)))

default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1)
tbl.eqtl.sub$start <- tbl.eqtl.sub$hg38 - 1
tbl.eqtl.sub$end <- tbl.eqtl.sub$hg38

ft$addRegionFeature(tbl.eqtl.sub, feature.genome="hg38", feature.guide, default.values)

    #----------------------------------
    # obtain GTEx brain tissue eQTLs
    #----------------------------------

library(ebi.eqtls)
ee <- ebi.eqtls$new()
tbl.cat <- ee$getCatalog()
gtex.brain.tissues <- grep("GTEx_V8.Brain_", tbl.cat$unique_id, v=TRUE)

    # choose a few
toi <- c("GTEx_V8.Brain_Amygdala",              "GTEx_V8.Brain_Anterior_cingulate_cortex_BA24",
         "GTEx_V8.Brain_Cerebellar_Hemisphere", "GTEx_V8.Brain_Cerebellum",
         "GTEx_V8.Brain_Cortex",                "GTEx_V8.Brain_Frontal_Cortex_BA9",
         "GTEx_V8.Brain_Hippocampus",           "GTEx_V8.Brain_Hypothalamus")

chrom <- tbl.fimo$chrom[1]
start <- min(tbl.fimo$start)
end   <- max(tbl.fimo$end)

chrom
start
end


tbls <- list()
for(tissue in toi){  # took 2-3 hours
   tbl.eqtls <- ee$fetch.eqtls.in.chunks(chrom=chrom,
                                         start=start,
                                         end=end,
                                         study=tissue,
                                         simplify=TRUE,
                                         chunk.size=10000)
   tbls[[tissue]] <- tbl.eqtls
   }

tbl.eqtls.all <- do.call(rbind, tbls)
rownames(tbl.eqtls.all) <- NULL
dim(tbl.eqtls.all)   # 1690734       8
if(!grepl("chr", tbl.eqtls.all$chrom[1]))
   tbl.eqtls.all$chrom <- paste0("chr", tbl.eqtls.all$chrom)
head(tbl.eqtls.all)

filename <- sprintf("gtex.brain.eqtls.%s-%d-%d.RData", chrom, start, end)
full.path <- file.path("~/github/tms-makeAndBreak/studies/rs4575098-ndufs2-region/shared", filename)
save(tbl.eqtls.all, file=full.path)

    #---------------------------------------------
    # load gtex brain eqtls to the feature table
    # these come without genomic coordinaes, so
    # load in rsid hg38 hg19 locs from ../shared
    #---------------------------------------------

tbl.snp.locs <- get(load("../shared/tbl.rsid.160684126-161688825.hg38.hg19.RData"))

for(brain.tissue in unique(tbl.eqtls.all$id)){
    short.name <- sub("V8.Brain_", "", brain.tissue)
    tbl.sub <- subset(tbl.eqtls.all, id==brain.tissue & gene == targetGene & pvalue < 0.1)
    printf("%s: %d", short.name, nrow(tbl.sub))
    tbl.load <- tbl.sub[, c("chrom", "hg38", "hg38", "rsid", "gene", "pvalue", "beta", "id")]
    colnames(tbl.load) <- c("chrom", "start", "end", "rsid", "gene", "pvalue", "beta", "tissue")
    tbl.load$start <- tbl.load$start - 1
    feature.guide <- list("rsid", "pvalue", "beta")
    feature.names <- sprintf("%s.%s", short.name, feature.guide)
    names(feature.guide) <- feature.names
    default.values <- list("", 1, 0)
    names(default.values) <- feature.names
    ft$addRegionFeature(tbl.load, feature.genome="hg38", feature.guide, default.values)
    }

tbl.ft <- ft$getTable()

threshold <- 0.025
dim(subset(tbl.ft, rosmap.eqtl.pvalue < threshold &
                   GTEx_Amygdala.pvalue < threshold &
                   GTEx_Anterior_cingulate_cortex_BA24.pvalue < threshold &
                   GTEx_Cerebellar_Hemisphere.pvalue < threshold &
                   GTEx_Cerebellum.pvalue < threshold &
                   GTEx_Cortex.pvalue < threshold &
                   GTEx_Frontal_Cortex_BA9.pvalue < threshold &
                   GTEx_Hippocampus.pvalue < threshold &
                   GTEx_Hypothalamus.pvalue < threshold))

checkTrue(all(as.character(feature.guide) %in% colnames(tbl.eqtl.sub)))

default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1)
tbl.eqtl.sub$start <- tbl.eqtl.sub$hg38 - 1
tbl.eqtl.sub$end <- tbl.eqtl.sub$hg38

ft$addRegionFeature(tbl.eqtl.sub, feature.genome="hg38", feature.guide, default.values)


    #-------------------------------
    #
    #-------------------------------


    #-------------------------------
    #
    #-------------------------------


    #-------------------------------
    #
    #-------------------------------


