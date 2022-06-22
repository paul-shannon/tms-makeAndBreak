library(TrenaProjectAD)
library(EndophenotypeExplorer)
library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
source("~/github/tms-makeAndBreak/R/tmsMB-class.R")
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

tag.snp <- "rs7384878"
tag.snp.chrom <- "chr7"
tag.snp.loc <- 100334426
rsid.oi <- tag.snp

gtex.brain.tissues <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
                        "GTEx_V8.Brain_Cerebellum",
                        "GTEx_V8.Brain_Cortex",
                        "GTEx_V8.Brain_Frontal_Cortex_BA9",
                        "GTEx_V8.Brain_Hippocampus",
                        "GTEx_V8.Brain_Hypothalamus")


if(!exists("igv") & FALSE){
    igv <- start.igv("all")
    tbl.track <- data.frame(chrom=tag.snp.chrom,
                            start=tag.snp.loc-1,
                            end=tag.snp.loc,
                            name=tag.snp,
                            stringsAsFactors=FALSE)
    track <- DataFrameAnnotationTrack(tag.snp, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    }

if(!exists("tbl.haploreg")){
    data.dir <- "../shared"
    tbl.haploreg <-
       read.table(file.path(data.dir, "haploreg.tsv"), sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
    tbl.haploreg$score <- 0.9999999999 - tbl.haploreg$Rsquared
    if(FALSE){
       track <- GWASTrack("LD", tbl.haploreg, chrom.col=1, pos.col=2, pval.col=8, trackHeight=100)
       displayTrack(igv, track)
       }
    }

if(!exists("tbl.bellenquez")){
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
    file.exists(data.dir)
    filename <- "bellenguez-2022-83variants-hg38.RData"
    full.path <- file.path(data.dir, filename)
    tbl.bellenquez <- get(load(full.path))
    subset(tbl.bellenquez, rsid==rsid.oi)
    filename <- "posthuma-2019-with-hg38-locs.RData"
    full.path <- file.path(data.dir, filename)
    tbl.posthuma <- get(load(full.path))
    subset(tbl.posthuma, rsid==rsid.oi)
    }


if(!exists("tbl.gwas3")){ # this 235 row, 215 rsid table, combines posthuma, schartzentruber, bellenquez
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
    file.exists(data.dir)
    filename <- "tbl.gwas.3studies.235x5.RData"
    full.path <- file.path(data.dir, filename)
    tbl.gwas3 <- get(load(full.path))
    stopifnot(colnames(tbl.gwas3)[c(1,2,5)] == c("chrom", "hg38", "pvalue"))
    head(tbl.gwas3)
    }

if(!exists("tbl.gwascat.ad")){
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
   file.exists(data.dir)
   file <- "alzheimerSubsetOfGWASCatatalog-12jun2022.RData"
   full.path <- file.path(data.dir, file)
   file.exists(full.path)
   tbl.gwascat.ad <- get(load(full.path))
   stopifnot(colnames(tbl.gwascat.ad)[c(1,2,33)] == c("seqnames", "start", "P.VALUE"))
   if(FALSE){
       track <- GWASTrack("gwascat.ad", tbl.gwascat.ad, chrom.col=1, pos.col=2, pval.col=33,
                          trackHeight=100)
       displayTrack(igv, track)
       }
   }

# if(!exists("tbl.rsids")){  # with mayo & rosmap scores for association of each variant to extreme braak score
#    data.dir <- "~/github/gwasExplorer/explore/newScore/doesLinkagePredictBraakscSeparatedEnrichment"
#    filename <- "tbl.rsids.summary.Tue.May.31.2022-13:26:55-hetAndHomSeparate.RData"
#    full.path <- file.path(data.dir, filename)
#    stopifnot(file.exists(full.path))
#    tbl.rsids <- get(load(full.path))
#    }

if(!exists("tbl.rosmap.nearby")){
   data.dir <- "../shared"
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.rosmap.nearby.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   dim(tbl.rosmap.nearby.raw)
   tbl.rosmap.nearby <- subset(tbl.rosmap.nearby.raw, pvalue < 0.001)
   dim(tbl.rosmap.nearby)   # 6262 12
   if(FALSE){
      track <- GWASTrack("rosmap eQTLs", tbl.rosmap.nearby, trackHeight=100,
                          chrom.col=1, pos.col=3, pval.col=5)
      displayTrack(igv, track)
      }
   } # tbl.rosmap.nearby

if(!exists("tbl.gtex.eqtls")){
   data.dir <- "../shared"
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   if(!grepl("chr", tbl.gtex.eqtls.raw$chrom[1]))
      tbl.gtex.eqtls.raw$chrom <- paste0("chr", tbl.gtex.eqtls.raw$chrom)
   dim(tbl.gtex.eqtls.raw)
   targetGene <- "ZSCAN21"
   tbl.track <- subset(tbl.gtex.eqtls.raw, pvalue < 0.001 & gene==targetGene)
   dim(tbl.track)
   if(FALSE){
      track <- GWASTrack("GTEx eQTLs", tbl.track, trackHeight=100,
                          chrom.col=7, pos.col=8, pval.col=2)
      displayTrack(igv, track)
      }
   } # tbl.rosmap.nearby


if(!exists("tbl.gwascat.ad")){
   data.dir <- system.file(package="igvR", "extdata", "gwas")
   file.exists(data.dir)
   file <- "alzheimerSubsetOfGWASCatatalog-29apr2022.RData"
   full.path <- file.path(data.dir, file)
   file.exists(full.path)
   tbl.gwascat.ad <- get(load(full.path))
   stopifnot(colnames(tbl.gwascat.ad)[c(1,2,33)] == c("seqnames", "start", "P.VALUE"))
   dim(tbl.gwascat.ad) # [1] 2141 43
   }

if(!exists("tbl.fimo")){
    tbl.fimo <- get(load(file.path("../shared", "tbl.fimo.PILRA.RData")))
    gr.fimo <- GRanges(tbl.fimo)
    }
#snpLocs.initialized <- FALSE
#----------------------------------------------------------------------------------------------------
initialize.snpLocs <- function()   # load into memory just once, saving time in motifbreakR prep
{
    t0 <- system.time(x0 <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, "rs7796006"))
    printf("--- snpLocs initialized in %d seconds", as.integer(t0[["elapsed"]]))
    return(TRUE)
}
#----------------------------------------------------------------------------------------------------
showTracks <- function()
{
    showGenomicRegion(igv, "all")
    tbl.track <- GWASTrack("AD gwascat", tbl.gwascat.ad, chrom.col=1, pos.col=2, pval.col=33)
    displayTrack(igv, tbl.track)

    rosmap.outliers <- rownames(subset(tbl.rsids, rosmap.f > 2))
    tbl.track <- subset(tbl.gwas3, rsid %in% rosmap.outliers)
    track <- GWASTrack("rosmap outliers", tbl.track,  chrom.col=1, pos.col=2, pval.col=5)
    displayTrack(igv, track)

    mayo.outliers <- rownames(subset(tbl.rsids, mayo.f > 2))
    tbl.track <- subset(tbl.gwas3, rsid %in% mayo.outliers)
    track <- GWASTrack("mayo outliers", tbl.track,  chrom.col=1, pos.col=2, pval.col=5)
    displayTrack(igv, track)

    track <- GWASTrack("rosmap eQTLs", tbl.rosmap.nearby, trackHeight=100,
                       chrom.col=1, pos.col=3, pval.col=5)
    displayTrack(igv, track)

} # showTracks
#----------------------------------------------------------------------------------------------------
identify.affected.genes <- function()
{
   roi <- getGenomicRegion(igv)  # 45k between 100,312,189  and  100,356,932

   fivenum(tbl.rosmap.nearby$pvalue)
   tbl.sub <-
     subset(tbl.rosmap.nearby,
            hg38 >= roi$start & hg38 <= roi$end &
            pvalue < 1e-5)
   goi <- unique(tbl.sub$genesymbol)
   length(goi)  # 13

   tbl.freq <- as.data.frame(sort(table(tbl.sub$genesymbol), decreasing=TRUE))
   colnames(tbl.freq) <- c("gene", "eqtl.count")

   scoreGene <- function(gene){
     #tbl.gene <- subset(tbl.rosmap.nearby, genesymbol==gene)
     tbl.gene <- subset(tbl.sub, genesymbol==gene)
     round(sum(abs(-log10(tbl.gene$pvalue) * tbl.gene$beta)), digits=2)
     }

   scores <- lapply(tbl.freq$gene, scoreGene)
   names(scores) <- tbl.freq$gene
   tbl.scores <- data.frame(gene=tbl.freq$gene, sum.sig.x.beta=unlist(scores), stringsAsFactors=FALSE)

   tbl <- merge(tbl.scores, tbl.freq, by="gene")
   tbl$mean <- round(tbl$sum.sig.x.beta/tbl$eqtl.count, digits=2)
   new.order <- order(tbl$sum.sig.x.beta, decreasing=TRUE)
   tbl <- tbl[new.order,]
   tbl$gene <- as.character(tbl$gene)
   rownames(tbl) <- NULL

   tbl

} # identify.affected.genes
#----------------------------------------------------------------------------------------------------
pilra.neighborhood <- function()
{
    showGenomicRegion(igv, "PILRA")
    for(i in 1:4) zoomOut(igv)

} # pilra.neighborhood
#----------------------------------------------------------------------------------------------------
build.model <- function(targetGene, gtex.tissue, roi,
                        eqtl.pval.threshold, rna.correlation.threshold,
                        tbl.oc, known.snps)

{
   printf("--- building model of %s in %s over %dk", targetGene, gtex.tissue, as.integer(roi$width/1000))
   trenaProject <- TrenaProjectAD()
   data.dir <- "~/github/tms-makeAndBreak/studies/rs7384878/shared"
   tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.PILRA.RData")))
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   tbl.gtex.eqtls.raw$chrom <- paste0("chr", tbl.gtex.eqtls.raw$chrom)
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.001)
   dim(tbl.gtex.eqtls)
       # is the requested tissue present in the filtered tbl.gtex.eqtls?
   if(!gtex.tissue %in% tbl.gtex.eqtls$id){
       printf("--- no eqtls for %s in tissue %s", targetGene, gtex.tissue)
       return(list())
       }
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.rosmap.eqtls.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   tbl.rosmap.eqtls <- subset(tbl.rosmap.eqtls.raw, pvalue < 0.001)
   dim(tbl.rosmap.eqtls)   # 6262 12

   tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls,
                    tbl.rosmap.eqtls, tbl.oc, known.snps)

   trenaProject <- TrenaProjectAD()
   tms$setStudyRegion(chrom=roi$chrom, start=roi$start, end=roi$end)

   study.region <- tms$getStudyRegion()
   checkTrue(gtex.tissue %in% gtex.brain.tissues)
   tms$set.current.GTEx.eqtl.tissue(gtex.tissue)

   tms$run.tms()
   tms$add.eqtls.toTmsTable()

   tbl.tms <- tms$get.tmsTable()
   tbl.tms.filtered <-
      subset(tbl.tms, ampad.eqtl.pval < 1 & gtex.eqtl.pval < eqtl.pval.threshold & abs(cor.all) >= rna.correlation.threshold )
   dim(tbl.tms.filtered)
   length(unique(tbl.tms.filtered$tf))
   tms$set.tmsFilteredTable(tbl.tms.filtered)
   # tms$get.tmsFilteredTable()
   tf.candidates <- unique(tbl.tms.filtered$tf)
   printf("tf.candidates: %d", length(tf.candidates))

   tbl.trena <- data.frame()
   tbl.breaks <- data.frame()

   if(length(tf.candidates) > 0){
      tms$run.trena(tf.candidates)
      tbl.trena <- tms$get.trenaTable()
      head(tbl.trena, n=10)
      tms$breakMotifs(head(tbl.trena, n=12), tbl.tms.filtered, tbl.rosmap.eqtls)
      tbl.breaks <- tms$get.breaksTable()
      if(nrow(tbl.breaks) > 0){
           # we are (mostly) interested in the strongest breaks, high negative pctDelta
         new.order <- order(tbl.breaks$pctDelta, decreasing=FALSE)
         tbl.breaks <- tbl.breaks[new.order,]
         } # tbl.breaks
      }

   new.known.snps <- tms$get.knownSnps()
   #browser()
   result <- list(#tms=tbl.tms,
                  tms.filtered=tbl.tms.filtered,
                  trena=tbl.trena,
                  tf.candidates=tf.candidates,
                  known.snps=new.known.snps,
                  tbl.breaks=tbl.breaks,
                  breaks=tms$get.motifBreaks(),
                  roi=roi,
                  eqtl.threshold=eqtl.pval.threshold,
                  rna.correlation.threshold=rna.correlation.threshold
                  #tbl.gtex.eqtls=tbl.gtex.eqtls,
                  #tbl.rosmap.eqtls=tbl.rosmap.eqtls
                  )
   invisible(result)

} # build.model
#----------------------------------------------------------------------------------------------------
test_build.model <- function()
{
   message(sprintf("--- test_build.model"))

   if(!snpLocs.initialized) snpLocs.initialized <- initialize.snpLocs()
         #-------------------------------------------------------
         # a 10k region with weak trena, good breaks of those tfs
         #-------------------------------------------------------

   roi.10k <- list(chrom="chr7", start=100048000, end=100058000, width=10001, string="chr7:100,048,000-100,058,000")
   printf("roi.10k width: %dk", as.integer((roi.10k$end - roi.10k$start)/1000, digits=1))

   targetGene <- "PMS2P1"
   x <- build.model(targetGene, "GTEx_V8.Brain_Cortex", roi.10k, eqtl.pval.threshold=1e-3,
                    rna.correlation.threshold=0.2)
   checkEquals(x, list())
   targetGene <- "STAG3"
   x <- build.model(targetGene, "GTEx_V8.Brain_Cortex", roi.10k, eqtl.pval.threshold=1e-3,
                    rna.correlation.threshold=0.2, tbl.oc=data.frame(), known.snps=known.snps)

   
   targetGene <- "ZSCAN21"
   #x <- build.model(targetGene, "GTEx_V8.Brain_Cortex", roi.10k, eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2)
   save(x, file="tmp-test.RData")
   checkTrue(nrow(x$trena) > 6)
   checkTrue(all(c("MLX", "MITF", "NEUROD1") %in% intersect(x$trena$gene, x$tbl.breaks$geneSymbol)))

         #-------------------------------------------------------
         #  120k region with weak trena, good breaks of those tfs
         #-------------------------------------------------------

   roi.120k <- list(chrom="chr7", start=100011184, end=100131671, width=120488, string="chr7:100,011,184-100,131,671")
   printf("roi.120k width: %dk", as.integer((roi.120k$end - roi.120k$start)/1000, digits=1))

   roi.528k <- list(chrom="chr7", start=100076962, end=100605694, width=528733, string="chr7:100,076,962-100,605,694")
   printf("roi.528k width: %dk", as.integer((roi.528k$end - roi.528k$start)/1000, digits=1))

   #x <- build.model(targetGene, "GTEx_V8.Brain_Cerebellum", roi.120k, eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2)
   save(x, file="tmp-test.RData")

   genes <- c("PILRB","STAG3L5P","PILRA","ZCWPW1","MEPCE","PMS2P1","ZSCAN21","NYAP1","STAG3","MBLAC1")
   known.snps <- GRanges()
   for(targetGene in genes[8]){
       models <- list()
       t3 <- system.time({
           for(tissue in gtex.brain.tissues[1:2]){
               #x <- build.model(targetGene, tissue, roi.528k, eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2)
               x <- build.model(targetGene, tissue, roi.10k,
                                eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2,
                                tbl.oc=data.frame(), known.snps=known.snps)
               title <- sprintf("%s-%s", targetGene, tissue)
               known.snps <- x$known.snps
               printf("--- %s %s, known.snps now total %d", targetGene, tissue, length(known.snps))
               models[[title]] <- x
           } # for tissue
       }) # system.time
       printf("%d models of %s in %d seconds", length(models), targetGene, as.integer(t3[["elapsed"]]))

       timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M")))
       filename <- sprintf("%s-models-%s.RData", targetGene, timestamp)
       save(models, file=filename)
       } # for targetGene

} # test_build.model
#----------------------------------------------------------------------------------------------------
test_zscan21_2.tissues <- function()
{
   roi.10k <- list(chrom="chr7", start=100048000, end=100058000, width=10001, string="chr7:100,048,000-100,058,000")
   genes <- "ZSCAN21"
   known.snps <- GRanges()
   for(targetGene in genes){
       models <- list()
       t3 <- system.time({
           for(tissue in gtex.brain.tissues[1:2]){
               #x <- build.model(targetGene, tissue, roi.528k, eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2)
               x <- build.model(targetGene, tissue, roi.10k,
                                eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2,
                                tbl.oc=data.frame(), known.snps=known.snps)
               title <- sprintf("%s-%s", targetGene, tissue)
               known.snps <- x$known.snps
               printf("--- %s %s, known.snps now total %d", targetGene, tissue, length(known.snps))
               models[[title]] <- x
           } # for tissue
       }) # system.time
       printf("%d models of %s in %d seconds", length(models), targetGene, as.integer(t3[["elapsed"]]))

       timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M")))
       filename <- sprintf("%s-models-%s.RData", targetGene, timestamp)
       save(models, file=filename)
       } # for targetGene
   lapply(models, function(model) head(model$trena, n=10))
   lapply(models, function(model) head(model$tbl.breaks, n=10))
   browser()
   xyz <- 99

} # test_zscan21_2.tissues
#----------------------------------------------------------------------------------------------------
run.all <- function()
{
   initialize.snpLocs()
   roi.10k <- list(chrom="chr7", start=100048000, end=100058000, width=10001, string="chr7:100,048,000-100,058,000")
   roi.120k <- list(chrom="chr7", start=100011184, end=100131671, width=120488, string="chr7:100,011,184-100,131,671")
   roi.528k <- list(chrom="chr7", start=100076962, end=100605694, width=528733, string="chr7:100,076,962-100,605,694")

   roi <- roi.528k
   genes <- c("PILRB","STAG3L5P","PILRA","ZCWPW1","MEPCE","PMS2P1","ZSCAN21","NYAP1","STAG3","MBLAC1")
   known.snps <- GRanges()

   for(targetGene in genes[9]){
       models <- list()
       t3 <- system.time({
          for(tissue in gtex.brain.tissues[3]){
             x <- build.model(targetGene, tissue, roi,
                              eqtl.pval.threshold=1e-5,
                              rna.correlation.threshold=0.2,
                              tbl.oc=data.frame(), known.snps=known.snps)
             title <- sprintf("%s-%s", targetGene, tissue)
             known.snps <- x$known.snps
             printf("--- %s %s, known.snps now total %d", targetGene, tissue, length(known.snps))
             models[[title]] <- x
            } # for tissue
        }) # system.time
       printf("%d models of %s in %d seconds", length(models), targetGene, as.integer(t3[["elapsed"]]))
       timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M")))
       filename <- sprintf("%s-models-%s-%s.RData", targetGene, "roi.528k", timestamp)
       save(models, file=filename)
       } # for targetGene

} # run.all
#----------------------------------------------------------------------------------------------------
view.gtex.eqtls <- function()
{

} # veiw.gtex.eqtls
#----------------------------------------------------------------------------------------------------
if(!interactive())
   run.all()
