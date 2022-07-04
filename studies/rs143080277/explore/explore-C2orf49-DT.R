library(TrenaProjectAD)
library(EndophenotypeExplorer)
library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
source("~/github/tms-makeAndBreak/R/tmsMB-class.R")
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
#------------------------------------------------------------------------------------------
tag.snp <- "rs143080277"
tag.snp.chrom <- "chr2"
rsid.oi <- tag.snp
tag.snp.hg38 <- 105749599
shoulder <- 500000
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
gtex.brain.tissues <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
                        "GTEx_V8.Brain_Cerebellum",
                        "GTEx_V8.Brain_Cortex",
                        "GTEx_V8.Brain_Frontal_Cortex_BA9",
                        "GTEx_V8.Brain_Hippocampus",
                        "GTEx_V8.Brain_Hypothalamus")


if(!exists("igv") & FALSE){
    igv <- start.igv("all")
    tbl.track <- data.frame(chrom=tag.snp.chrom,
                            start=tag.snp.hg38-1,
                            end=tag.snp.hg38,
                            name=tag.snp,
                            stringsAsFactors=FALSE)
    track <- DataFrameAnnotationTrack(tag.snp, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    showGenomicRegion(igv, "C2orf49")
    zoomOut(igv)
    }

if(!exists("tbl.fimo.raw")){
   data.dir <- "~/github/tms-makeAndBreak/studies/rs143080277/shared"
   full.path <- file.path(data.dir, "tbl.fimo.NCK2.RData")
   stopifnot(file.exists(full.path))
   tbl.fimo.raw <- get(load(full.path))
   }

if(!exists("tbl.boca")){
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.boca <- get(load(file.path(data.dir, filename)))
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

if(!exists("tbl.eqtl.rosmap")){
   data.dir <- "../shared"
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.eqtl.rosmap.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   dim(tbl.eqtl.rosmap.raw)
   "ECRG4" %in% tbl.eqtl.rosmap.raw$genesymbol
   dim(tbl.eqtl.rosmap.raw)
   tbl.eqtl.rosmap <- subset(tbl.eqtl.rosmap.raw, pvalue < 0.001)
   dim(tbl.eqtl.rosmap)   # 6262 12
   tbl.tmp <- subset(tbl.eqtl.rosmap, genesymbol=="ECRG4")
   if(FALSE){
      track <- GWASTrack("rosmap eQTLs", tbl.tmp, trackHeight=100,
                          chrom.col=1, pos.col=3, pval.col=5)
      displayTrack(igv, track)
      }
   "C2orf40" %in% tbl.eqtl.rosmap.raw$genesymbol
   } # tbl.eqtl.rosmap

if(!exists("tbl.eqtl.gtex")){
   data.dir <- "../shared"
   tbl.eqtl.gtex.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   if(!grepl("chr", tbl.eqtl.gtex.raw$chrom[1]))
      tbl.eqtl.gtex.raw$chrom <- paste0("chr", tbl.eqtl.gtex.raw$chrom)
   dim(tbl.eqtl.gtex.raw)
   targetGene <- "ZSCAN21"
   tbl.track <- subset(tbl.eqtl.gtex.raw, pvalue < 0.001 & gene==targetGene)
   dim(tbl.track)
   if(FALSE){
      track <- GWASTrack("GTEx eQTLs", tbl.track, trackHeight=100,
                          chrom.col=7, pos.col=8, pval.col=2)
      displayTrack(igv, track)
      }
   } # tbl.eqtl.gtex


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
    tbl.fimo <- get(load(file.path("../shared", "tbl.fimo.NCK2.RData")))
    gr.fimo <- GRanges(tbl.fimo)
    }

if(!exists("tbl.eqtl.c2orf49dt")){
    tbl.eqtl.c2orf49dt <- subset(tbl.eqtl.rosmap, genesymbol=="C2orf49-DT")
    tbl.eqtl.c2orf49dt$score <- with(tbl.eqtl.c2orf49dt, -log10(pvalue) * beta)
    tbl.track <- tbl.eqtl.c2orf49dt[, c("chrom", "hg38", "hg38", "score")]
    colnames(tbl.track) <- c("chrom", "start", "end", "score")
    tbl.track$start <- tbl.track$start - 1
    if(FALSE){
       track <- DataFrameQuantitativeTrack("rosmap eqtl", tbl.track, autoscale=FALSE, color="black", min=-7, max=7)
       displayTrack(igv, track)
       }
    }

if(!exists("tbl.eqtl.c2orf49dt.gtex")){
    tbl.eqtl.c2orf49dt.gtex <- subset(tbl.eqtl.gtex, gene=="RP11-332H14.2")
    dim(tbl.eqtl.c2orf49dt.gtex)
    tbl.eqtl.c2orf49dt.gtex$score <- with(tbl.eqtl.c2orf49dt.gtex, -log10(pvalue) * beta)
    tbl.track <- tbl.eqtl.c2orf49dt.gtex[, c("chrom", "hg38", "hg38", "score")]
    colnames(tbl.track) <- c("chrom", "start", "end", "score")
    tbl.track$start <- tbl.track$start - 1
    if(FALSE){
       track <- DataFrameQuantitativeTrack("gtex eqtl", tbl.track, autoscale=FALSE, color="black", min=-7, max=7)
       displayTrack(igv, track)
       }
    }




# see which genes have strong pval eqtls in the region covered by rs143080277
#------------------------------------------------------------------------------------------
# restrict region to 82kb region chr7:100,333,519-100,416,094, where LD r^2 >= 0.5, and gwascat.ad
# provides support.  this is all downstream of the tag snp
as.data.frame(sort(table(subset(tbl.eqtl.gtex.raw, chrom==tag.snp.chrom &
                                                    hg38 >= tag.snp.hg38 - shoulder &
                                                    hg38 <= tag.snp.hg38 + shoulder &
                                                    pvalue < 1e-3)$gene), decreasing=TRUE))
#             Var1 Freq
# 1        C2orf40  179
# 2     AC012360.4   79
# 3  RP11-332H14.2   41
# 4           FHL2   13
# 5           NCK2    5
# 6          GPR45    2
# 7      LINC01158    1
# 8          MRPS9    1
# 9         POU3F3    1
# 10 RP11-332H14.1    1
as.data.frame(sort(table(subset(tbl.eqtl.rosmap, chrom==tag.snp.chrom &
                                                  hg38 >= tag.snp.hg38 - shoulder &
                                                  hg38 <= tag.snp.hg38 + shoulder &
                                                  pvalue < 1e-3)$genesymbol), decreasing=TRUE))
#         Var1 Freq
# 1 C2orf49-DT.GTEX  159
# 2      ECRG4   88
# 3       NCK2   53
# 4       UXS1   17
# 5       FHL2    8
# 6   TGFBRAP1    5
# 7      MRPS9    2


#snpLocs.initialized <- FALSE
#----------------------------------------------------------------------------------------------------
initialize.snpLocs <- function()   # load into memory just once, saving time in motifbreakR prep
{
    t0 <- system.time(x0 <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, tag.snp))
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

    track <- GWASTrack("rosmap eQTLs", tbl.eqtl.rosmap, trackHeight=100,
                       chrom.col=1, pos.col=3, pval.col=5)
    displayTrack(igv, track)

} # showTracks
#----------------------------------------------------------------------------------------------------
identify.affected.genes <- function()
{
   roi <- getGenomicRegion(igv)  # 45k between 100,312,189  and  100,356,932

   fivenum(tbl.eqtl.rosmap$pvalue)
   tbl.sub <-
     subset(tbl.eqtl.rosmap,
            hg38 >= roi$start & hg38 <= roi$end &
            pvalue < 1e-5)
   goi <- unique(tbl.sub$genesymbol)
   length(goi)  # 13

   tbl.freq <- as.data.frame(sort(table(tbl.sub$genesymbol), decreasing=TRUE))
   colnames(tbl.freq) <- c("gene", "eqtl.count")

   scoreGene <- function(gene){
     #tbl.gene <- subset(tbl.eqtl.rosmap, genesymbol==gene)
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
build.model <- function(targetGene, chromosome=NA, gtex.tissue, roi,
                        eqtl.pval.threshold, rna.correlation.threshold,
                        tbl.oc, known.snps)

{
   printf("--- building model of %s in %s over %dk",
          targetGene, gtex.tissue, as.integer(roi$width/1000))

   trenaProject <- TrenaProjectAD()
   tbl.fimo <- subset(tbl.fimo.raw, chrom==roi$chrom & start >= roi$start & end <= roi$end)

   tbl.eqtl.gtex <- subset(tbl.eqtl.gtex.raw,
                            pvalue <= eqtl.pval.threshold &
                            gene==targetGene &
                            hg38 >= roi$start & hg38 <= roi$end &
                            id == gtex.tissue)
   message(sprintf("--- %d %s eqtls for %s in region", nrow(tbl.eqtl.gtex), gtex.tissue, targetGene))

   if(nrow(tbl.eqtl.gtex) == 0){
       message(sprintf("--- no gtex eqtls for %s in tissue %s", targetGene, gtex.tissue))
       return(list())
       }

   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.eqtl.rosmap <- subset(tbl.eqtl.rosmap.raw,
                             pvalue < eqtl.pval.threshold &
                             genesymbol==targetGene &
                             hg38 >= roi$start & hg38 <= roi$end)
   if(nrow(tbl.eqtl.rosmap) == 0){
       message(sprintf("--- no rosmap eqtls for %s", targetGene))
       return(list())
       }

   message(sprintf("--- %d rosmap eqtls for %s in region", nrow(tbl.eqtl.rosmap), targetGene))
   message(sprintf("--- eQTL rsids gtex: %d   rosmap: %d   shared: %d",
                   length(tbl.eqtl.gtex$rsid),
                   length(tbl.eqtl.rosmap$rsid),
                   length(intersect(tbl.eqtl.gtex$rsid, tbl.eqtl.rosmap$rsid))))

   tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.eqtl.gtex,
                    tbl.eqtl.rosmap, tbl.oc, known.snps, chromosome=chromosome)

   trenaProject <- TrenaProjectAD()
   tms$setStudyRegion(chrom=roi$chrom, start=roi$start, end=roi$end)

   study.region <- tms$getStudyRegion()
   checkTrue(gtex.tissue %in% gtex.brain.tissues)
   tms$set.current.GTEx.eqtl.tissue(gtex.tissue)

   tms$run.tms()
   tms$add.eqtls.toTmsTable()

   tbl.tms <- tms$get.tmsTable()
   tbl.tms <- subset(tbl.tms, fimo_pvalue <= 1e-4)
   ampad.and.gtex.eqtls <- with(tbl.tms,
                                ampad.eqtl.pval < 1 & gtex.eqtl.pval < eqtl.pval.threshold)
   ampad.or.gtex.eqtls <- with(tbl.tms,
                                ampad.eqtl.pval < 1 | (gtex.eqtl.pval < eqtl.pval.threshold))
   print(table(ampad.and.gtex.eqtls))
   print(table(ampad.or.gtex.eqtls))
   strong.ampad.eqtl <- tbl.tms$ampad.eqtl.score > 1
   print(table(strong.ampad.eqtl))
   openChromatin <- tbl.tms$oc
   print(table(openChromatin))
   chip <- tbl.tms$chip
   if(is.null(chip)) chip <- rep(FALSE, length(openChromatin))
   print(table(chip))
   correlated.expression <- with(tbl.tms, abs(cor.all) >= rna.correlation.threshold)
   genehancer <- tbl.tms$gh > 1
   table(genehancer)

   filter <- correlated.expression &
       (ampad.or.gtex.eqtls | (openChromatin & strong.ampad.eqtl & genehancer))
   print(table(filter))

       #subset(tbl.tms, (ampad.eqtl.pval < 1 & gtex.eqtl.pval < eqtl.pval.threshold) &
       #                abs(cor.all) >= rna.correlation.threshold)

   tbl.tms.filtered <- tbl.tms[which(filter),]
   print(dim(tbl.tms.filtered))

   tms$set.tmsFilteredTable(tbl.tms.filtered)
   tf.candidates <- unique(tbl.tms.filtered$tf)
   printf("tf.candidates: %d", length(tf.candidates))

   tbl.trena <- data.frame()
   tbl.breaks <- data.frame()

   if(length(tf.candidates) > 0){
      tms$run.trena(tf.candidates)
      tbl.trena <- tms$get.trenaTable()
      head(tbl.trena, n=10)
      tms$breakMotifs(head(tbl.trena, n=12), tbl.tms.filtered, tbl.eqtl.rosmap)
      tbl.breaks <- tms$get.breaksTable()
      if(is.null(tbl.breaks)) tbl.breaks <- data.frame()
      if(!is.null(tbl.breaks)){
          if(nrow(tbl.breaks) > 0){
             # we are (mostly) interested in the strongest breaks, high negative pctDelta
             new.order <- order(tbl.breaks$pctDelta, decreasing=FALSE)
             tbl.breaks <- tbl.breaks[new.order,]
          } # nrow tbl.breaks
          } # is.null
      } # tf.candidates

   new.known.snps <- tms$get.knownSnps()
   result <- list(#tms=tbl.tms,
                  tms.filtered=tbl.tms.filtered,
                  trena=head(tbl.trena, n=20),
                  tf.candidates=tf.candidates,
                  known.snps=new.known.snps,
                  tbl.breaks=tbl.breaks,
                  breaks=tms$get.motifBreaks(),
                  roi=roi,
                  eqtl.threshold=eqtl.pval.threshold,
                  rna.correlation.threshold=rna.correlation.threshold
                  #tbl.eqtl.gtex=tbl.eqtl.gtex,
                  #tbl.rosmap.eqtls=tbl.rosmap.eqtls
                  )
   invisible(result)

} # build.model
#----------------------------------------------------------------------------------------------------
# C2orf49-DT, sometimes knonw as RP11-332H14.2, or ENSG00000272994
# this function tests the use of modified GTEx expression and eQTL tables to use the
# name favored by amp-ad.
test_c2orf49dt <- function()
{

   message(sprintf("--- test_c2orf49dt"))
   old.name <- "RP11-332H14.2"   # used by gtex
   new.name <- "C2orf49-DT"      # used by genehancer and rosmap
   ensg <-  "ENSG00000272994"       #

   targetGene <- new.name
   checkGeneNameFoundInAllSources(targetGene, chromosome="chr2")

   roi.cherry <- list(chrom="chr2", start=105243579, end=105410755)
   roi.cherry$width <- with(roi.cherry, 1+end-start)
   roi.cherry$width  # 167k

   roi.224k <- list(chrom="chr2", start=105236169, end=105460720)
   roi.224k$width <- with(roi.224k, 1+end-start)
   roi.224k$width  # 224,552

   tissue <- "GTEx_V8.Brain_Cerebellar_Hemisphere"

   known.snps <- GRanges()

   x <- build.model(targetGene, tissue, roi.224k,
                    eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps, chromosome="chr2")

   known.snps <- x$known.snps  # keep these for running again

} # test_c2orf49dt
#----------------------------------------------------------------------------------------------------
# ECRG4: Augurin, neuropeptide hormone activity, cellular senescence, central nervous system development
# used to be called "C2orf40".   make sure the modern name can be used throughout
test_c2orf40 <- function()
{
   message(sprintf("--- test_c2orf40"))

   old <- "C2orf40"         # older alias
   new <- "ECRG4"           # preferred, new HUGO

   checkGeneNameFoundInAllSources(new)

   targetGene <- new

   roi.229k <- list(chrom="chr2", start=106055752, end=106285069)
   roi.229k$width <- with(roi.229k, 1+end-start)
   roi.229k$width  # 167k

   tissue <- "GTEx_V8.Brain_Cerebellar_Hemisphere"

   known.snps <- GRanges()

   targetGene <- new
   x <- build.model(targetGene, tissue, roi.229k,
                    eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps)

   known.snps <- x$known.snps  # keep these for running again

} # test_c2orf40
#----------------------------------------------------------------------------------------------------
combine.tables <- function(tbl.trena, tbl.tms, tbl.breaks, targetGene, TF)
{
    tbl.trena <- subset(tbl.trena, gene==TF)[, c("gene", "betaLasso", "spearmanCoeff", "rfScore", "rfNorm", "tfbs", "tf.rank")]
    colnames(tbl.trena) <- c("tf", "betaLasso", "spearman", "rfScore", "rfNorm", "tfbs", "tf.rank")
    tbl.tms <- subset(tbl.tms, tf==TF)  # only some of these will be broken
    tbl.breaks <- subset(tbl.breaks, geneSymbol==TF)
    colnames(tbl.breaks)[3:6] <- c("hg38", "rsid", "tf", "motif")
    tbl.breaks <- tbl.breaks[, c("chrom", "hg38", "rsid", "tf", "motif", "pctRef", "pctDelta")]
    #printf("nrow tbl.breaks for %s -> %s: %d", TF, targetGene, nrow(tbl.breaks))
         #--------------------------------------------------
         # intersect breaks with tfbs in tbl.tms
         #--------------------------------------------------
    if(nrow(tbl.breaks) == 0){

        return(list(trenaScore=0, tbl.tbt=data.frame()))
        }

    gr.breaks <- with(tbl.breaks, GRanges(seqname=chrom[1], IRanges(hg38)))
    gr.tms    <- GRanges(tbl.tms)
    tbl.ov <- as.data.frame(findOverlaps(gr.breaks, gr.tms))
    break.indices <- tbl.ov[, 1]
    tms.indices <- tbl.ov[, 2]
    tbl.bt <- cbind(tbl.breaks[break.indices,], tbl.tms[tms.indices,])
    dup.cols <- c(grep("chrom", colnames(tbl.bt))[2],
                      grep("tf", colnames(tbl.bt))[2])
    if(length(dup.cols) > 0)
       tbl.bt <- tbl.bt[, -(dup.cols)]

    tbl.tbt <- merge(tbl.trena, tbl.bt, by="tf")
    motifBreak.score <- with(tbl.tbt, abs(pctDelta) * 10)
    eqtl.score <- with(tbl.tbt, -log10(ampad.eqtl.pval)* abs(ampad.eqtl.beta) * 1)
    trena.score <- with(tbl.tbt, (abs(betaLasso) * 10) + (rfNorm * 10))
    tbl.tbt$breakScore <- motifBreak.score
    tbl.tbt$eqtlScore <- eqtl.score
    tbl.tbt$trenaScore <- trena.score
    tbl.tbt$compositeScore <- as.integer(0.5 + (motifBreak.score * eqtl.score * trena.score))

    coi <- c("targetGene", "tf", "compositeScore", "breakScore", "eqtlScore", "trenaScore","tf.rank", "rsid", "chrom", "hg38",
             "start", "end", "pctRef","pctDelta",
             "tfbs","motif", "betaLasso","spearman","rfScore","rfNorm",
             "fimo_pvalue",
             "phast7", "phast100","gh","oc","tss","motif_id","ampad.eqtl.pval",
             "ampad.eqtl.beta","ampad.eqtl.score", "gtex.eqtl.pval","gtex.eqtl.beta","gtex.eqtl.score")
    setdiff(coi, colnames(tbl.tbt))

    tbl.tbt <- tbl.tbt[, coi]
    trenaScore <- as.integer(0.5 + mean(tbl.tbt$compositeScore))

    return(list(trenaScore=trenaScore, tbl.tbt=tbl.tbt))

} # combine.tables
#----------------------------------------------------------------------------------------------------
test_combine.tables <- function()
{
   load("PILRB-models-roi.528k-2022.jun.24-13:02.RData")
   targetGene <- "PILRB"
   model <- models[[1]]
   tbl.trena <- head(model$trena, n=10)
   tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
   tbl.breaks <- model$tbl.breaks
   tbl.tms <- model$tms.filtered
   tbl.tms$targetGene <- targetGene

   subset(tbl.trena, gene %in% tbl.breaks$geneSymbol)
   scores <- unlist(lapply(tbl.trena$gene, function(tf)
       combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF=tf)$trenaScore))
   tbl.trena$breakageScore <- scores
   breakage.total <- sum(tbl.trena$breakageScore)
   checkEqualsNumeric(breakage.total, 7755, tol=500)

   unbroken <- setdiff(tbl.trena$gene, tbl.breaks$geneSymbol)
       # these should all have a zero breakage score
   checkEquals(sum(subset(tbl.trena, gene %in% unbroken)$breakageScore), 0)

       # now examine the calculation, step by step, only TEAD1, the top-ranked TF
   tbl.tbt <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF="TEAD1")$tbl.tbt
   tbt <- as.list(tbl.tbt[1,])

      #  motifBreak.score <- with(tbl.tbt, abs(pctDelta) * 10)
   checkEqualsNumeric(tbt$pctDelta, -0.09885357, tolerance=1e-7)
   checkEqualsNumeric(tbt$breakScore, 0.9885357, tolerance=1e-6)  # abs(pctDelta) * 10

      # eqtl.score <- with(tbl.tbt, -log10(ampad.eqtl.pval)* abs(ampad.eqtl.beta) * 1)
    checkEqualsNumeric(tbt$ampad.eqtl.pval, 1.366627e-116, tol=1e-100)
    checkEqualsNumeric(tbt$ampad.eqtl.beta, 1.425089, tol=1e-5)
    checkEqualsNumeric(tbt$eqtlScore, 165.117, tol=1e-2)

      # trena.score <- with(tbl.tbt, (abs(betaLasso) * 10) + (rfNorm * 10))
    checkEqualsNumeric(tbt$betaLasso, -0.30, tol=0.01)
    checkEquals(tbt$rfNorm, 1)
    checkEquals(tbt$trenaScore, 13)

    compositeScore <- round(0.5 + (165.117 * 13 * 0.9885 ))
    checkEqualsNumeric(compositeScore, 2122, tol=0.05)

    tbl.tbt$breakScore <- motifBreak.score
    tbl.tbt$eqtlScore <- eqtl.score
    tbl.tbt$trenaScore <- trena.score
    tbl.tbt$compositeScore <- as.integer(0.5 + (motifBreak.score * eqtl.score * trena.score))

   # as.data.frame(t(tbl.tbt[1,]))
        # targetGene                                    PILRB
        # tf                                            TEAD1
        # compositeScore                                 2122
        # breakScore                                0.9885357
        # eqtlScore                                   165.117
        # trenaScore                                       13
        # tf.rank                                           1
        # rsid                                      rs3087502
        # chrom                                          chr7
        # hg38                                      100357741
        # start                                     100357734
        # end                                       100357743
        # pctRef                                    0.9015386
        # pctDelta                                -0.09885357
        # tfbs                                              1
        # motif                                      MA0090.3
        # betaLasso                                      -0.3
        # spearman                                      -0.56
        # rfScore                                       26.21
        # rfNorm                                            1
        # fimo_pvalue                                1.54e-05
        # phast7                                            0
        # phast100                                          0
        # gh                                           258.32
        # oc                                            FALSE
        # tss                                            -281
        # motif_id         Hsapiens-jaspar2018-TEAD1-MA0090.2
        # ampad.eqtl.pval                       1.366627e-116
        # ampad.eqtl.beta                            1.425089
        # ampad.eqtl.score                            165.117
        # gtex.eqtl.pval                          3.66317e-43
        # gtex.eqtl.beta                             0.850749
        # gtex.eqtl.score                            36.10251


} # test_combine.tables
#----------------------------------------------------------------------------------------------------
summarize.models <- function()
{

   models.all <- c(PVRIG="PVRIG-models-roi.528k-2022.jun.24-13:23.RData",
                   MEPCE="MEPCE-models-roi.528k-2022.jun.24-13:22.RData",
                   STAG3="STAG3-models-roi.528k-2022.jun.24-13:20.RData",
                   ZCWPW1="ZCWPW1-models-roi.528k-2022.jun.24-13:17.RData",
                   PILRA="PILRA-models-roi.528k-2022.jun.24-13:12.RData",
                   STAG3L5P="STAG3L5P-models-roi.528k-2022.jun.24-13:07.RData",
                   PILRB="PILRB-models-roi.528k-2022.jun.24-13:02.RData")

   tbls <- list()
   i <- 0
   for(targetGene in (names(models.all))){
      models <- get(load(models.all[[targetGene]]))
      model.names <- names(models)
      for(model.name in model.names){
          i <- i + 1
          tbl <- data.frame(targetGene=targetGene, tissue=model.name, score=0, stringsAsFactors=FALSE)
          model <- models[[model.name]]
          if(length(model) >= 0){
              tbl.trena <- model$trena
              if(!is.null(tbl.trena)){
                  tbl.trena <- head(tbl.trena, n=10)
                  printf("targetGene %s model %s", targetGene, model.name)
                  tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
                  tbl.breaks <- model$tbl.breaks
                  if(nrow(tbl.breaks) > 0){
                      tbl.tms <- model$tms.filtered
                      tbl.tms$targetGene <- targetGene
                      scores <- unlist(lapply(tbl.trena$gene, function(tf)
                          combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF=tf)$trenaScore))
                      tbl.trena$breakageScore <- scores
                      breakage.total <- sum(tbl.trena$breakageScore)
                      printf("%s %s: %d", targetGene, model.name, breakage.total)
                      tbl$score <- breakage.total
                     } # breaks found
                  } # !null tbl.trena
              } # model found
          tbls[[i]] <- tbl
          } # for model.name
      } # for targetGene

    tbl <- do.call(rbind, tbls)
    tbl$tissue <- sub(".*Brain_", "", tbl$tissue)
    tissues <- sort(unique(tbl$tissue))
    genes <- sort(unique(tbl$targetGene))
    mtx <- matrix(rep(0,42), nrow=length(genes), dimnames=list(genes, tissues))
    for(r in seq_len(nrow(tbl))){
       gene <- tbl$targetGene[r]
       tissue <- tbl$tissue[r]
       score <- tbl$score[r]
       mtx[gene, tissue] <- score
       }

   tbl.summary <- as.data.frame(mtx)
   tbl.summary$sum <- as.integer(rowSums(tbl.summary))
   new.order <- order(tbl.summary$sum, decreasing=TRUE)
   tbl.summary[new.order,]

} # summarize.models
#----------------------------------------------------------------------------------------------------
run.all <- function()
{
   initialize.snpLocs()

   roi.1m <- list(chrom="chr2",  start=105224001, end=106224000)
   roi.1m$width <- with(roi.1m, 1 + end - start)
   roi.1m
   roi <- roi.1m

        # as.data.frame(sort(table(subset(tbl.eqtl.rosmap,
        #                                 chrom==tag.snp.chrom &
        #                                 hg38 >= tag.snp.hg38 - shoulder &
        #                                 hg38 <= tag.snp.hg38 + shoulder &
        #                                 pvalue < 1e-3)$genesymbol), decreasing=TRUE))

    #          Var1 Freq
    #  1 C2orf49-DT  159    # name not recognized by all parts.  needs work
    #  2      ECRG4   88    # gtex and ampad prefer C2orf40
    #  3       NCK2   53
    #  4       UXS1   17
    #  5       FHL2    8
    #  6   TGFBRAP1    5
    #  7      MRPS9    2

     # c("C2orf49-DT", "ECRG4",

   genes <- c("C2orf49-DT", "ECRG4", "NCK2", "UXS1", "FHL2", "TGFBRAP1", "MRPS9")
   length(genes)
   genes <- intersect(genes, tbl.eqtl.gtex.raw$gene)
   length(genes)

   known.snps <- GRanges()


   for(targetGene in genes){
       models <- list()
       t3 <- system.time({
          for(tissue in gtex.brain.tissues){
             x <- build.model(targetGene, tissue, roi,
                              eqtl.pval.threshold=1e-5,
                              rna.correlation.threshold=0.2,
                              tbl.oc=tbl.boca, known.snps=known.snps)
             title <- sprintf("%s-%s", targetGene, tissue)
             if("known.snps" %in% names(x))
                known.snps <- x$known.snps  # so they may accumulate
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
