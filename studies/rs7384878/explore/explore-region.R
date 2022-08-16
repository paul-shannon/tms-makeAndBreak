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
# see which genes have strong pval eqtls in the region covered by rs7384878
#------------------------------------------------------------------------------------------
# restrict region to 82kb region chr7:100,333,519-100,416,094, where LD r^2 >= 0.5, and gwascat.ad
# provides support.  this is all downstream of the tag snp
#  tbl.xtab <- table(subset(tbl.gtex.eqtls.raw, chrom=="chr7" &
#                                               hg38 >= 100333519 & hg38 <= 100416094 &
#                                               pvalue < 1e-5)$gene)
#        Var1 Freq
# 1     PILRB  476
# 2  STAG3L5P  464
# 3     PILRA  458
# 4    ZCWPW1  330
# 5     STAG3  108
# 6     MEPCE   60
# 7     PVRIG   57
# 8    PMS2P1   36
# 9     CNPY4   26
# 10  C7orf61   19
# 11     GPC2   16
# 12  LAMTOR4   11
# 13    AZGP1    9
# 14    SAP25    9
# 15    AGFG2    8
# 16   MBLAC1    7
# 17  ZSCAN21    7
# 18  GAL3ST4    2
# 19     MCM7    2
# 20      EPO    1
# 21     GATS    1
#------------------------------------------------------------------------------------------
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

if(!exists("tbl.boca")){
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.boca <- get(load(file.path(data.dir, filename)))
   }

if(!exists("tbl.haploreg")){
    data.dir <- "../shared"
    haploreg.file <- "haploreg-rs7384878.tsv"
    tbl.haploreg <-
       read.table(file.path(data.dir, haploreg.file), sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
    tbl.haploreg$score <- 0.9999999999 - tbl.haploreg$Rsquared
    if(FALSE){
       tbl.track <- tbl.haploreg[, c("chrom", "hg38", "hg38", "Rsquared")]
       colnames(tbl.track) <- c("chrom", "start", "end", "score")
       tbl.track$start <- tbl.track$start - 1
       tbl.track$chrom <- as.character(tbl.track$chrom)
       track <- DataFrameQuantitativeTrack("LD", tbl.track, color="red", trackHeight=25, autoscale=TRUE)
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
   goi <- "PILRB"
   data.dir <- "../shared"
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.rosmap.nearby.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   dim(tbl.rosmap.nearby.raw)
   tbl.rosmap.nearby.raw$score <- -log10(tbl.rosmap.nearby.raw$pvalue)
   tbl.rosmap.nearby <- subset(tbl.rosmap.nearby.raw, pvalue < 0.001 & genesymbol==goi)
   dim(tbl.rosmap.nearby)   # 6262 12
   if(FALSE){
      #track <- GWASTrack("rosmap eQTLs", tbl.rosmap.nearby.raw, trackHeight=100,
       #                    chrom.col=1, pos.col=3, pval.col=5)
      coi <- c("chrom", "hg38", "hg38", "score", "genesymbol")
      tbl.track <- tbl.rosmap.nearby.raw[, coi]
      colnames(tbl.track) <- c("chrom", "start", "end", "score")
      tbl.track$start <- tbl.track$start - 1
      dim(tbl.track)
      track <- DataFrameQuantitativeTrack("rosmap eqtl", tbl.track, autoscale=TRUE, min=0, max=4)
      displayTrack(igv, track)
         # now the top affected genes
      target.genes <- names(sort(table(subset(tbl.track, score > 50)$genesymbol), decreasing=TRUE))
      for(gene in target.genes){
         tbl.track.gene <- subset(tbl.track, genesymbol==gene)
         title <- sprintf("%s eqtl", gene)
         track <- DataFrameQuantitativeTrack(title, tbl.track.gene, autoscale=TRUE, min=0, max=4, color="random")
         displayTrack(igv, track)
         } # for gene
      } # FALSE
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
    deleters <- grep("HOCOMOCO", tbl.fimo$motif_id)
    if(length(deleters) > 0)
       tbl.fimo <- tbl.fimo[-deleters,]
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
   tbl.fimo.raw <- get(load(file.path(data.dir, "tbl.fimo.PILRA.RData")))
   tbl.fimo <- subset(tbl.fimo.raw, chrom==roi$chrom & start >= roi$start & end <= roi$end)

   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   if(!grepl("chr", tbl.gtex.eqtls.raw$chrom[1]))
      tbl.gtex.eqtls.raw$chrom <- paste0("chr", tbl.gtex.eqtls.raw$chrom)

   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw,
                            pvalue <= eqtl.pval.threshold &
                            gene==targetGene &
                            hg38 >= roi$start & hg38 <= roi$end &
                            id == gtex.tissue)
   message(sprintf("--- %d %s eqtls for %s in region", nrow(tbl.gtex.eqtls), gtex.tissue, targetGene))

   if(nrow(tbl.gtex.eqtls) == 0){
       message(sprintf("--- no gtex eqtls for %s in tissue %s", targetGene, gtex.tissue))
       return(list())
       }

   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.rosmap.eqtls.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   tbl.rosmap.eqtls <- subset(tbl.rosmap.eqtls.raw,
                              pvalue < eqtl.pval.threshold &
                              genesymbol==targetGene &
                              hg38 >= roi$start & hg38 <= roi$end)
   if(nrow(tbl.rosmap.eqtls) == 0){
       message(sprintf("--- no rosmap eqtls for %ss", targetGene))
       return(list())
       }

   message(sprintf("--- %d rosmap eqtls for %s in region", nrow(tbl.rosmap.eqtls), targetGene))
   message(sprintf("--- eQTL rsids gtex: %d   rosmap: %d   shared: %d",
                   length(tbl.gtex.eqtls$rsid),
                   length(tbl.rosmap.eqtls$rsid),
                   length(intersect(tbl.gtex.eqtls$rsid, tbl.rosmap.eqtls$rsid))))
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
   tbl.tms <- subset(tbl.tms, fimo_pvalue <= 1e-4)
   ampad.and.gtex.eqtls <- with(tbl.tms,
                                ampad.eqtl.pval < 1 & gtex.eqtl.pval < eqtl.pval.threshold)
   table(ampad.and.gtex.eqtls)
   strong.ampad.eqtl <- tbl.tms$ampad.eqtl.score > 1
   openChromatin <- tbl.tms$oc
   chip <- tbl.tms$chip
   correlated.expression <- with(tbl.tms, abs(cor.all) >= rna.correlation.threshold)
   genehancer <- tbl.tms$gh > 1


   filter <- correlated.expression &
       (ampad.and.gtex.eqtls | (openChromatin & chip & strong.ampad.eqtl & genehancer))
   table(filter)

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
      tms$breakMotifs(head(tbl.trena, n=12), tbl.tms.filtered, tbl.rosmap.eqtls)
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
   #browser()
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
                  #tbl.gtex.eqtls=tbl.gtex.eqtls,
                  #tbl.rosmap.eqtls=tbl.rosmap.eqtls
                  )
   invisible(result)

} # build.model
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
# bud31 is 931 kb from the rs7384878, has few rosmap eqtls nearby, fails to model due to
# no gtex eqtls in region & tissue,  below threshold, then with no candidate tfs
# should model w/breakage nearly null
test_bud31 <- function()
{
   roi.82k  <- list(chrom="chr7", start=100333519, end=100416094, width=82575)
   roi.25k  <- list(chrom="chr7", start=99394115, end=99419550, width=25436, string="chr7:99,394,115-99,419,550")
   roi.120k <- list(chrom="chr7", start=100011184, end=100131671, width=120488, string="chr7:100,011,184-100,131,671")
   roi.528k <- list(chrom="chr7", start=100076962, end=100605694, width=528733, string="chr7:100,076,962-100,605,694")

   gene <- "BUD31"
   tissue <- "GTEx_V8.Brain_Cortex"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not BUD31 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.528k,
                    eqtl.pval.threshold=1e-3, rna.correlation.threshold=0.2,
                    tbl.oc=data.frame(), known.snps=known.snps)
   checkEquals(x, list())

      #------------------------------------------------------------
      # relax the eqtl pval threshold, but no tf candidates found
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.528k,
                    eqtl.pval.threshold=5e-3, rna.correlation.threshold=0.2,
                    tbl.oc=data.frame(), known.snps=known.snps)
   checkEquals(x$trena, data.frame())

   x <- build.model(gene, tissue, roi.528k,
                    eqtl.pval.threshold=5e-3, rna.correlation.threshold=0.0,
                    tbl.oc=data.frame(), known.snps=known.snps)
   checkEquals(x$trena, data.frame())

} # test_bud31
#----------------------------------------------------------------------------------------------------
test_stag3l5p <- function()
{
   roi.82k  <- list(chrom="chr7", start=100333519, end=100416094, width=82575)

   gene <- "STAG3L5P"
   tissue <- "GTEx_V8.Brain_Cortex"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not BUD31 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.82k,
                    eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps)
   checkTrue(is.list(x))
   checkTrue(all(c("TGIF2", "TCF7L2","MEF2D","TEAD2") %in%
                 intersect(head(x$trena$gene, n=8), x$tbl.breaks$geneSymbol)))

} # test_stag3l
#----------------------------------------------------------------------------------------------------
test_stag3l5p <- function()
{
   roi.82k  <- list(chrom="chr7", start=100333519, end=100416094, width=82575)

   gene <- "STAG3L5P"
   tissue <- "GTEx_V8.Brain_Cortex"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not BUD31 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.82k,
                    eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps)
   checkTrue(is.list(x))
   checkTrue(all(c("TGIF2", "TCF7L2","MEF2D","TEAD2") %in%
                 intersect(head(x$trena$gene, n=8), x$tbl.breaks$geneSymbol)))

} # test_stag3l
#----------------------------------------------------------------------------------------------------
# bud31 is 931 kb from the rs7384878, has few rosmap eqtls nearby, fails to model due to
# no gtex eqtls in region & tissue,  below threshold, then with no candidate tfs
# should model w/breakage nearly null
test_pilrb <- function()
{
   roi.82k  <- list(chrom="chr7", start=100333519, end=100416094, width=82575)
   roi.25k  <- list(chrom="chr7", start=99394115, end=99419550, width=25436, string="chr7:99,394,115-99,419,550")
   roi.120k <- list(chrom="chr7", start=100011184, end=100131671, width=120488, string="chr7:100,011,184-100,131,671")
   roi.528k <- list(chrom="chr7", start=100076962, end=100605694, width=528733, string="chr7:100,076,962-100,605,694")

   gene <- "PILRB"
   tissue <- "GTEx_V8.Brain_Cerebellar_Hemisphere"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not BUD31 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.82k,
                    eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps)

} # test_pil4b
#----------------------------------------------------------------------------------------------------
# spdye3 is 13 kb from the rs7384878, has quite a few rosmap eqtls nearby
#
test_spdye3 <- function()
{
   roi.115k <- list(chrom="chr7", start=99394115, end=99419550, width=25436,
                    string="chr7:99,394,115-99,419,550")

   gene <- "SPDYE3"
   tissue <- "GTEx_V8.Brain_Cerebellar_Hemisphere"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not SPDYE3 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.115k,
                    eqtl.pval.threshold=0.01, rna.correlation.threshold=0.2,
                    tbl.oc=data.frame(), known.snps=known.snps)
   checkEquals(x, list())

      #------------------------------------------------------------
      # relax the eqtl pval threshold, but no tf candidates found
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.528k,
                    eqtl.pval.threshold=5e-3, rna.correlation.threshold=0.2,
                    tbl.oc=data.frame(), known.snps=known.snps)
   checkTrue(is.list(x))
   checkEquals(x$trena, data.frame())

} # test_spdye3
#----------------------------------------------------------------------------------------------------
# only 3 strong (1e-10) gtex eqtls in haploreg region for zscan21
test_zscan21 <- function()
{
   message(sprintf("--- test_zscan21"))
   roi.454k <- list(chrom="chr7", start=100127439, end=100582263, width=454825,
                    string="chr7:100,127,439-100,582,263")
   gene <- "ZSCAN21"
   tissue <- "GTEx_V8.Brain_Cerebellar_Hemisphere"

   known.snps <- GRanges()

      #------------------------------------------------------------
      # will fail since there are not SPDYE3 gtex eqtls pval <= 1e-3
      #------------------------------------------------------------

   x <- build.model(gene, tissue, roi.454k,
                    eqtl.pval.threshold=1e-5, rna.correlation.threshold=0.2,
                    tbl.oc=tbl.boca, known.snps=known.snps)
   checkEquals(x, list())


} # test_zscan21
#----------------------------------------------------------------------------------------------------
run.all <- function()
{
   initialize.snpLocs()
   roi.10k <- list(chrom="chr7", start=100048000, end=100058000, width=10001, string="chr7:100,048,000-100,058,000")
   roi.120k <- list(chrom="chr7", start=100011184, end=100131671, width=120488, string="chr7:100,011,184-100,131,671")
   roi.528k <- list(chrom="chr7", start=100076962, end=100605694, width=528733, string="chr7:100,076,962-100,605,694")
   roi.82k  <- list(chrom="chr7", start=100333519, end=100416094, width=82575)

   roi <- roi.82k
     # 1     PILRB  476
     # 2  STAG3L5P  464
     # 3     PILRA  458
     # 4    ZCWPW1  330
     # 5     STAG3  108
     # 6     MEPCE   60
     # 7     PVRIG   57

   #genes <- c("STAG3", "MEPCE", "PVRIG")
   genes <- c("PILRB","STAG3L5P","PILRA","ZCWPW1", "STAG3", "MEPCE", "PVRIG")
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
eqtl.to.phenotype <- function()
{
   source("~/github/gwasExplorer/explore/newScore/doesLinkagePredictBraakscSeparatedEnrichment/eqtl-to-phenotype.R")

   targetGene <- "PILRB"
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")

   RSID <- tag.snp
   x <- rosmap.ad.ctl.separation(RSID)
   checkEquals((names(x)), c("pval.t", "pval.fisher", "tbl.geno", "tbl.pt", "mtx.geno", "mtx.geno.study", "pt.ad", "pt.ctl"))

   tbl.rosmap.nearby.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   roi <- getGenomicRegion(igv)
   goi <- "PILRB"
   pval.threshold <- 1e-8
   tbl.eqtl <- subset(tbl.rosmap.nearby.raw, hg38 >= roi$start & hg38 <= roi$end &
                                             pvalue < pval.threshold)
                                             # & genesymbol == goi)
   dim(tbl.eqtl)
   tbl.track <- tbl.eqtl[, c("chrom", "hg38", "hg38", "pvalue")]
   colnames(tbl.track) <- c("chrom", "start", "end", "score")
   tbl.track$start <- tbl.track$start - 1
   tbl.track$score <- -log10(tbl.track$score)
   track <- DataFrameQuantitativeTrack("rosmap eQTL", tbl.track, autoscale=TRUE, color="blue")
   displayTrack(igv, track)

   rsids <- unique(tbl.eqtl$rsid)
   length(rsids)
   xx <- lapply(rsids, rosmap.ad.ctl.separation)
   length(xx)
   names(xx) <- rsids

   pval.fisher <- unlist(lapply(xx, function(el) el$pval.fisher))
   pval.t <- unlist(lapply(xx, function(el) el$pval.t))
   tbl.rsids.braak <- data.frame(rsid=rsids, pval.fisher=pval.fisher, pval.t=pval.t)
   save(xx, tbl.rsids.braak, file="braakScoreEvaluated-643-rsids.RData")
   tbl.track <- unique(merge(tbl.eqtl[, c("chrom", "hg38", "hg38", "rsid")], tbl.rsids.braak, by="rsid"))
   tbl.track$score <- -log10(tbl.track$pval.fisher)
   colnames(tbl.track) <- c("rsid", "chrom", "start", "end", "pval.fisher", "pval.t", "score")
   tbl.track <- tbl.track[, c("chrom", "start", "end", "score", "pval.fisher", "pval.t", "rsid")]
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameQuantitativeTrack("braak sep", tbl.track, color="magenta", autoscale=FALSE,
                                       min=0, max = 1.05 * max(tbl.track$score))
   displayTrack(igv, track)



} # eqtl.to.phenotype
#----------------------------------------------------------------------------------------------------
view.gtex.eqtls <- function()
{

} # veiw.gtex.eqtls
#----------------------------------------------------------------------------------------------------
if(!interactive())
   run.all()
