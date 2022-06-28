library(GenomicRanges)
library(RUnit)


targetGene <- "NCK2"

gtex.brain.tissues <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
                        "GTEx_V8.Brain_Cerebellum",
                        "GTEx_V8.Brain_Cortex",
                        "GTEx_V8.Brain_Frontal_Cortex_BA9",
                        "GTEx_V8.Brain_Hippocampus",
                        "GTEx_V8.Brain_Hypothalamus")


if(!exists("tbl.bellenguez")){
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
    file.exists(data.dir)
    filename <- "bellenguez-2022-83variants-hg38.RData"
    full.path <- file.path(data.dir, filename)
    tbl.bellenguez <- get(load(full.path))
    subset(tbl.bellenguez, rsid==rsid.oi)
    filename <- "posthuma-2019-with-hg38-locs.RData"
    full.path <- file.path(data.dir, filename)
    tbl.posthuma <- get(load(full.path))
    subset(tbl.posthuma, rsid==rsid.oi)
    } # tbl.bellenguez


if(!exists("tbl.gwas3")){ # this 235 row, 215 rsid table, combines posthuma, schartzentruber, bellenguez
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
    file.exists(data.dir)
    filename <- "tbl.gwas.3studies.235x5.RData"
    full.path <- file.path(data.dir, filename)
    tbl.gwas3 <- get(load(full.path))
    stopifnot(colnames(tbl.gwas3)[c(1,2,5)] == c("chrom", "hg38", "pvalue"))
    head(tbl.gwas3)
    } # tbl.gwas3

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
   } # tbl.gwascat.ad



if(!exists("igv")){
   igv <- start.igv(targetGene, "hg38")
   zoomOut(igv)
   zoomOut(igv)
   }
tag.snp <- "rs143080277"
tag.snp.chrom <- "chr2"
tag.snp.hg38 <- 105749599
roi.1m <- "chr2:105249599-106249599"

if(!exists("tbl.hap")){
    tbl.hap <- read.table("../shared/haploreg.tsv", header=TRUE, sep="\t", as.is=TRUE)
    tbl.hap$score <- (1.00001 - tbl.hap$Rsquared)/100
    track <- GWASTrack("tag.snp and LD", tbl.hap, chrom.col=1, pos.col=2, pval.col=8)
    displayTrack(igv, track)
    }

if(!exists("tbl.eqtl.rosmap")){
    require(RPostgreSQL)
    geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
    tag.snp.loc <- tag.snp.hg38
    shoulder <- 500000
    start <- tag.snp.loc - shoulder
    end   <- tag.snp.loc + shoulder
    dbGetQuery(geneRegDB, "select * from eqtl2 limit 3")
    query <- sprintf("select * from eqtl2 where chrom='%s' and hg38 >= %d and hg38 <= %d",
                     tag.snp.chrom, start, end)
    tbl <- dbGetQuery(geneRegDB, query)
    dim(tbl) # 30725    12
        # rename C2or40 to ECRG4
    c2orf40.rows <- grep("C2orf40", tbl$genesymbol)
    length(c2orf40.rows)
    tbl$genesymbol[c2orf40.rows] <- "ECRG4"
        # fill in empty genesymbol for  ENSG00000272994: C2orf49-DT
    ENSG00000272994.rows <- grep("ENSG00000272994", tbl$ensg)
    length(ENSG00000272994.rows)  #88
    tbl$genesymbol[ENSG00000272994.rows] <- "C2orf49-DT"

    tbl.eqtl.rosmap <- subset(tbl, pvalue < 1e-3)
    dim(tbl.eqtl.rosmap)  # 332 at 1e-3, 169 at 1e-5
    unique(tbl.eqtl.rosmap$genesymbol)
    as.data.frame(sort(table(tbl.eqtl.rosmap$genesymbol), decreasing=TRUE))
    tbl.eqtl.rosmap$score <- -log10(tbl.eqtl.rosmap$pvalue) * tbl.eqtl.rosmap$beta
    save(tbl.eqtl.rosmap, file="tbl.eqtls.rosmap.RData")  # -17.996464  -5.583823  -1.865126   2.332745   3.857272
    if(FALSE){
       geneTargets <- unique(tbl.eqtl.rosmap$genesymbol)
       min.score <- min(tbl.eqtl.rosmap$score)
       max.score <- max(tbl.eqtl.rosmap$score)
       for(gene in geneTargets){
          tbl.sub <- subset(tbl.eqtl.rosmap, genesymbol==gene)
          tbl.track <- tbl.sub[, c("chrom", "hg38", "hg38", "score")]
          colnames(tbl.track) <- c("chrom", "start", "end", "score")
          tbl.track$start <- tbl.track$start -1
          track <- DataFrameQuantitativeTrack(gene, tbl.track, autoscale=FALSE,
                                              min=min.score, max=max.score, color="darkblue")
          displayTrack(igv, track)
          } # for gene
    } # if FALSE


#lapply(models, function(model) head(model$trena, n=5))
#lapply(models, function(model) head(model$tbl.breaks, n=5))

# prior formula
# motifBreak.score <- with(tbl.tf, abs(pctDelta) * 100)
#      eqtl.score <- with(tbl.tf, -log10(gtex.eqtl.pval)* abs(gtex.eqtl.beta) * 100)
#     trena.score <- with(tbl.tf, (abs(betaLasso) * 100) + (rfNorm * 10))
#      tfbs.score <- 1/tfbs.count
#score <- trena.score * eqtl.score * motifBreak.score * tfbs.score

#names(models)


library(TrenaProjectAD)
library(EndophenotypeExplorer)
library(RUnit)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
source("~/github/tms-makeAndBreak/R/tmsMB-class.R")
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

targetGene <- "PILRB"

tag.snp <- "rs7384878"
tag.snp.chrom <- "chr7"
tag.snp.loc <- 100334426
rsid.oi <- tag.snp

pilrb.snp.1 <- list(rsid="rs3087502",   hg38=100357741, chrom="chr7",  tf="TEAD1")
pilrb.snp.2 <- list(rsid="rs11765869",  hg38=100367166, chrom="chr7",  tf="ZEB1")
pilrb.snp.3 <- list(rsid="rs6946768",   hg38=100356770, chrom="chr7",  tf="HSF1")
pilrb.snp.4 <- list(rsid="rs113261079", hg38=100340180, chrom="chr7", tf="MEF2C")

if(FALSE){
    tbl.track <- with(pilrb.snp.1, data.frame(chrom=chrom,
                                              start=hg38-1,
                                              end=hg38,
                                              name=rsid,
                                              stringsAsFactors=FALSE))
    track <- DataFrameAnnotationTrack(pilrb.snp.1$rsid, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    tbl.track <- with(pilrb.snp.2, data.frame(chrom=chrom,
                                              start=hg38-1,
                                              end=hg38,
                                              name=rsid,
                                              stringsAsFactors=FALSE))
    track <- DataFrameAnnotationTrack(pilrb.snp.2$rsid, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    tbl.track <- with(pilrb.snp.3, data.frame(chrom=chrom,
                                              start=hg38-1,
                                              end=hg38,
                                              name=rsid,
                                              stringsAsFactors=FALSE))
    track <- DataFrameAnnotationTrack(pilrb.snp.3$rsid, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    tbl.track <- with(pilrb.snp.4, data.frame(chrom=chrom,
                                              start=hg38-1,
                                              end=hg38,
                                              name=rsid,
                                              stringsAsFactors=FALSE))
    track <- DataFrameAnnotationTrack(pilrb.snp.4$rsid, tbl.track, color="red", trackHeight=25)
    displayTrack(igv, track)
    }

if(FALSE){  # create motifbreakR logos
    pilrb.snp <- pilrb.snp.3
    motifs <- query(MotifDb, c("jaspar2022", "sapiens", pilrb.snp$tf))
    snps.gr <- snps.from.rsid(rsid = pilrb.snp$rsid,
                              dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                              search.genome=BSgenome.Hsapiens.UCSC.hg38)
    results <- motifbreakR(snpList = snps.gr,
                           filterp = TRUE,
                           pwmList = motifs,
                           show.neutral=FALSE,
                           method = c("ic", "log", "notrans")[1],
                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose=TRUE)
    plotMB(results[1], rsid=pilrb.snp$rsid, effect=c("strong")) # takes a minute or two
    } # plotMB

gtex.brain.tissues <- c("GTEx_V8.Brain_Cerebellar_Hemisphere",
                        "GTEx_V8.Brain_Cerebellum",
                        "GTEx_V8.Brain_Cortex",
                        "GTEx_V8.Brain_Frontal_Cortex_BA9",
                        "GTEx_V8.Brain_Hippocampus",
                        "GTEx_V8.Brain_Hypothalamus")


if(!exists("igv") & FALSE){
    igv <- start.igv(targetGene)
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
    tbl.haploreg.strong <- subset(tbl.haploreg, Rsquared >= 0.5)
    dim(tbl.haploreg)  # 109 7
    dim(tbl.haploreg.strong) # 10
    tbl.track <- tbl.haploreg.strong[, c("chrom", "hg38", "hg38", "Rsquared")]
    colnames(tbl.track) <- c("chrom", "start", "end", "score")
    tbl.track$chrom <- paste0("chr", tbl.track$chrom)
    tbl.track$start <- tbl.track$start - 1
    if(FALSE){
       track <- DataFrameQuantitativeTrack("LD", tbl.track, color="black", autoscale=TRUE)
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


if(!exists("tbl.gwas3")){ # this 235 row, 215 rsid table, combines posthuma, schwartzentruber, bellenquez
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
test_combine.tables.pilrb <- function()
{
    load("PILRB-models-roi.528k-2022.jun.24-13:02.RData")
    targetGene <- "PILRB"
    names(models)

    # x <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF="TEAD1")

    for(name in names(models)){
        model <- models[[name]]
        tbl.trena <- head(model$trena, n=10)
        tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
        tbl.breaks <- model$tbl.breaks
        tbl.tms <- model$tms.filtered
        tbl.tms$targetGene <- targetGene
        scores <- unlist(lapply(tbl.trena$gene, function(gene) {
            x <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF=gene)
            x$trenaScore
        }))
        printf("%50s: %5d", name, sum(scores))
        } # for name

    checkTrue(sum(scores) > 8000 & sum(scores) < 10000)
    tbl.trena$score <- scores

    checkEquals(ncol(tbl.out), 30)
    browser(); abc <- 99

} # test_combine.tables
#----------------------------------------------------------------------------------------------------
test_combine.tables.stag3 <- function()
{
    message(sprintf("--- test_combine.tables.stag3"))

    load("STAG3-models-roi.528k-2022.jun.22-13:49.RData")
    names(models)
    lapply(models, function(model) dim(model$trena))
    model <- models[[1]]
    targetGene <- "STAG3"
    tbl.trena <- head(model$trena, n=10)
    tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
    tbl.breaks <- model$tbl.breaks
    tbl.tms <- model$tms.filtered
    tbl.tms$targetGene <- targetGene

    x <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF="TFAP4")

    scores <- unlist(lapply(tbl.trena$gene, function(gene) {
        x <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene=targetGene, TF=gene)
        x$trenaScore
        }))

    tbl.trena$score <- scores

    checkEquals(ncol(tbl.out), 30)
    browser(); abc <- 99


} # test_combine.tables.stag3
#------------------------------------------------------------
viz <- function()
{
    load("../pilrb/PILRB-models-roi.528k-2022.jun.24-13:02.RData")
    names(models)
    model <- models[[3]]
    tbl.trena <- head(model$trena, n=10)
    tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
    tbl.breaks <- model$tbl.breaks
    tbl.tms <- model$tms.filtered
    tbl.tms$targetGene <- targetGene
    track.names <- getTrackNames(igv)
    removeTracksByName(igv, track.names[7:length(track.names)])
    for(TF in tbl.trena$gene){
        tbl.track <- subset(tbl.tms, tf==TF)[, c("chrom", "start", "end")]
        track <- DataFrameAnnotationTrack(TF, tbl.track, color="random", trackHeight=25)
        displayTrack(igv, track)
        }


} # viz
#------------------------------------------------------------------------------------------------------------------------
