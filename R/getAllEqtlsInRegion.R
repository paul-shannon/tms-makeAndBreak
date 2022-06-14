library(EndophenotypeExplorer)
library(ADvariantExplorer)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(ebi.eqtls)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_build.snpLocs.fromScratch()
    test_get.ampad.eqtls()
    test_get.ebi.eqtls()

} # runTests
#----------------------------------------------------------------------------------------------------
build.snpLocs.fromScratch <- function(chromosome, center.hg38, shoulder)
{
   message("")
   message(sprintf("--- about to build tbl.snpLocs from scratch, width=%d", shoulder))
   message("")

       #--------------
       # hg38 first
       #--------------

   gr <- GRanges(seqnames=chromosome, #sub("chr", "", chromosome),
                 IRanges(start=center.hg38-shoulder, end=center.hg38+shoulder))
   seqlevelsStyle(gr) <- "NCBI"
   gr.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
   seqlevelsStyle(gr.snps) <- "UCSC"   # needed for rtracklayer liftover

   chain.file <- "hg38ToHg19.over.chain.gz"
   if(!file.exists(chain.file)){
      system(sprintf("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/%s",
                     chain.file))
      system(sprintf("gunzip %s", chain.file))
      }
   chain <- import.chain(sub(".gz", "", chain.file, fixed=TRUE))
   x <- liftOver(gr.snps, chain)
   gr.hg19 <- unlist(x)
   tbl.snpLocs.hg19 <- as.data.frame(gr.hg19)[, c("seqnames", "start", "RefSNP_id")]
   colnames(tbl.snpLocs.hg19) <- c("chrom", "hg19", "rsid")

   tbl.snpLocs.hg38 <- as.data.frame(gr.snps)[, c(1,2,4)]
   colnames(tbl.snpLocs.hg38) <- c("chrom", "hg38", "rsid")
   tbl.snpLocs.hg38$chrom <- as.character(tbl.snpLocs.hg38$chrom)
   dim(tbl.snpLocs.hg38)

   tbl.snpLocs <- merge(tbl.snpLocs.hg38, tbl.snpLocs.hg19[, c("rsid", "hg19")], by=c("rsid"), all=TRUE)
   rownames(tbl.snpLocs) <- tbl.snpLocs$rsid
   coi <- c("chrom", "hg19", "hg38", "rsid")
   tbl.snpLocs <- tbl.snpLocs[, coi]
   rownames(tbl.snpLocs) <- NULL

   tbl.snpLocs

} # function
#----------------------------------------------------------------------------------------------------
test_build.snpLocs.fromScratch <- function()
{
    message(sprintf("--- test_build.snpLocs.fromScratch"))
      # rs7904826
      # 10:58385197 (GRCh38)
      # 10:60144957 (GRCh37)
    tbl.snpLocs <-
        build.snpLocs.fromScratch(chromosome="chr10",
                                  center.hg38=58385197,
                                  shoulder=3)

    checkEquals(dim(tbl.snpLocs), c(3, 4))
    checkEquals(as.list(subset(tbl.snpLocs, rsid=="rs7904826")),
                list(chrom="chr10", hg19=60144957, hg38=58385197, rsid="rs7904826"))

       # rs1000000404
       # 10:57663008 (GRCh38)
       # 10:59422768 (GRCh37)

    tbl.snpLocs <-
        build.snpLocs.fromScratch(chromosome="chr10",
                                  center.hg38=57663008,
                                  shoulder=1)
    checkEquals(as.list(subset(tbl.snpLocs, rsid=="rs1000000404")),
                list(chrom="chr10", hg19=59422768, hg38=57663008, rsid="rs1000000404"))


      #--------------------------------
      # now a larger region: +/- 1000
      #--------------------------------

    tbl.snpLocs <-
        build.snpLocs.fromScratch(chromosome="chr10",
                                  center.hg38=57663008,
                                  shoulder=1000)
        # previous test should still pass
    checkEquals(as.list(subset(tbl.snpLocs, rsid=="rs1000000404")),
                list(chrom="chr10", hg19=59422768, hg38=57663008, rsid="rs1000000404"))

      #-----------------------------
      # and larger still: +/- 10k
      #-----------------------------

    tbl.snpLocs <-
        build.snpLocs.fromScratch(chromosome="chr10",
                                  center.hg38=57663008,
                                  shoulder=10000)
    checkEquals(as.list(subset(tbl.snpLocs, rsid=="rs1000000404")),
                list(chrom="chr10", hg19=59422768, hg38=57663008, rsid="rs1000000404"))


      #-----------------------------
      # and larger still: +/- 1M
      #-----------------------------

    tbl.snpLocs <-
        build.snpLocs.fromScratch(chromosome="chr10",
                                  center.hg38=57663008,
                                  shoulder=1e6)
    checkEquals(as.list(subset(tbl.snpLocs, rsid=="rs1000000404")),
                list(chrom="chr10", hg19=59422768, hg38=57663008, rsid="rs1000000404"))
    checkTrue(nrow(tbl.snpLocs) > 650000)
    hg19.missing <- length(which(is.na(tbl.snpLocs$hg19)))
    hg38.missing <- length(which(is.na(tbl.snpLocs$hg38)))
    checkEquals(hg19.missing, 0)
    checkEquals(hg38.missing, 0)

} # test_build.snpLocs.fromScratch
#----------------------------------------------------------------------------------------------------
test_get.ampad.eqtls <- function()
{
   message(sprintf("--- test_get.ampad.eqtls"))

   tbl.eqtls <- get.ampad.eqtls()
   checkTrue(nrow(tbl.eqtls) > 13000)
   expected <- c("chrom", "hg19", "hg38", "rsid", "pvalue", "ensg", "genesymbol", "study", "tissue", "assay")
   checkEquals(colnames(tbl.eqtls), expected)
   tbl.freq <- as.data.frame(table(tbl.eqtls$study))
   checkTrue(subset(tbl.freq, Var1=="ampad-mayo")$Freq > 8000)
   checkTrue(subset(tbl.freq, Var1=="ampad-rosmap")$Freq > 4000)

} # test_get.ampad.eqtls
#----------------------------------------------------------------------------------------------------
test_get.ebi.eqtls <- function()
{
   message(sprintf("--- test_get.ebi.eqtls"))

   chrom <- "chr10"
   center.hg38 <- 58385409
   tbl.small.no.chunks <- fetch.eqtls(chrom,
                                      start=center.hg38 - 100000,
                                      end=center.hg38 + 100000,
                                      study="GTEx_V8.Brain_Cerebellum",
                                      simplify=TRUE, chunk.size=10000)

   viz <- FALSE
   if(viz){
       tbl.2 <- merge(tbl.small.no.chunks, tbl.snpLocs, by="rsid")
       tbl.2 <- subset(tbl.2, gene=="TFAM")
       sprintf("chr10:%d-%d", min(tbl.2$hg38), max(tbl.2$hg38))
       tbl.2$score <- with(tbl.2, -log10(pvalue) * beta)
       tbl.track <- tbl.2[, c("chrom", "hg38", "hg38", "score")]
       colnames(tbl.track)[2:3] <- c("start", "end")
       tbl.track$start <- tbl.track$start-1
       track <- DataFrameQuantitativeTrack("new2 snp", tbl.track, color="darkgreen", autoscale=TRUE)
       displayTrack(igv, track)
       }

   checkEquals(colnames(tbl.small.no.chunks), c("rsid", "pvalue", "gene", "total.alleles", "beta", "id"))
   checkTrue(nrow(tbl.small.no.chunks) > 10)

   tbl.small.2k.chunks <- fetch.eqtls(chrom,
                                      start=center.hg38 - 10000,
                                      end=center.hg38 + 10000,
                                      study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                      simplify=TRUE, chunk.size=2000)
   dim(tbl.small.2k.chunks)
   checkTrue(nrow(tbl.small.2k.chunks) > 600)

   tbl.small.2k.chunks.2studies <- fetch.eqtls(chrom,
                                               start=center.hg38 - 10000,
                                               end=center.hg38 + 10000,
                                               study=c("GTEx_V8.Brain_Cerebellar_Hemisphere",
                                                       "GTEx_V8.Brain_Hippocampus"),
                                               simplify=TRUE, chunk.size=2000)
   dim(tbl.small.2k.chunks.2studies)
   checkTrue(nrow(tbl.small.2k.chunks.2studies) > 1000)
   checkEquals(sort(unique(tbl.small.2k.chunks.2studies$id)),
               c("GTEx_V8.Brain_Cerebellar_Hemisphere", "GTEx_V8.Brain_Hippocampus"))

   run.big <- FALSE
   if(run.big){
       tbl.big <- fetch.eqtls(chrom,
                              start=center.hg38 - 100000,
                              end=center.hg38 + 100000,
                              study=c("GTEx_V8.Brain_Cerebellar_Hemisphere", "GTEx_V8.Brain_Hippocampus"),
                              simplify=TRUE, chunk.size=5000)
      dim(tbl.big)
      checkTrue(nrow(tbl.big) > 10000)
      checkEquals(sort(unique(tbl.big$id)),
                   c("GTEx_V8.Brain_Cerebellar_Hemisphere", "GTEx_V8.Brain_Hippocampus"))
      } # if run.big

} # test_get.ebi.eqtls
#----------------------------------------------------------------------------------------------------
get.ampad.eqtls <- function()
{
    tbl.eqtls <- etx$get.ampad.EQTLsForGene()
    viz <- FALSE
    if(viz){
       dim(tbl.eqtls)
       head(tbl.eqtls, n=20)
       tail(tbl.eqtls, n=20)
       tbl.track <- subset(tbl.eqtls, pvalue < 1e-2 & study != "GTEx")[, c("chrom", "hg19", "hg19", "rsid", "pvalue")]
       dim(tbl.track)
       colnames(tbl.track)[2:3] <- c("start", "end")
       tbl.track$start <- tbl.track$start - 1
       colnames(tbl.track) <- NULL
       track <- GWASTrack("ampad", tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
       displayTrack(igv, track)
       }
   invisible(tbl.eqtls)

} # ampad.eqtls
#----------------------------------------------------------------------------------------------------
fetch.eqtls <- function(chrom, start, end, study, simplify, chunk.size)
{
  roi.width <- 1 + end - start

  if(roi.width <= chunk.size){
     message(sprintf("--- fetch.eqtls, just one chunk"))
     tbl <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, study, simplify=simplify)

  } else {
     boundaries.needed <- 1 + (roi.width %/% chunk.size)
     starts <- as.integer(seq(from=start, to=end, length.out=boundaries.needed))
     ends <- starts[-1]
     starts <- starts[-(length(starts))]
     tbls <- list()
     intervals <- length(starts)
     message(sprintf("==== fetch.eqtls, %d chunks", intervals))
     for(i in seq_len(intervals)){
        message(sprintf("--- fetching chunk %2d/%d for %s", i, intervals, study))
        tbl.chunk <- avx$geteQTLsByLocationAndStudyID(chrom,
                                                      as.integer(starts[i]),
                                                      as.integer(ends[i]),
                                                      study,
                                                      targetGene.only=FALSE,
                                                      simplify=simplify)
        tbls[[i]] <- tbl.chunk
        } # for i
     tbl <- do.call(rbind, tbls)
     } # else

  invisible(tbl)

} # fetch.eqtls
#----------------------------------------------------------------------------------------------------
fetch.all.gtex.brain.eqtls <- function(chrom, center.hg38, shoulder, chunk.size)
{
   loc <- list(chrom=chrom,
               start=center.hg38 - shoulder,
               end=center.hg38 + shoulder,
               string = sprintf("%s:%d-%d", chrom, center.hg38 - shoulder, center.hg38 + shoulder))

   message(sprintf("--- fetching gtex brain eqtls, region size: %dk", as.integer(loc$width/1000)))

   selected.studies <- c("GTEx_V8.Brain_Amygdala",
                         "GTEx_V8.Brain_Anterior_cingulate_cortex_BA24",
                         "GTEx_V8.Brain_Caudate_basal_ganglia",
                         "GTEx_V8.Brain_Cerebellar_Hemisphere",
                         "GTEx_V8.Brain_Cerebellum",
                         "GTEx_V8.Brain_Cortex",
                         "GTEx_V8.Brain_Frontal_Cortex_BA9",
                         "GTEx_V8.Brain_Hippocampus",
                         "GTEx_V8.Brain_Hypothalamus",
                         "GTEx_V8.Brain_Nucleus_accumbens_basal_ganglia",
                         "GTEx_V8.Brain_Putamen_basal_ganglia",
                         "GTEx_V8.Brain_Spinal_cord_cervical_c-1",
                         "GTEx_V8.Brain_Substantia_nigra")[7] #c(5,6,7,8,9)]

   filename <- sprintf("tbl.eqtls.gtex.%d.brain.tissues.%s.RData", length(selected.studies), loc$string)

   for(study in selected.studies){
      tbl.eqtl <- with(loc, ebi.eqtls$fetch.eqtls.in.chunks(chrom, start, end,
                                                     study=study,
                                                     simplify=TRUE,
                                                     chunk.size=chunk.size))
      study.title <- sub("GTEx_V8.Brain_", "", study)
      study.filename <- sprintf("tbl.eqtls.gtex.%s.%s.RData", study.title, loc$string)
      dim(tbl.eqtl)
      deleters <- which(is.na(tbl.eqtl$rsid))
      length(deleters)
      if(length(deleters) > 0)
          tbl.eqtl <- tbl.eqtl[-deleters,]
      tbl.eqtl <- subset(tbl.eqtl, pvalue < 0.05)
      tbl.eqtl <- merge(tbl.eqtl, tbl.snpLocs, by="rsid", all.x=TRUE)
      tbl.eqtl$score <- -log10(tbl.eqtl$pvalue) * abs(tbl.eqtl$beta)
      dim(tbl.eqtl)
      new.order <- order(tbl.eqtl$score, decreasing=TRUE)
      tbl.eqtl <- tbl.eqtl[new.order,]
      tbl.eqtl$name <- sprintf("%s-%s", tbl.eqtl$rsid, tbl.eqtl$gene)
      rownames(tbl.eqtl) <- NULL
      deleters <- grep("RP11", tbl.eqtl$gene)
      if(length(deleters) > 0)
          tbl.eqtl <- tbl.eqtl[-deleters,]
      # tbls.eqtl[[study.title]] <- tbl.eqtl
      message("")
      message(sprintf("---- saving %d %s eqtls to %s", nrow(tbl.eqtl), study.title, filename))
      message("")
      save(tbl.eqtl, file=study.filename)
      examine <- FALSE
      if(examine){
         tbl.eqtl.strong <- subset(tbl.eqtl, score > 0.2)
         gois <- sort(unique(tbl.eqtl.strong$gene))
         for(goi in gois){
             tbl.track <- subset(tbl.eqtl.strong, gene==goi)[, c("chrom", "hg38", "hg38", "score", "rsid")]
             deleters <- which(is.na(tbl.track$hg38))
             if(length(deleters) > 0)
             tbl.track <- tbl.track[-deleters,]
             colnames(tbl.track) <- c("chrom", "start", "end", "score", "name")
             tbl.track$chrom <- sprintf("chr%s", tbl.track$chrom)
             tbl.track$start <- tbl.track$start - 1
             dim(tbl.track)
             track.title <- sprintf("%s %s", study.title, goi)
             track <- DataFrameQuantitativeTrack(track.title, tbl.track, autoscale=FALSE,
                                                 min=0, max=5, color="random",
                                                 trackHeight=25)
             displayTrack(igv, track)
             } # for goi
          } # if examine
      } # for each brain study

   tbls.eqtl <- list()
   fs <- list.files(".", "tbl.*RData")

   for(f in fs){
      tissue <- sub("tbl.eqtls.gtex\\.", "", f)
      tissue <- sub("\\.chr.*RData", "", tissue)
      tbl.eqtl <- get(load(f))
      tbl.eqtl$tissue <- tissue
      tbls.eqtl[[tissue]] <- tbl.eqtl
      } # for f


   tbl.all <- do.call(rbind, tbls.eqtl)
   save(tbl.all, file=filename)
   invisible(tbl.all)

} # fetch.all.gtex.brain.eqtls
#----------------------------------------------------------------------------------------------------
ebi.eqtls <- ebi.eqtls$new()

if(!interactive()){
    args <- commandArgs(trailingOnly=TRUE)
    stopifnot(length(args) == 5)
    targetGene <- args[1]
    chrom <- args[2]
    center.hg38 <- as.numeric(args[3])
    shoulder <-  as.numeric(args[4])
    chunk.size <- as.numeric(args[5])
    snpLocs.filename <- sprintf("snpLocs-%s.RData", targetGene)
    if(!exists("tbl.snpLocs")){
        if(!file.exists(snpLocs.filename)){
            msg <- "need to create tbl.snpLocs for region, save as .RData"
            tbl.snpLocs <- build.snpLocs.fromScratch(chrom, center.hg38, shoulder)
            message(sprintf("saving tbl.snpLocs with %d rows", nrow(tbl.snpLocs)))
            save(tbl.snpLocs, file=snpLocs.filename)
            }
        load(snpLocs.filename)
        } # exists tbl.snpLocs
    if(!exists("etx"))
        etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
    if(!exists("avx"))
        avx <- ADvariantExplorer$new(targetGene, chrom, center.hg38-shoulder, center.hg38+shoulder)
    tbl.ampad.eqtls <- get.ampad.eqtls()
    save(tbl.ampad.eqtls, file=sprintf("ampad.eqtls.%s.RData", targetGene))
    tbl.eqtl <- fetch.all.gtex.brain.eqtls(chrom, center.hg38, shoulder=shoulder, chunk.size=chunk.size)
} else {
    targetGene <- "TFAM"
    chrom <- "chr10"
    center.hg38 <- 58385407
    shoulder <- 1000
    chunk.size <- 2000
    snpLocs.filename <- sprintf("snpLocs-%s.RData", targetGene)
    if(!exists("tbl.snpLocs")){
        if(!file.exists(snpLocs.filename)){
            msg <- "need to create tbl.snpLocs for region, save as .RData"
            tbl.snpLocs <- build.snpLocs.fromScratch(chrom, center.hg38, shoulder)
            message(sprintf("saving tbl.snpLocs with %d rows", nrow(tbl.snpLocs)))
            save(tbl.snpLocs, file=snpLocs.filename)
            }
        load(snpLocs.filename)
        } # exists tbl.snpLocs
    if(!exists("etx"))
        etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
    if(!exists("avx"))
        avx <- ADvariantExplorer$new(targetGene, chrom, center.hg38-shoulder, center.hg38+shoulder)
    tbl.eqtl <- fetch.all.gtex.brain.eqtls(chrom, center.hg38, shoulder, chunk.size)
    }
