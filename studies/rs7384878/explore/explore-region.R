library(TrenaProjectAD)
library(EndophenotypeExplorer)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)

tag.snp <- "rs7384878"
tag.snp.chrom <- "chr7"
tag.snp.loc <- 100334426
rsid.oi <- tag.snp


igv <- start.igv("all", "hg38")
if(!exists("igv")){
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
    track <- GWASTrack("LD", tbl.haploreg, chrom.col=1, pos.col=2, pval.col=8, trackHeight=100)
    displayTrack(igv, track)
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
   data.dir <- system.file(package="igvR", "extdata", "gwas")
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

if(!exists("tbl.rsids")){  # with mayo & rosmap scores for association of each variant to extreme braak score
   data.dir <- "~/github/gwasExplorer/explore/newScore/doesLinkagePredictBraakscSeparatedEnrichment"
   filename <- "tbl.rsids.summary.Tue.May.31.2022-13:26:55-hetAndHomSeparate.RData"
   full.path <- file.path(data.dir, filename)
   stopifnot(file.exists(full.path))
   tbl.rsids <- get(load(full.path))
   }

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
build.model <- function(gene, gtex.tissue, roi)
{
   tbl.eqtls <- subset(tbl.rosmap.nearby, genesymbol==gene & hg38 > roi$start & hg38 < roi$end)
   gr.eqtl <- with(tbl.eqtls, GRanges(seqnames=chrom[1], IRanges(start=hg38)))
   tbl.ov <- as.data.frame(findOverlaps(gr.eqtl, gr.fimo))
   dim(tbl.ov)
   tbl.fimo.sub <- tbl.fimo[tbl.ov[,2],]
   tms <- TMS$new(TrenaProjectAD(), gene, tbl.fimo.sub, quiet=FALSE)
   tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
   etx <- EndophenotypeExplorer$new(gene, "hg38", vcf.project="AMPAD")
   names(etx$get.rna.matrix.codes())
   mtx.rna <- etx$get.rna.matrix(gtex.tissue))
   tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

   tbl.tms <- tms$getTfTable()
   candidate.tfs <- subset(tbl.tms, abs(cor.all) > 0.3)$tf
   length(candidate.tfs)
   tbl.trena <- tms$build.trena.model(candidate.tfs, mtx.rna)
   tbl.trena$rfNorm <- with(tbl.trena, rfScore/max(rfScore))
   tfs.oi <- subset(tbl.trena, abs(betaLasso) > 0.1 | rfNorm > 0.5)$gene
   rsids <- tbl.eqtls$rsid
   tbls.breaks <- list()
   for(tf in tfs.oi){
      printf("--- starting on tf %s for targetGene %s, %d rsids", tf, gene, length(rsids))
      motifs <- query(MotifDb, c(tf, "sapiens", "jaspar2022"))
      snps.gr <- snps.from.rsid(rsids,
                                dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                search.genome=BSgenome.Hsapiens.UCSC.hg38)
      results <- motifbreakR(snpList = snps.gr,
                             filterp = TRUE,
                             pwmList = motifs,
                             show.neutral=FALSE,
                             method = c("ic", "log", "notrans")[1],
                             bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                             BPPARAM = SerialParam(),
                             verbose=TRUE)
      tbl.breaks <- as.data.frame(results, row.names=NULL)
      tbl.breaks <- subset(tbl.breaks, effect=="strong")
      tbl.breaks$pctDelta <- tbl.breaks$pctAlt - tbl.breaks$pctRef
      tbls.breaks[[tf]] <- tbl.breaks
      print(tbl.breaks[, c("SNP_id", "geneSymbol", "pctRef", "pctAlt", "pctDelta")])
      browser(); xyz <- 99
      } # for tf

   return(tbls.breaks)


} # build.model
#----------------------------------------------------------------------------------------------------
run.all <- function()
{
    # set igv's region to include only tag.snp's LD region, also considering the gwascat.ad
    # variants, and the rosmapeQTLS
    # for rs7384878 (pilra-region) it is chr7:100,088,993-100,600,689
   roi <- getGenomicRegion(igv)
   tbl.affected <- identify.affected.genes()
   median <- median(tbl$sum.sig.x.beta)
   tbl.affected.top <- subset(tbl.affected, sum.sig.x.beta >= median)
   goi <- as.character(tbl.affected.top$gene)
   for(gene in goi){
      build.model(gene, roi)

#----------------------------------------------------------------------------------------------------
test_build.model <- function()
{
   build.model("ZSCAN21")
   x <- build.model("PILRB")

} # test_build.model
#----------------------------------------------------------------------------------------------------
view.gtex.eqtls <- function()
{

} # veiw.gtex.eqtls
#----------------------------------------------------------------------------------------------------

