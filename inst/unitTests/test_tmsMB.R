library(RUnit)
library(TrenaProjectAD)
library(tmsMB)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_setStudyRegion()
   test_tissue()
   test_run.tms()
   test_addEqtlsToTMStable()
   test_run.trena()
   test_run.small()
   test_run.2k.break.motifs()
   test_run.large()

} # runTests
#----------------------------------------------------------------------------------------------------
targetGene <- "PILRA"
tag.snp <- "rs7384878"
tag.snp.chrom <- "chr7"
tag.snp.loc <- 100334426

#----------------------------------------------------------------------------------------------------
trenaProject <- TrenaProjectAD()
data.dir <- "~/github/tms-makeAndBreak/studies/rs7384878/shared"
if(!exists("tbl.fimo")){
   tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.PILRA.RData")))
   #tbl.ampad.eqtls.raw <- get(load(file.path(data.dir, "ampad.eqtls.pilra.plusMinus.1M.RData")))
   #tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, study=="ampad-rosmap" & pvalue < 0.001)
   tbl.ampad.eqtls.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, genesymbol==targetGene & pvalue < 0.001)
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.001)
   #data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   #filename <- "mayoAllPeaks.merged.96064x4.RData"
   #tbl.mayoAtac <- get(load(file.path(data.dir, filename)))
    #   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    # filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.oc <- data.frame() #tbl.mayoAtac
   tbl.haploreg <- read.table(file.path(data.dir, "haploreg.tsv"), sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
   }

#----------------------------------------------------------------------------------------------------
overview.viz <- function()
{
if(!exists("igv")){
    igv <- start.igv("PILRA", "hg38")
    showGenomicRegion(igv, "chr7:100,304,414-100,343,410") #   good region for TEAD1 regulating PILRA
    }

   tbl.track <- data.frame(chrom=tag.snp.chrom,
                           start=tag.snp.loc-1,
                           end=tag.snp.loc,
                           name=tag.snp,
                           stringsAsFactors=FALSE)
   track <- DataFrameAnnotationTrack(tag.snp, tbl.track, color="red", trackHeight=25)
   displayTrack(igv, track)

  tbl.track <- tbl.haploreg[, c("chrom", "hg38", "hg38", "Rsquared")]
  colnames(tbl.track) <- c("chrom", "start", "end", "score")
  tbl.track$chrom <- paste0("chr", tbl.track$chrom)
  tbl.track$start <- tbl.track$start - 1
  track <- DataFrameQuantitativeTrack("haploreg", tbl.track, autoscale=TRUE, color="darkgray")
  displayTrack(igv, track)

  ghdb <- GeneHancerDB()
  tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
  #tbl.gh$score <- asinh(tbl.gh$combinedscore)
  tbl.gh$score <- tbl.gh$combinedscore
  track <- DataFrameQuantitativeTrack("GH-all", tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
  displayTrack(igv, track)

  tbl.track <- subset(tbl.ampad.eqtls, genesymbol==targetGene)
  track <- GWASTrack("rosmap eqtls", tbl.track, chrom.col=1, pos.col=3, pval.col=5)
  displayTrack(igv, track)

} # overview.viz
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    checkTrue(all(c("R6", "tmsMB") %in% class(tms)))

    region <- tms$getFimoGenomicRegion()
    checkEquals(region$start, 57391973)
    checkEquals(region$end,   59392030)
    checkEquals(region$width.kb, 2000.06)

    region <- tms$getGTEx.eqtl.genomicRegion()
    checkEquals(region$start, 57386057)
    checkEquals(region$end,   59384352)
    checkEquals(region$width.kb, 1998.3)

    region <- tms$getAMPAD.eqtl.genomicRegion()
    checkEquals(region$start, 57391973)
    checkEquals(region$end,   59392030)
    checkEquals(region$width.kb, 2000.06)

    #study.region <- tms$getStudyRegion()
    #checkEquals(study.region$start, 58385244)
    #checkEquals(study.region$end, 58385416)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_setStudyRegion <- function()
{
    message(sprintf("--- test_setStudyRegion"))

    study.region <- tms$getStudyRegion()
    checkEquals(study.region$chrom, "chr10")
    checkEquals(study.region$start, 57391973)
    checkEquals(study.region$end,   59392030)

    new.start <- 58384216
    new.end   <- 58390229
    tms$setStudyRegion(chrom="chr10", start=new.start, end=new.end)
    new.region <- tms$getStudyRegion()

    checkEquals(new.region$chrom, "chr10")
    checkEquals(new.region$start, new.start)
    checkEquals(new.region$end, new.end)
    checkEquals(new.region$width.kb, 6.01)

} # test_setStudRegion
#----------------------------------------------------------------------------------------------------
test_tissue <- function()
{
    message(sprintf("--- test_tissue"))
    tissues <- tms$getGTEx.eqtl.tissues()
    checkEquals(length(tissues), 6)
    for(tissue in tissues){
       tms$set.current.GTEx.eqtl.tissue(tissue)
       checkEquals(tms$get.current.GTEx.eqtl.tissue(), tissue)
       } # for tissue

    TRUE

} # test_tissue
#----------------------------------------------------------------------------------------------------
# add pvalue, beta, and -log10(pvalue) * beta for gtex and ampad eqtls
test_addEqtlsToTMStable <- function()
{
    message(sprintf("--- test_addEqtlsToTMStable"))

    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])
    new.start <- 58383120
    new.end   <- 58383200
    tms$setStudyRegion(chrom="chr10", start=new.start, end=new.end)

    study.region <- tms$getStudyRegion()
    tms$run.tms()
    tbl.tms <- tms$get.tmsTable()
    dim(tbl.tms)

    tms$add.eqtls.toTmsTable()
    tbl.tms <- tms$get.tmsTable()
    expected.new.columns <- c("ampad.eqtl.pval", "ampad.eqtl.beta", "gtex.eqtl.pval", "gtex.eqtl.beta",
                              "ampad.eqtl.score", "gtex.eqtl.score")
    checkTrue(all(expected.new.columns %in% colnames(tbl.tms)))

    checkTrue(is.numeric(tbl.tms$ampad.eqtl.pval))
    checkTrue(is.numeric(tbl.tms$ampad.eqtl.beta))
    checkTrue(is.numeric(tbl.tms$ampad.eqtl.score))

    checkTrue(is.numeric(tbl.tms$gtex.eqtl.pval))
    checkTrue(is.numeric(tbl.tms$gtex.eqtl.beta))
    checkTrue(is.numeric(tbl.tms$gtex.eqtl.score))


    checkEquals(ncol(tbl.tms), 22)
    checkTrue(nrow(tbl.tms) > 100 & nrow(tbl.tms) < 200)
    new.cols <- setdiff(colnames(tbl.tms), colnames(tbl.fimo))
    checkTrue(all(c("chip","phast7","phast100","gh","oc","tss","cor.all") %in% new.cols))

} # test_addEqtlsToTMStable
#----------------------------------------------------------------------------------------------------
test_run.tms <- function()
{
    message(sprintf("--- test_run.tms"))

    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])

    new.start <- 58385244
    new.end   <- 58385416
    tms$setStudyRegion(chrom="chr10", start=new.start, end=new.end)

    tms$run.tms()
    tbl.tms <- tms$get.tmsTable()
    checkEquals(ncol(tbl.tms), 16)
    checkTrue(nrow(tbl.tms) > 500 & nrow(tbl.tms) < 650)
    new.cols <- setdiff(colnames(tbl.tms), colnames(tbl.fimo))
    checkTrue(all(c("chip","phast7","phast100","gh","oc","tss","cor.all") %in% new.cols))

    tms$add.eqtls.toTmsTable()
    tbl.tms <- tms$get.tmsTable()
    checkEquals(ncol(tbl.tms), 22)

} # test_run.tms
#----------------------------------------------------------------------------------------------------
test_run.trena <- function()
{
    message(sprintf("--- test_run.trena"))
    tbl.tms <- tms$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms,
                               abs(ampad.eqtl.score) > 1 &
                               abs(gtex.eqtl.score) > 1 & abs(cor.all) > 0.3)
    dim(tbl.tms.filtered)
    tms$set.tmsFilteredTable(tbl.tms.filtered)
    tf.candidates <- unique(tbl.tms.filtered$tf)
    if(length(tf.candidates) > 0)
        tms$run.trena(tf.candidates)
    tbl.trena <- tms$get.trenaTable()
    checkEquals(dim(tbl.trena), c(4,8))

} # test_run.trena
#----------------------------------------------------------------------------------------------------
run.small <- function(shoulder=86)
{
    tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    chrom <- "chr10"
    center <- 58385330
    shoulder <- 2000
    new.start <- center - shoulder
    new.end   <- center + shoulder
    tms$setStudyRegion(chrom=chrom, start=new.start, end=new.end)

    study.region <- tms$getStudyRegion()
    checkEquals(study.region$width.kb, 4)
    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])

    tms$run.tms()
    tms$add.eqtls.toTmsTable()

    tbl.tms <- tms$get.tmsTable()
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl.pval < 1 & gtex.eqtl.pval < 1 & abs(cor.all) > 0.4)
    tms$set.tmsFilteredTable(tbl.tms.filtered)
    tms$get.tmsFilteredTable()
    tf.candidates <- unique(tbl.tms.filtered$tf)
    tms$run.trena(tf.candidates)
    tbl.trena <- tms$get.trenaTable()
    checkTrue(all(c("SP4", "GABPA") %in% tbl.trena$gene))

    return(tms)

} # run.small
#----------------------------------------------------------------------------------------------------
# the above tests were incremental.  this one runs all steps together
test_run.small <- function()
{
    message(sprintf("--- test_run.small"))

    tms.small <- run.small()
    study.region <- tms.small$getStudyRegion()
    checkTrue(study.region$width.kb < 5)

    tbl.tms <- tms.small$get.tmsFilteredTable()
    checkEquals(dim(tbl.tms), c(6, 22))

    tbl.trena <- tms.small$get.trenaTable()
    checkEquals(dim(tbl.trena), c(5, 8))
    checkTrue(all(c("SP4", "HES7", "GABPA") %in% tbl.trena$gene))
    checkEquals(tbl.trena$gene, c("SP4", "GABPA", "PLAG1"))
    checkEquals(tbl.trena$tfbs, c(2,1,1,1,1))

} # test_run.small
#----------------------------------------------------------------------------------------------------
# the above tests were incremental.  this one runs all steps together
test_run.2k.break.motifs <- function()
{
    message(sprintf("--- test_run.2k.break.motifs"))

    tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    new.start <- 58384447
    new.end   <- 58386447
    tms$setStudyRegion(chrom="chr10", start=new.start, end=new.end)

    study.region <- tms$getStudyRegion()
    checkTrue(study.region$width.kb >= 2)
    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])

    tms$run.tms()
    tms$add.eqtls.toTmsTable()

    tbl.tms <- tms$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms,
                               abs(ampad.eqtl.score) > 1 &
                               abs(gtex.eqtl.score) > 1 & abs(cor.all) > 0.3)


    dim(tbl.tms.filtered)
    tms$set.tmsFilteredTable(tbl.tms.filtered)
    tms$get.tmsFilteredTable()
    tf.candidates <- unique(tbl.tms.filtered$tf)
    checkTrue(length(tf.candidates) == 5)
    tms$run.trena(tf.candidates)
    tbl.trena <- tms$get.trenaTable()
    checkEquals(dim(tbl.trena), c(5, 8))
    checkEquals(all(tbl.trena$gene %in% c("HES7","ZEB1","TGIF1","E2F4","USF2")))

    checkEquals(tbl.trena$tfbs, c(1,1,1,1,1))
    tms$breakMotifs(tbl.trena, tbl.tms)
    checkEquals(length(tms$get.motifBreaks()), 12)
    tbl.breaks <- tms$get.breaksTable()
    checkEquals(dim(tbl.breaks), c(5, 9))
    checkTrue("HES7" %in% tbl.breaks$geneSymbol)

} # test_run.2k.break.motifs
#----------------------------------------------------------------------------------------------------
test_run.large <- function()
{
    message(sprintf("--- test_run.large"))

    tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    study.region <- tms$getStudyRegion()
    checkTrue(study.region$width.kb > 2000)  # 2M
    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])

    tms$run.tms()
    tms$add.eqtls.toTmsTable()

    tbl.tms <- tms$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    dim(tbl.tms.filtered)
    tms$set.tmsFilteredTable(tbl.tms.filtered)
    tf.candidates <- unique(tbl.tms.filtered$tf)
    length(tf.candidates)
    checkTrue(length(tf.candidates) > 30)
    tms$run.trena(tf.candidates)
    tbl.trena <- tms$get.trenaTable()
    checkEquals(dim(tbl.trena), c(42, 8))
    checkTrue(all(tbl.trena$tfbs > 0))

    tbl.trena.strong <- subset(tbl.trena, abs(betaLasso) > .05 | rfNorm > 0.3)
    rownames(tbl.trena.strong) <- NULL
    tbl.tfbs <- subset(tbl.tms, tf %in% tbl.trena.strong$gene)
    save(tbl.trena.strong, tbl.tfbs, file="full-run-large.RData")

} # test_run.large
#----------------------------------------------------------------------------------------------------
test_viz <- function()
{
   message(sprintf("--- test_viz"))

   ts <- run.small(shoulder=300)

   tbl.trena <- ts$get.trenaTable()
   tbl.tms   <- ts$get.tmsFilteredTable()
   ts$breakMotifs(tbl.trena, tbl.tms)
   tbl.breaks <- ts$get.breaksTable()

   gc <- sort(intersect(tbl.tms$tf, tbl.breaks$geneSymbol))
   tbls <- list()
   for(gene in gc){
       tbl.breaks.gene <- subset(tbl.breaks, geneSymbol==gene)
       gr.breaks <- GRanges(tbl.breaks.gene)
       tbl.tms.gene <- subset(tbl.tms, tf==gene)
       gr.tms <- GRanges(tbl.tms.gene)
       tbl.ov <- as.data.frame(findOverlaps(gr.breaks, gr.tms))
       tbl.new <- data.frame()
       if(nrow(tbl.ov) > 0){
          tbl.new <- tbl.breaks.gene[unique(tbl.ov[,1]),]
          dups <- which(duplicated(tbl.new[, c("chrom", "start", "SNP_id", "geneSymbol")]))
          if(length(dups) > 0)
              tbl.new <- tbl.new[-dups,]
          }  # if tbl.ov
       tbls[[gene]] <- tbl.new
       browser(); xyz <-99
       }
   tbl.breaks.tfbs <- do.call(rbind, tbls)
   rownames(tbl.breaks.tfbs) <- NULL


   igv <- start.igv(targetGene, "hg38")
   ts$viz(igv)


} # test_viz
#----------------------------------------------------------------------------------------------------
test_ZCWPW1 <- function()
{
   message(sprintf("--- test_ZCWPW1"))

   targetGene <- "ZCWPW1"
   trenaProject <- TrenaProjectAD()

   data.dir <- "~/github/tms-makeAndBreak/studies/rs7384878/shared"
   tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.PILRA.RData")))
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.001)
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.rosmap.eqtls.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   tbl.rosmap.eqtls <- subset(tbl.rosmap.eqtls.raw, pvalue < 0.001)
   dim(tbl.rosmap.eqtls)   # 6262 12
   tbl.oc <- data.frame() #tbl.mayoAtac

   tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.rosmap.eqtls, tbl.oc)

   roi <- list(chrom="chr7", start=100088993, end=100600689, width=511697,
               string="chr7:100,088,993-100,600,689")

   tms$setStudyRegion(chrom=roi$chrom, start=roi$start, end=roi$end)

   study.region <- tms$getStudyRegion()
   checkTrue(study.region$width.kb > 500)
   tissues <- tms$getGTEx.eqtl.tissues()
   tms$set.current.GTEx.eqtl.tissue(tissues[1])

   tms$run.tms()
   tms$add.eqtls.toTmsTable()

   tbl.tms <- tms$get.tmsTable()
   pval.threshold <- 1e-6
   tbl.tms.filtered <-
      subset(tbl.tms, ampad.eqtl.pval < 1 & gtex.eqtl.pval < pval.threshold & abs(cor.all) > 0.4)
   dim(tbl.tms.filtered)
   length(unique(tbl.tms.filtered$tf))
   tms$set.tmsFilteredTable(tbl.tms.filtered)
   # tms$get.tmsFilteredTable()
   tf.candidates <- unique(tbl.tms.filtered$tf)
   print(length(tf.candidates))
   tms$run.trena(tf.candidates)
   tbl.trena <- tms$get.trenaTable()
   head(tbl.trena, n=10)

   tms$breakMotifs(tbl.trena, tbl.tms.filtered, tbl.rosmap.eqtls)
   tbl.breaks <- tms$get.breaksTable()

} # test_ZCWPW1
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
