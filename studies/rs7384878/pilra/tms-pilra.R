library(TrenaProjectAD)
library(ghdb)


trenaProject <- TrenaProjectAD()
f <- "~/github/tms-makeAndBreak/R/tmsMB-class.R"
stopifnot(file.exists(f))
source(f)
data.dir <- "../shared"
stopifnot(file.exists(data.dir))
targetGene <- "PILRA"
trenaProject <- TrenaProjectAD()


if(!exists("tbl.fimo")){
   tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.PILRA.RData")))
   #tbl.ampad.eqtls.raw <- get(load(file.path(data.dir, "ampad.eqtls.pilra.plusMinus.1M.RData")))
   #tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, study=="ampad-rosmap" & pvalue < 0.001)
   tbl.ampad.eqtls.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, genesymbol==targetGene & pvalue < 0.001)
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.001)
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "mayoAllPeaks.merged.96064x4.RData"
   tbl.mayoAtac <- get(load(file.path(data.dir, filename)))
    #   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    # filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.oc <- tbl.mayoAtac
   tbl.haploreg <- read.table("../shared/haploreg.tsv", sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
   }

tag.snp <- "rs7384878"
tag.snp.chrom <- "chr7"
tag.snp.loc <- 100334426
#----------------------------------------------------------------------------------------------------
run <- function()
{
    browser()
    chrom <- tag.snp.chrom
    center <- tag.snp.loc
    shoulder <- 100000 #000
    #roi <- getGenomicRegion(igv)
    new.start <- center - shoulder
    new.end   <- center + shoulder

    tms$setStudyRegion(chrom=chrom, start=new.start, end=new.end)
    #tms$setStudyRegion(chrom=roi$chrom, start=roi$start, end=roi$end)

    study.region <- tms$getStudyRegion()
    study.region$width.kb
    tissues <- tms$getGTEx.eqtl.tissues()
    tms$set.current.GTEx.eqtl.tissue(tissues[1])  # "GTEx_V8.Brain_Frontal_Cortex_BA9"

    tms$run.tms()
    tms$add.eqtls.toTmsTable()
    tbl.tms <- tms$get.tmsTable()
    dim(tbl.tms)

     psql.method <- function(){
        require(RPostgreSQL)
        roi <- getGenomicRegion(igv)
        db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
        dbGetQuery(db, "select * from eqtl2 limit 3")
        query <- sprintf("select * from eqtl2 where chrom='chr7' and hg38 > %d and hg38 < %d",
                         roi$start, roi$end)
        tbl.db <- dbGetQuery(db, query)
        tbl.db <- subset(tbl.db, genesymbol=="PILRA")
        tbl.track <- tbl.db[, c("chrom", "hg38", "hg38", "pvalue", "beta")]
        colnames(tbl.track) <- c("chrom", "start", "end", "pvalue", "beta")
        tbl.track$start <- tbl.track$start - 1
        tbl.track$score <- with(tbl.track, -log10(pvalue) * beta)
        fivenum(tbl.track$score)
        track <- DataFrameQuantitativeTrack("rosmap eqtl score",
                                            tbl.track[, c("chrom", "start", "end", "score")],
                                            color="blue", autoscale=TRUE)
        displayTrack(igv, track)
        }

        #     select * from eqtl2 where chrom='chr7' and hg38 > 100357719 and hg38 < 100357756 and genesymbol ='PILRA';


    tbl.tms <- tms$get.tmsTable()
    eqtl.score.threshold <- 25
    tbl.tms.filtered <- subset(tbl.tms,
                               abs(ampad.eqtl.score) > eqtl.score.threshold &
                               abs(gtex.eqtl.score)  > eqtl.score.threshold &
                               abs(cor.all) > 0.3)# &
                                    #fimo_pvalue < 5e-4)

    dim(tbl.tms.filtered)
    nrow(subset(tbl.tms.filtered, tf=="TEAD1"))
    tms$set.tmsFilteredTable(tbl.tms.filtered)
    x <- tms$get.tmsFilteredTable()
    tf.candidates <- unique(tbl.tms.filtered$tf)
    length(tf.candidates)
    tms$run.trena(tf.candidates)
    tbl.trena <- tms$get.trenaTable()
    head(tbl.trena, n=20)

    tms$breakMotifs(tbl.trena[1:10,], tbl.tms)
    breaks <- tms$get.motifBreaks()
    tbl.breaks <- tms$get.breaksTable()

    timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M")))

    filename <- sprintf("%s-%s-results-%s.RData", targetGene, tms$get.current.GTEx.eqtl.tissue(), timestamp)
    save(tbl.trena, tbl.tms, tbl.tms.filtered, breaks, tbl.breaks, file=filename)

} # run
#----------------------------------------------------------------------------------------------------
viz <- function()
{
  filename <- "PILRA-GTEx_V8.Brain_Cerebellar_Hemisphere-results.RData"
  print(load(filename))
  igv <- start.igv(targetGene, "hg38")
  showGenomicRegion(igv, "chr7:99,915,322-100,857,434")

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
  track <- DataFrameQuantitativeTrack("haploreg", tbl.track, autoscale=TRUE, color="red")
  displayTrack(igv, track)

  ghdb <- GeneHancerDB()
  tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
  #tbl.gh$score <- asinh(tbl.gh$combinedscore)
  tbl.gh$score <- tbl.gh$combinedscore
  track <- DataFrameQuantitativeTrack("GH-all", tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
  displayTrack(igv, track)


  tbl.gh.brain <- tbl.gh[grep("brain", tbl.gh$tissue, ignore.case=TRUE),]
  track <- DataFrameQuantitativeTrack("GH-brain", tbl.gh.brain[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
  displayTrack(igv, track)


  tbl.trena.sub <- subset(tbl.trena, abs(betaLasso) >= 0.05 | rfScore > 10)
  tbl.tms.filtered <- subset(tbl.tms.filtered, tf %in% tbl.trena.sub$gene)

  tms$viz(igv,
          current.tissue="GTEx_V8.Brain_Cerebellar_Hemisphere",
          tbl.trena.sub,
          tbl.tms,
          tbl.tms.filtered,
          tbl.gtex.eqtls,
          tbl.ampad.eqtls,
          tbl.oc,
          tbl.breaks)

} # viz
#----------------------------------------------------------------------------------------------------
