library(TrenaProjectAD)
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
   tbl.ampad.eqtls.raw <- get(load(file.path(data.dir, "ampad.eqtls.pilra.plusMinus.1M.RData")))
   tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, study=="ampad-rosmap" & pvalue < 0.001)
   tbl.gtex.eqtls.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.001)
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "mayoAllPeaks.merged.96064x4.RData"
   tbl.mayoAtac <- get(load(file.path(data.dir, filename)))
    #   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    # filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.oc <- tbl.mayoAtac
   }

tms <- tmsMB$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)

chrom <- "chr7"
center <- 100334426  # tag snp
shoulder <- 100000
new.start <- center - shoulder
new.end   <- center + shoulder
tms$setStudyRegion(chrom=chrom, start=new.start, end=new.end)

study.region <- tms$getStudyRegion()
tissues <- tms$getGTEx.eqtl.tissues()
tms$set.current.GTEx.eqtl.tissue(tissues[1])

tms$run.tms()
tms$add.eqtls.toTmsTable()

tbl.tms <- tms$get.tmsTable()
tbl.tms.filtered <- subset(tbl.tms, (ampad.eqtl & gtex.eqtl) & abs(cor.all) > 0.4)
dim(tbl.tms.filtered)
tms$set.tmsFilteredTable(tbl.tms.filtered)
x <- tms$get.tmsFilteredTable()
tf.candidates <- unique(tbl.tms.filtered$tf)
tms$run.trena(tf.candidates)
tbl.trena <- tms$get.trenaTable()
head(tbl.trena, n=20)

tms$breakMotifs(tbl.trena[1:10,], tbl.tms)
breaks <- tms$get.motifBreaks()
tbl.breaks <- tms$get.breaksTable()

filename <- sprintf("%s-results.RData", targetGene)
save(tbl.trena, tbl.tms, tbl.tms.filtered, breaks, tbl.breaks, file=filename)


