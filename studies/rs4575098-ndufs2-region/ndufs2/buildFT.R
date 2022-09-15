library(trena)
library(ghdb)
library(RUnit)
library(RPostgreSQL)
library(EndophenotypeExplorer)
#----------------------------------------------------------------------------------------------------
tbl.fimo <- get(load("../shared/tbl.fimo.NDUFS2.RData"))
p.value.col <- grep("^p\\.value$", colnames(tbl.fimo))
if(length(p.value.col) == 1)
    colnames(tbl.fimo)[p.value.col] <- "fimo.pval"
tbl.braakSep <- get(load("../shared/braakScore.ac.ctl.separation.22816x4.RData"))
tbl.rsids <- get(load("../shared/tbl.rsid.160684126-161688825.hg38.hg19.RData"))
checkTrue(all(rownames(tbl.braakSep) %in% tbl.rsids$rsid))
#----------------------------------------------------------------------------------------------------
buildFT = R6Class("buildFT",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   gtex.tissue=NULL,
                   ft=NULL,
                   rosmap.eqtls=NULL,
                   fimo.file=NULL,
                   chrom=NULL,
                   start=NULL,
                   end=NULL,
                   etx=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

        initialize = function(targetGene, gtex.tissue, fimo.file,
                              study.region.start=NA, study.region.end=NA){
           private$targetGene <- targetGene
           private$gtex.tissue <- gtex.tissue
           private$ft <- FeatureTable$new(target.gene=targetGene, reference.genome="hg39")
           private$etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
           stopifnot(file.exists(fimo.file))
           private$fimo.file <- fimo.file
           tbl.fimo <- get(load(fimo.file))
           if(!is.na(study.region.start))
               tbl.fimo <- subset(tbl.fimo, start >= study.region.start)
           if(!is.na(study.region.end))
               tbl.fimo <- subset(tbl.fimo, end <= study.region.end)
           private$start <- min(tbl.fimo$start)
           private$end   <- max(tbl.fimo$end)
           private$chrom <- tbl.fimo$chrom[1]
           p.value.col <- grep("^p\\.value$", colnames(tbl.fimo))
           if(length(p.value.col) == 1)
               colnames(tbl.fimo)[p.value.col] <- "fimo.pval"
           score.col <- grep("^score$", colnames(tbl.fimo))
           if(length(score.col) == 1)
               colnames(tbl.fimo)[score.col] <- "fimo.score"
           private$ft$setFundamentalRegions(tbl.fimo)
           },

       #------------------------------------------------------------
       getTable = function(){
           return(private$ft$getTable())
           },

       #------------------------------------------------------------
       getRegion = function(){
           return(list(chrom=private$chrom,
                       start=private$start,
                       end=private$end,
                       width=1 + private$end - private$start))
           },

       #------------------------------------------------------------
       add.rosmap.eqtls = function(){
           if(file.exists("rosmap.eqtls.RData")){
               printf("adding rosmap eqtls from prior run")
               tbl.eqtl <- get(load("rosmap.eqtls.RData"))
           } else {
              geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
              query <- "select * from eqtl2 where chrom='chr1' and hg38 > 160684126 and hg38 < 161688825  and genesymbol ='NDUFS2';"
              tbl.eqtl <- dbGetQuery(geneRegDB, query)
              tbl.eqtl$start <- tbl.eqtl$hg38 - 1
              tbl.eqtl$end <- tbl.eqtl$hg38
              }
           if(!"score" %in% colnames(tbl.eqtl))
               tbl.eqtl$score <- with(tbl.eqtl, -log10(pvalue) * beta)
           feature.guide <- list(rosmap.eqtl.rsid="rsid", rosmap.eqtl.pvalue="pvalue",
                                 rosmap.eqtl.beta="beta", rosmap.eqtl.score="score")
           checkTrue(all(as.character(feature.guide) %in% colnames(tbl.eqtl)))
           default.values <- list(rosmap.eqtl.rsid="", rosmap.eqtl.pvalue=1,
                                  rosmap.eqtl.beta=0, rosmap.eqtl.score=0)
           private$ft$addRegionFeature(tbl.eqtl, feature.genome="hg38", feature.guide, default.values)
           private$rosmap.eqtls <- tbl.eqtl
           }, # add.rosmap.eqtls

       #------------------------------------------------------------
       add.gtex.eqtls = function(){
           gtex.eqtl.file <- "gtex.brain.eqtls.chr1-160684126-161688825.RData"
           stopifnot(file.exists(gtex.eqtl.file))
           tbl.eqtl <- get(load(gtex.eqtl.file))
           tbl.eqtl.sub <- subset(tbl.eqtl, hg38 >= private$start & hg38 <= private$end &
                                            chrom==private$chrom & id == private$gtex.tissue &
                                            gene==private$targetGene)
           tbl.eqtl.sub$start <- tbl.eqtl.sub$hg38 - 1
           tbl.eqtl.sub$end <- tbl.eqtl.sub$hg38
           if(!"score" %in% colnames(tbl.eqtl.sub))
               tbl.eqtl.sub$score <- with(tbl.eqtl.sub, -log10(pvalue) * beta)
           feature.guide <- list("rsid", "pvalue", "beta", "score")
           short.name <- sub("V8.Brain_", "", private$gtex.tissue)
           feature.names <- sprintf("%s.eqtl.%s", short.name, feature.guide)
           names(feature.guide) <- feature.names
           default.values <- list("", 1, 0, 0)
           names(default.values) <- feature.names
           private$ft$addRegionFeature(tbl.eqtl.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.gtex.eqtls

       #------------------------------------------------------------
       add.mayo.atac = function(){

           data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
           file <- "mayoAllPeaks.merged.96064x4.RData"
           full.path <- file.path(data.dir, file)
           stopifnot(file.exists(full.path))
           tbl.atac <- get(load(full.path))
           tbl.atac.sub <- subset(tbl.atac, start >= private$start & end <= private$end &
                                            chrom==private$chrom)
           tbl.atac.sub$status <- TRUE
           feature.guide <-  list(mayoAtac="status")
           default.values <- list(mayoAtac=FALSE)
           private$ft$addRegionFeature(tbl.atac.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.mayo.atac

       #------------------------------------------------------------
       add.boca.atac = function(){

           data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
           file <- "boca-hg38-consensus-ATAC.RData"
           full.path <- file.path(data.dir, file)
           stopifnot(file.exists(full.path))
           tbl.atac <- get(load(full.path))
           tbl.atac.sub <- subset(tbl.atac, start >= private$start & end <= private$end &
                                            chrom==private$chrom)
           tbl.atac.sub$status <- TRUE
           feature.guide <-  list(bocaAtac="status")
           default.values <- list(bocaAtac=FALSE)
           private$ft$addRegionFeature(tbl.atac.sub, feature.genome="hg38", feature.guide, default.values)
           }, # add.boca.atac

       #------------------------------------------------------------
       add.gtex.expression.correlation = function(){
           stopifnot(private$gtex.tissue %in% names(private$etx$get.rna.matrix.codes()))
           mtx <- private$etx$get.rna.matrix(private$gtex.tissue)
           short.name <- sub("GTEx_V8.Brain_", "", private$gtex.tissue)
           x <- lapply(rownames(mtx), function(gene) cor(mtx[private$targetGene,], mtx[gene,]))
           tbl.cor <- data.frame(gene=rownames(mtx), cor=as.numeric(x))
           new.colname <- sprintf("rna.cor.%s", short.name)
           colnames(tbl.cor)[2] <- new.colname
           private$ft$addGeneFeature(tbl.cor, feature.name=new.colname, default.value=0)
           },

       #------------------------------------------------------------
       add.genehancer = function(tissue){

           stopifnot(tissue %in% c("brain", "all"))
           ghdb <- GeneHancerDB()
           tbl.gh <- retrieveEnhancersFromDatabase(ghdb, private$targetGene, tissues=tissue)
           feature.name <- sprintf("gh.%s", tissue)
           feature.guide <- list("combinedscore")
           default.values <- list(0)
           names(feature.guide) <- feature.name
           names(default.values) <- feature.name
           private$ft$addRegionFeature(tbl.gh, feature.genome="hg38", feature.guide, default.values)
           } # add.genehancer

       ) # public


    ) # class
#----------------------------------------------------------------------------------------------------
# find a tiny region with both rosmap and gtex eqtls
choose.test.region <- function()
{
   igv <- start.igv("NDUFS2", "hg38")

   tbl.eqtl.rosmap <- get(load("rosmap.eqtls.RData"))
   tbl.eqtl.rosmap$score <- with(tbl.eqtl.rosmap, -log10(pvalue) * beta)
   dim(tbl.eqtl.rosmap)

   tbl.track <- tbl.eqtl.rosmap[, c("chrom", "start", "end", "score")]
   track <- DataFrameQuantitativeTrack("rosmap", tbl.track, color="random", autoscale=TRUE)
   displayTrack(igv, track)

   gtex.eqtl.file <- "gtex.brain.eqtls.chr1-160684126-161688825.RData"
   stopifnot(file.exists(gtex.eqtl.file))
   tbl.eqtl.gtex <- get(load(gtex.eqtl.file))
   tbl.eqtl.gtex$start <- tbl.eqtl.gtex$hg38 - 1
   tbl.eqtl.gtex$end <- tbl.eqtl.gtex$hg38
   tbl.eqtl.gtex <- subset(tbl.eqtl.gtex, gene=="NDUFS2")
   if(!"score" %in% colnames(tbl.eqtl.gtex))
      tbl.eqtl.gtex$score <- with(tbl.eqtl.gtex, -log10(pvalue) * beta)

   tbl.track <- tbl.eqtl.gtex[, c("chrom", "start", "end", "score")]
   track <- DataFrameQuantitativeTrack("gtex", tbl.track, color="random", autoscale=TRUE)
   displayTrack(igv, track)

   getGenomicRegion(igv) # $chrom  [1] "chr1"
                         # $start  [1] 161214266
                         # $end    [1] 161214354
                         # $width  [1] 89
                         # $string [1] "chr1:161,214,266-161,214,354"

} # choose.test.region
#----------------------------------------------------------------------------------------------------
targetGene <- "NDUFS2"
bft <- buildFT$new(targetGene=targetGene,
                   gtex.tissue="GTEx_V8.Brain_Frontal_Cortex_BA9",
                   fimo.file="../shared/tbl.fimo.NDUFS2.RData")

run <- function(){
    bft$add.rosmap.eqtls()
    bft$add.genehancer("brain")
    bft$add.genehancer("all")
    bft$add.gtex.eqtls()
    bft$add.mayo.atac()
    bft$add.boca.atac()
    bft$add.gtex.expression.correlation()
    tbl <- bft$getTable()
    subset(tbl, abs(GTEx_Frontal_Cortex_BA9.score) > 0.4 &
                abs(rna.cor.Frontal_Cortex_BA9) > 0.4 &
                (mayoAtac | bocaAtac))$tf
    table(subset(tbl, abs(GTEx_Frontal_Cortex_BA9.score) > 0.2 &
                      abs(rna.cor.Frontal_Cortex_BA9) > 0.4 &
#                      abs(rosmap.eqtl.score) > 0.3)$tf)
                      (mayoAtac | bocaAtac))$tf)
    table(subset(tbl,
                 abs(rna.cor.Frontal_Cortex_BA9) > 0.3 &
                 abs(GTEx_Frontal_Cortex_BA9.score) > 0.2 &
                 (mayoAtac | bocaAtac | gh > 10)
                 )$tf)
    subset(tbl, tf=="SOX21" &
                abs(GTEx_Frontal_Cortex_BA9.score) > 0.2 &
                abs(rna.cor.Frontal_Cortex_BA9) > 0.4)  # see slide 26 of slideset "rs4575098"
       #   1  SOX21    -0.073    -0.025        -0.560       -0.545  13.036   0.135    1 NDUFS2    2   1.00          15.63
       #   2   EBF1    -0.071    -0.040        -0.554       -0.564  11.814   0.317    2 NDUFS2   18   0.91           7.40
       #   3   NFIA    -0.217    -0.044        -0.590       -0.598  10.888   0.032    3 NDUFS2    3   0.84           4.49
       #   4   ZEB1    -0.066    -0.025        -0.592       -0.572   8.691   0.016    4 NDUFS2   12   0.67           4.49
       #   5  FOXL1    -0.117    -0.069        -0.488       -0.485   7.462   0.012    5 NDUFS2    1   0.57           6.13
       #   6  ASCL1     0.000    -0.028        -0.520       -0.524   6.225   0.007    6 NDUFS2   13   0.48           2.99
    }
