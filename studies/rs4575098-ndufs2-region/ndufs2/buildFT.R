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
                   etx=NULL,
                   tbl.filtered=NULL,
                   tbl.filtered.collapsed=NULL,
                   tbl.trena=NULL
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
       getRnaMatrixCode = function(){
           names(private$etx$get.rna.matrix.codes())
           },

       #------------------------------------------------------------
       setFilteredTable = function(tbl.filtered){
           private$tbl.filtered <- tbl.filtered
           },

       #------------------------------------------------------------
       getFilteredTable = function(){
           invisible(private$tbl.filtered)
           },

       #------------------------------------------------------------
       getCollapsedFilteredTable = function(){
           invisible(private$tbl.filtered.collapsed)
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
           }, # add.genehancer

       #------------------------------------------------------------
       run.trena = function(tfs){
          solver <- EnsembleSolver(private$etx$get.rna.matrix(private$gtex.tissue),
                                   targetGene=private$targetGene,
                                   candidateRegulators=tfs,
                                   solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
                                   geneCutoff=1.0)
          tbl.trena <- run(solver)
          new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
          tbl.trena <- tbl.trena[new.order,]
          tbl.fft <-  self$getCollapsedFilteredTable()
          tbl.trena$tfbs <- as.integer(lapply(tbl.trena$gene, function(gene) nrow(subset(tbl.fft, tf==gene))))
          gtex.eqtl.column.name <- sub("V8.Brain_", "", sprintf("%s.eqtl.rsid", private$gtex.tissue))
          x <- lapply(tbl.trena$gene, function(gene) subset(tbl.fft, tf==gene)[, ..gtex.eqtl.column.name])
          rsids <- as.character(lapply(x, function(rsids) paste(unlist(rsids), collapse=";")))
          tbl.trena$rsids <- rsids
          rownames(tbl.trena) <- NULL
          private$tbl.trena <- tbl.trena
          invisible(private$tbl.trena)
          }, # run.trena

       #------------------------------------------------------------
          # if a tf has overlapping binding sites, and the same rsid, reduce it to a single tfbs
       reduce.tbl.filtered = function(){
          stopifnot(!is.null(private$tbl.filtered))
          tbl.fft <- private$tbl.filtered
          tf.xtab <- as.list(table(tbl.fft$tf))
          tfs.mult <- names(tf.xtab[tf.xtab > 1])
          tfs.single <- names(tf.xtab[tf.xtab == 1])
          tbl.fftC <- subset(tbl.fft, tf %in% tfs.single)  # "C":  collapsed
          gtex.eqtl.column.name <- sub("V8.Brain_", "", sprintf("%s.eqtl.rsid", private$gtex.tissue))
          for(tf.mult in tfs.mult){
              tbl.sub <- subset(tbl.fft, tf==tf.mult)
              all.rsids <- as.character(unlist(tbl.sub[, ..gtex.eqtl.column.name]))
                 # match returns only the first hit
              keeper.rows <- match(unique(all.rsids), all.rsids)
              tbl.sub.collapsed <- tbl.sub[keeper.rows,]
              printf("adding %d rows for %s", nrow(tbl.sub.collapsed), tf.mult)
              tbl.fftC <- rbind(tbl.fftC, tbl.sub.collapsed)
              }
          #browser()
          private$tbl.filtered.collapsed <- tbl.fftC
          xyz <- 99
          },
       findMotifBreaks = function(){
           snps.gr <- snps.from.rsid(rsid=rsids.sig,
                                     dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                                     search.genome=BSgenome.Hsapiens.UCSC.hg38)

           motifs <- query(MotifDb, c("sapiens"), c("jaspar2018", "HOCOMOCOv11-core-A"))

           results <- motifbreakR(snpList = snps.gr,
                                  filterp = TRUE,
                                  pwmList = motifs,
                                  show.neutral=FALSE,
                                  method = c("ic", "log", "notrans")[1],
                                  bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                  BPPARAM = BiocParallel::bpparam(),
                                  verbose=TRUE)
           tbl.breaks <- as.data.frame(results, row.names=NULL)
           tbl.breaks <- subset(tbl.breaks, effect=="strong")
           } # findMotifBreaks

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
                   gtex.tissue="GTEx_V8.Brain_Hippocampus",
                   #gtex.tissue="GTEx_V8.Brain_Frontal_Cortex_BA9",
                   fimo.file="../shared/tbl.fimo.NDUFS2.RData")


go <- function(){
    bft$add.rosmap.eqtls()
    bft$add.genehancer("brain")
    bft$add.genehancer("all")
    bft$add.gtex.eqtls()
    bft$add.mayo.atac()
    bft$add.boca.atac()
    bft$add.gtex.expression.correlation()
    tbl <- bft$getTable()
    tbl.filtered <- subset(tbl, abs(rna.cor.Hippocampus) > 0.3 &
                                abs(GTEx_Hippocampus.eqtl.score > 0.4))
                                # (mayoAtac | bocaAtac))
    table(tbl.filtered$tf)
    bft$setFilteredTable(tbl.filtered)
    bft$reduce.tbl.filtered()
    tbl.fftC <- bft$getCollapsedFilteredTable()
    table(tbl.fftC$tf)
    tfs <- sort(unique(tbl.fftC$tf))

    tbl.trena <- bft$run.trena(tfs)
    dim(tbl.trena)
    }
