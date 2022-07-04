library(R6)
library(EndophenotypeExplorer)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(BiocParallel)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#----------------------------------------------------------------------------------------------------
#' @title tfam
#' @description A template for building documented, tested R6 classes
#' @name tfam

#' @field id identifier for a class object
#'
#' @examples
#'   rt <- R6Template$new(id="abc")
#' @export

tmsMB = R6Class("tmsMB",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   chromosome=NA,       # only sometimes needed
                   trenaProject=NULL,
                   study.region=NULL,
                   tbl.fimo=NULL,
                   tbl.oc=NULL,
                   tbl.ampad.eqtls=NULL,
                   tbl.gtex.eqtls=NULL,
                   tbl.trena=NULL,
                   gtex.eqtl.tissues=NULL,
                   current.tissue=NULL,
                   tbl.rosmap.eqtls=NULL,
                   tms=NULL,
                   tbl.tms=NULL,
                   tbl.tmsFiltered=NULL,   # has just the TFs and regions used by trena
                   etx=NULL,
                   mtx.rna=NULL,
                   known.snps=GRanges(),
                   motifBreaks=NULL,
                   tbl.breaks=data.frame()
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param id character, an indentifier for this object
         #' @return a new instance of the class
        initialize = function(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls,
                              tbl.oc, known.snps, chromosome=NA){
            private$targetGene <- targetGene
            private$trenaProject <- trenaProject
            private$tbl.fimo <- tbl.fimo
            private$tbl.gtex.eqtls <- tbl.gtex.eqtls
            private$tbl.ampad.eqtls <- tbl.ampad.eqtls
            private$gtex.eqtl.tissues <- sort(unique(tbl.gtex.eqtls$id))
            private$current.tissue <- private$gtex.eqtl.tissues[1] # default, likely changed
            private$tbl.oc <- tbl.oc
            private$study.region <- self$getFimoGenomicRegion()
            private$etx <- EndophenotypeExplorer$new(targetGene, "hg38", vcf.project="ADNI",
                                                     chromosome=chromosome,
                                                     initialize.snpLocs=FALSE)
            private$known.snps <- known.snps
            private$chromosome <- chromosome
            },
        #------------------------------------------------------------
        setStudyRegion = function(chrom, start, end){
            width.kb <- round(((1 + end - start)/1000), digits=2)
            private$study.region <- list(chrom=chrom, start=start, end=end, width.kb=width.kb)
            },
        #------------------------------------------------------------
        getStudyRegion = function(){
            private$study.region
            },
        #------------------------------------------------------------
        getFimoGenomicRegion = function(){
            chrom <- private$tbl.fimo$chrom[1]
            start <- min(private$tbl.fimo$start)
            end   <- max(private$tbl.fimo$end)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(chrom=chrom, start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        get.gtex.eqtls = function(){
            invisible(tbl.gtex.eqtls)
            },
        #------------------------------------------------------------
        get.ampad.eqtls = function(){
            invisible(tbl.ampad.eqtls)
            },
        #------------------------------------------------------------
        getGTEx.eqtl.genomicRegion = function(){
            start <- min(tbl.gtex.eqtls$hg38)
            end   <- max(tbl.gtex.eqtls$hg38)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        getAMPAD.eqtl.genomicRegion = function(){
            start <- min(private$tbl.fimo$start)
            end   <- max(private$tbl.fimo$end)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        getGTEx.eqtl.tissues = function(){
            private$gtex.eqtl.tissues
            },
        #------------------------------------------------------------
        get.current.GTEx.eqtl.tissue = function(){
            private$current.tissue
            },
        #------------------------------------------------------------
        set.current.GTEx.eqtl.tissue = function(new.tissue){
            stopifnot(new.tissue %in% private$gtex.eqtl.tissues)
            private$current.tissue <- new.tissue
            },
        #------------------------------------------------------------
        run.tms = function(){
            current.region <- self$getStudyRegion()
            tbl.fimo.sub <- subset(tbl.fimo, start >= current.region$start & end <= current.region$end)
            private$tms <- TMS$new(private$trenaProject,
                                   private$targetGene,
                                   tbl.fimo.sub,
                                   private$tbl.oc,
                                   quiet=FALSE)

            private$tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
            current.tissue <- self$get.current.GTEx.eqtl.tissue()
            private$mtx.rna <- private$etx$get.rna.matrix(current.tissue)
            private$tms$add.tf.mrna.correlations(private$mtx.rna, featureName="cor.all")
            private$tbl.tms <- private$tms$getTfTable()
            },
        #------------------------------------------------------------
        add.eqtls.toTmsTable = function(){
            tbl.tms <- private$tbl.tms

               # give default values, corrected below when overlaps are found
            pval <- rep(1, nrow(tbl.tms))
            beta <- rep(0, nrow(tbl.tms))
            score <- rep(0, nrow(tbl.tms))

            private$tbl.tms$ampad.eqtl.pval <- pval
            private$tbl.tms$ampad.eqtl.beta <- beta
            private$tbl.tms$ampad.eqtl.scorel <- score

            private$tbl.tms$gtex.eqtl.pval <- pval
            private$tbl.tms$gtex.eqtl.beta <- beta
            private$tbl.tms$gtex.eqtl.score <- score

            gr.tms <- GRanges(tbl.tms)
            tbl.ae <- subset(private$tbl.ampad.eqtls, hg38 > private$study.region$start &
                                                      hg38 < private$study.region$end)
            if(nrow(tbl.ae) == 0) return()
            gr.ampad <- with(tbl.ae, GRanges(seqnames=chrom[1], IRanges(hg38)))
            tbl.ge <- subset(private$tbl.gtex.eqtls, id==private$current.tissue &
                                                     hg38 > private$study.region$start &
                                                     hg38 < private$study.region$end)
            if(nrow(tbl.ge) == 0) return()
            #gr.gtex <- with(tbl.ge, GRanges(seqnames=sprintf("chr%s", chrom[1]), IRanges(hg38)))
            gr.gtex <- with(tbl.ge, GRanges(seqnames=chrom[1], IRanges(hg38)))

            tbl.ov <- as.data.frame(findOverlaps(gr.ampad, gr.tms))
            if(nrow(tbl.ov) == 0) return()
                #------------------------------------
                # first the ampad eqtl pvalue & beta
                #------------------------------------

            pval <- rep(1, nrow(tbl.tms))
            beta <- rep(0, nrow(tbl.tms))
            if(nrow(tbl.ov) > 0){
               tms.indices <- tbl.ov[,2]
               eqtl.indices <- tbl.ov[,1]
               pval[tms.indices] <- tbl.ae[eqtl.indices, "pvalue"]
               beta[tms.indices] <- tbl.ae[eqtl.indices, "beta"]
               }
            private$tbl.tms$ampad.eqtl.pval <- pval
            private$tbl.tms$ampad.eqtl.beta <- beta
            private$tbl.tms$ampad.eqtl.score <- -log10(pval) * beta

                #---------------------------------
                # now the gtex eqtl pvalue & beta
                #---------------------------------

            tbl.ov <- as.data.frame(findOverlaps(gr.gtex, gr.tms))
            dim(tbl.ov)
            pval <- rep(1, nrow(tbl.tms))
            beta <- rep(0, nrow(tbl.tms))
            if(nrow(tbl.ov) > 0){
               tms.indices <- tbl.ov[,2]
               eqtl.indices <- tbl.ov[,1]
               pval[tms.indices] <- tbl.ge[eqtl.indices, "pvalue"]
               beta[tms.indices] <- tbl.ge[eqtl.indices, "beta"]
               }
            private$tbl.tms$gtex.eqtl.pval <- pval
            private$tbl.tms$gtex.eqtl.beta <- beta
            private$tbl.tms$gtex.eqtl.score <- -log10(pval) * beta
            },
        #------------------------------------------------------------
        get.tmsTable = function(){
            private$tbl.tms
            },
        #------------------------------------------------------------
        get.knownSnps = function(){
            private$known.snps
            },

        #------------------------------------------------------------
        get.motifBreaks = function(){
            private$motifBreaks
            },

        #------------------------------------------------------------
        get.breaksTable = function(){
            private$tbl.breaks
            },
        #------------------------------------------------------------
        run.trena = function(tf.candidates){
            solvers=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest")
            targetGene <- private$targetGene
            #if(targetGene == "C2orf49-DT")
            #    targetGene <- "RP11-332H14.2"  # gtex needs this gene symbol
            solver <- EnsembleSolver(private$mtx.rna,
                                     targetGene=targetGene,
                                     candidateRegulators=tf.candidates,
                                     geneCutoff=1.0,
                                     solverNames=solvers)
            tbl.trena <- trena::run(solver)
            new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
            tbl.trena <- tbl.trena[new.order,]
            rownames(tbl.trena) <- NULL
            tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
            tbl.trena <- cbind(tbl.trena[1], as.data.frame(lapply(tbl.trena[-1], function(x) round(x, digits=2))))
            tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene,
                                            function(tf) length(grep(tf, private$tbl.tmsFiltered$tf))))
            private$tbl.trena <- tbl.trena

            },
        #------------------------------------------------------------
        get.trenaTable = function(){
           private$tbl.trena
           },
        #------------------------------------------------------------
        set.tmsFilteredTable = function(tbl.filtered){
           private$tbl.tmsFiltered <- tbl.filtered
           },
        #------------------------------------------------------------
        get.tmsFilteredTable = function(){
           private$tbl.tmsFiltered
           },
        #------------------------------------------------------------
        breakMotifs = function(tbl.trena, tbl.tms, tbl.eqtls){
            #browser()
            tbl.tfbs <- subset(tbl.tms, tf %in% tbl.trena$gene)
            dim(tbl.tfbs)

            gr.tfbs <- GRanges(tbl.tfbs[, c("chrom", "start", "end")])
            #tbl.eqtls <- self$get.ampad.eqtls()
            tbl.eqtls.tmp <- tbl.eqtls[, c("chrom", "hg38", "hg38", "rsid")]
            colnames(tbl.eqtls.tmp) <- c("chrom", "start", "end", "rsid")
            #tbl.eqtls.tmp$chrom <- sprintf("chr%s", tbl.eqtls.tmp$chrom)
            gr.eqtls <- GRanges(tbl.eqtls.tmp)
            tbl.ov <- as.data.frame(findOverlaps(gr.eqtls, gr.tfbs))
            dim(tbl.ov)
            tbl.eqtls.ov <- tbl.eqtls.tmp[tbl.ov$queryHits,]
            tbl.tfbs.ov <- tbl.tfbs[tbl.ov$subjectHits,]
            tfs.oi <- unique(tbl.tfbs.ov$tf)
            tfs.oi <- tfs.oi[match(tbl.trena$gene, tfs.oi)]
            deleters <- which(is.na(tfs.oi))
            if(length(deleters) > 0)
                tfs.oi <- tfs.oi[-deleters]
            rsids.oi <- unique(tbl.eqtls.ov$rsid)
            length(rsids.oi)

            #mdb.human <- query(MotifDb, "sapiens", c("jaspar2022", "hocomoco-core-A"))
            mdb.human <- query(MotifDb, c("sapiens", "jaspar2022"))
            motifs.selected <- query(mdb.human, andStrings="sapiens", orStrings=tfs.oi)
            message(sprintf("--- about to look up %d rsids", length(rsids.oi)))
                # get fast reproducible results by doing a throw-away call beforehand:
                #   t0 <- system.time(x0 <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, "rs7796006"))
            lookup.rsid <- function(rsid){
                snps.from.rsid(rsid,
                               dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                               search.genome=BSgenome.Hsapiens.UCSC.hg38)
                }
            if(length(private$known.snps) == 0){
               new.rsids <- rsids.oi
            } else{
                new.rsids <- setdiff(rsids.oi, names(private$known.snps))
                }
            printf("--- looking up %d new snps, with %d already known",
                   length(new.rsids), length(private$known.snps))
            t2 <- system.time({
                snps.gr.list <-lapply(new.rsids, lookup.rsid)
                snps.gr <- unlist(as(snps.gr.list, "GRangesList"))
                })
            private$known.snps <- c(private$known.snps, snps.gr)
            if(length(private$known.snps) > 0)
               snps.gr <- subset(private$known.snps, SNP_id %in% rsids.oi)
            attributes(snps.gr)$genome.package <- attributes(BSgenome.Hsapiens.UCSC.hg38)$pkgname
            message(sprintf("--- snps.gr for %d snps obtained, elapsed time: %f",
                            length(rsids.oi), t2[["elapsed"]]))

            bpparam <- MulticoreParam(workers=20)
            message(sprintf("--- calling motifbreakR"))

            results <- motifbreakR(snpList = snps.gr,
                                   filterp = TRUE,
                                   pwmList = motifs.selected,
                                   show.neutral=FALSE,
                                   method = c("ic", "log", "notrans")[1],
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = bpparam,
                                   verbose=TRUE)
            private$motifBreaks <- results
            tbl.breaks.tfbs <- data.frame()
            if(length(results) > 0){
               tbl.breaks <- as.data.frame(results, row.names=NULL)
               colnames(tbl.breaks)[1] <- "chrom"
               tbl.breaks$chrom <- as.character(tbl.breaks$chrom)
               tbl.breaks <- subset(tbl.breaks, effect=="strong")
               tbl.breaks$start <- tbl.breaks$start - 1
               tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt-pctRef)
               new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
               tbl.breaks <- tbl.breaks[new.order,]
               motifbreakR.coi <- c("chrom", "start", "end", "SNP_id", "geneSymbol", "providerName", "pctRef", "pctAlt", "pctDelta")
               tbl.breaks <- tbl.breaks[, motifbreakR.coi]

                # motifbreakR will find fresh binding sites, possibly for genes for which we did not
                # so traverse these break regions, make sure they overlap with regions we found for
                # those genes
               gc <- sort(intersect(tbl.tms$tf, tbl.breaks$geneSymbol))
               tbls <- list()
               for(gene in gc){
                   tbl.breaks.gene <- subset(tbl.breaks, geneSymbol==gene)
                   gr.breaks <- GRanges(tbl.breaks.gene)
                   tbl.tms.gene <- subset(tbl.tms, tf==gene)
                   gr.tms <- GRanges(tbl.tms.gene)
                   tbl.ov <- as.data.frame(findOverlaps(gr.breaks, gr.tms))
                   tbl.new <- data.frame()
                   message(sprintf("--- found %d breaks in tfbs for %s", nrow(tbl.ov), gene))
                   if(nrow(tbl.ov) > 0){
                       tbl.new <- tbl.breaks.gene[unique(tbl.ov[,1]),]
                       dups <- which(duplicated(tbl.new[, c("chrom", "start", "SNP_id", "geneSymbol")]))
                       if(length(dups) > 0)
                           tbl.new <- tbl.new[-dups,]
                   }  # if tbl.ov
                   tbls[[gene]] <- tbl.new
               } # for gene
               tbl.breaks.tfbs <- do.call(rbind, tbls)
               rownames(tbl.breaks.tfbs) <- NULL
               } # if results > 0
            message(sprintf("--- saving tbl.breaks, %d x %d", nrow(tbl.breaks.tfbs), ncol(tbl.breaks.tfbs)))
            private$tbl.breaks <- tbl.breaks.tfbs
            },  # breakMotifs
        #------------------------------------------------------------
        viz = function(igv, current.tissue, tbl.trena, tbl.tms, tbl.tms.filtered,
                       tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc, tbl.breaks) {
            shoulder <- 1000
            roi <- self$getStudyRegion()
            roi.string <- with(roi, sprintf("%s:%d-%d", chrom, start-1000, end+1000))
            showGenomicRegion(igv, roi.string)

                #----------------------------------------------------
                # motifbreakR results
                #-----------------------------------------------------

            if(nrow(tbl.breaks) > 0){
                tbl.breaks$score <- with(tbl.breaks, pctRef * pctDelta * 100)
                dups <- which(duplicated(tbl.breaks[, c("SNP_id", "geneSymbol")]))
                if(length(dups) > 0)
                    tbl.breaks <- tbl.breaks[-dups,]
                tbl.track <- tbl.breaks[, c("chrom", "start", "end", "score", "SNP_id", "geneSymbol")]
                track <- DataFrameQuantitativeTrack("motif disruption", tbl.track, colo="red", autoscale=TRUE)
                displayTrack(igv, track)
                }  # tbl.breaks

                #----------------------------------------------------
                # gtex tissue-specific eqtls, alredy filterd for pval
                #-----------------------------------------------------

            tbl.eqtl <- tbl.gtex.eqtls
            tbl.eqtl <- subset(tbl.eqtl, hg38 >= roi$start & hg38 <= roi$end)
            tbl.eqtl <- tbl.eqtl[, c("chrom", "hg38", "hg38", "rsid", "pvalue", "beta")]
            colnames(tbl.eqtl)[c(2,3)] <- c("start", "end")
            tbl.eqtl$start <- tbl.eqtl$start - 1
            tbl.eqtl$score <- with(tbl.eqtl, -log10(pvalue) * beta)
            coi <- c("chrom", "start", "end", "score", "rsid")
            tbl.eqtl <- tbl.eqtl[, coi]
            dups <- which(duplicated(tbl.eqtl[, c("chrom", "start", "end", "rsid")]))
            if(length(dups) > 0)
                tbl.eqtl <- tbl.eqtl[-dups,]
            title <- current.tissue
            title <- sub("Brain_", "", title)
            track <- DataFrameQuantitativeTrack(title, tbl.eqtl, autoscale=TRUE, color="black")
            displayTrack(igv, track)

                #----------------------------------------------------
                # ampad rosmap eqtls, alredy filterd for pval
                #-----------------------------------------------------
            tbl.track <- subset(tbl.ampad.eqtls, genesymbol==private$targetGene)
            track <- GWASTrack("rosmap eqtls", tbl.track, chrom.col=1, pos.col=3, pval.col=5)
            displayTrack(igv, track)


            tbl.tms <- tbl.tms.filtered
            stopifnot(all(tbl.tms$tf %in% tbl.trena$gene))
                #-----------------------------------
                # draw the tfbs, one track per tf
                #-----------------------------------
            for(TF in tbl.trena$gene){
                tbl.track <- subset(tbl.tms, tf==TF)[, c("chrom", "start", "end")]
                track <- DataFrameAnnotationTrack(TF, tbl.track, color="random", trackHeight=25)
                displayTrack(igv, track)
                }
                #-----------------------------------
                # draw the open chromatin
                #-----------------------------------
            extra.margin <- 2000
            tbl.track <- subset(tbl.oc, chrom==roi$chrom &
                                        start >= roi$start-extra.margin &
                                        end <= roi$end+extra.margin)[, c("chrom", "start", "end")]
            track <- DataFrameAnnotationTrack("OC", tbl.track, color="darkgrey")
            displayTrack(igv, track)



            } # viz
        #------------------------------------------------------------
       ) # public

    ) # class

#----------------------------------------------------------------------------------------------------
# some genes not yet found in TxDb.Hsapiens.UCSC.hg38.knownGene, 3.15.0.  lncRNA C2orf49-DT
# for instance.  so chromosome must be supplied.
checkGeneNameFoundInAllSources <- function(targetGene, chromosome=NA)
{
      # TrenaProjectAD on top of TrenaProjectHG38, transcripts table providing chrom start end

   path <- system.file(package="TrenaProjectHG38", "extdata", "geneInfoTable.RData")
   #path <- "~/github/TrenaProjectHG38/inst/extdata/geneInfoTable.RData"
   tbl.hg38 <- get(load(path))
   printf("hg38 transcripts: %s", targetGene %in% tbl.hg38$geneSymbol)

     # GeneHancer

   require(ghdb)
   ghdb <- GeneHancerDB()
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
   printf("GeneHancer: %s", nrow(tbl.gh) > 0)

    # ampad (rosmap) eQTLs
   data.dir <- "../shared"
   full.path <- file.path(data.dir, "tbl.eqtls.rosmap.RData")
   file.exists(full.path)
   tbl.eqtl.rosmap.raw <- get(load(file.path(data.dir, "tbl.eqtls.rosmap.RData")))
   printf("ampad (ROSMAP) eQTLs: %s", targetGene %in% tbl.eqtl.rosmap.raw$genesymbol)

    # GTEx eQTLs
   data.dir <- "../shared"
   tbl.eqtl.gtex.raw <- get(load(file.path(data.dir, "gtex-eqtls-tbl.RData")))
   printf("GTEx eQTLs: %s", targetGene %in% tbl.eqtl.gtex.raw$gene)

    # gtex expression matrices
    etx <- EndophenotypeExplorer$new(targetGene, "hg38", vcf.project="AMPAD",
                                     chromosome=chromosome)
    codes.full <- etx$get.rna.matrix.codes()
    codes <- names(codes.full)
    checkTrue(length(codes) >= 10)
    for(code in codes){
       if(!grepl("GTEx_V8.Brain", code)) next;
       mtx <- etx$get.rna.matrix(code)
       printf("expression matrix %s: %s", code, targetGene %in% rownames(mtx))
       }



} # checkGeneNameFoundInAllSources
#----------------------------------------------------------------------------------------------------
