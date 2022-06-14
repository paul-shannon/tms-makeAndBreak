rsid.oi <- "rs7384878"
igv <- start.igv("all", "hg38")

if(!exists("igv"))
   igv <- start.igv("all")

if(!exists("tbl.haploreg")){
   tbl.haploreg <- read.table("haploreg.tsv", sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
   dim(tbl.haploreg)
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
   }

if(!exists("tbl.rsids")){  # with mayo & rosmap scores for association of each variant to extreme braak score
   data.dir <- "~/github/gwasExplorer/explore/newScore/doesLinkagePredictBraakscSeparatedEnrichment"
   filename <- "tbl.rsids.summary.Tue.May.31.2022-13:26:55-hetAndHomSeparate.RData"
   full.path <- file.path(data.dir, filename)
   stopifnot(file.exists(full.path))
   tbl.rsids <- get(load(full.path))
   }

if(!exists("tbl.gwascat.ad")){
   data.dir <- system.file(package="igvR", "extdata", "gwas")
   file.exists(data.dir)
   file <- "alzheimerSubsetOfGWASCatatalog-29apr2022.RData"
   full.path <- file.path(data.dir, file)
   file.exists(full.path)
   tbl.gwascat.ad <- get(load(full.path))
   stopifnot(colnames(tbl.gwascat.ad)[c(1,2,33)] == c("seqnames", "start", "P.VALUE"))
   dim(tbl.gwascat.ad) # [1] 44 43
}
#----------------------------------------------------------------------------------------------------
pilra.neighborhood <- function()
{
    showGenomicRegion(igv, "PILRA")
    lapply(list(1:4), zoomOut(igv))

} # pilra.neighborhood
#----------------------------------------------------------------------------------------------------
some.standard.tracks <- function()
{
   tbl.track <- data.frame(chrom="chr7", start=100334425, end=100334426, stringsAsFactors=FALSE,
                           trackHeight=25)
   track <- DataFrameAnnotationTrack(rsid.oi, tbl.track, color="red")
   displayTrack(igv, track)

   tbl.track <- tbl.haploreg[, c("chrom", "hg38", "hg38", "Rsquared", "rsid")]
   colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
   tbl.track$chrom <- paste0("chr", tbl.track$chrom)
   tbl.track$start <- tbl.track$start - 1
   track <- DataFrameQuantitativeTrack("haploreg R^2", tbl.track, autoscale=FALSE, min=0, max=1, color="red")
   displayTrack(igv, track)

} # some.standard.tracks
#----------------------------------------------------------------------------------------------------
# results saved ~/github/slideSets/AD/non-codingVariantStudies/rs7384878-pilra.key
show.gtex.eqtls <- function()
{
    # see the 10 most frequent
  tbl.top <- head(as.data.frame(sort(table(subset(tbl.eqtl, pvalue < 1e-10)$gene), decreasing=TRUE)), n=10)
  genes <- as.character(tbl.top$Var1)

  tbl.eqtl <- get(load("gtex-eqtls-tbl.RData"))
  dim(tbl.eqtl)
  tissues <- head(unique(tbl.eqtl$id))
  for(this.gene in genes[3:10]){
    for(tissue in tissues){
       tbl.sub <- subset(tbl.eqtl, id==tissue & pvalue < 1e-10 & gene==this.gene)
       dim(tbl.sub)
       title <- sprintf("%s-%s", this.gene, tissue)
       title <- sub("GTEx_V8.Brain_", "", title, fixed=TRUE)
       track <- GWASTrack(title, tbl.sub, chrom.col=7, pos.col=8, pval.col=2)
       displayTrack(igv, track)
    } # for tissue
    browser(); xyz <- 99
  } # for this.gene

} # show.gtex.eqtls
#----------------------------------------------------------------------------------------------------

