library(GenomicRanges)
print(load("PILRB-models-roi.528k-21jun2022.RData"))
lapply(models, function(model) head(model$trena, n=5))
lapply(models, function(model) head(model$tbl.breaks, n=5))

motifBreak.score <- with(tbl.tf, abs(pctDelta) * 100)
      eqtl.score <- with(tbl.tf, -log10(gtex.eqtl.pval)* abs(gtex.eqtl.beta) * 100)
     trena.score <- with(tbl.tf, (abs(betaLasso) * 100) + (rfNorm * 10))
      tfbs.score <- 1/tfbs.count

score <- trena.score * eqtl.score * motifBreak.score * tfbs.score

names(models)


#----------------------------------------------------------------------------------------------------
combine.tables <- function(tbl.trena, tbl.tms, tbl.breaks, targetGene, TF)
{
    browser()
    tbl.trena <- subset(tbl.trena, gene==TF)[, c("gene", "betaLasso", "spearmanCoeff", "rfScore", "rfNorm", "tfbs", "tf.rank")]
    colnames(tbl.trena) <- c("tf", "betaLasso", "spearman", "rfScore", "rfNorm", "tfbs", "tf.rank")
    tbl.tms <- subset(tbl.tms, tf==TF)  # only some of these will be broken
    tbl.breaks <- subset(tbl.breaks, geneSymbol==TF)
    colnames(tbl.breaks)[3:6] <- c("hg38", "rsid", "tf", "motif")
    tbl.breaks <- tbl.breaks[, c("chrom", "hg38", "rsid", "tf", "motif", "pctRef", "pctDelta")]

         #--------------------------------------------------
         # intersect breaks with tfbs in tbl.tms
         #--------------------------------------------------
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
    coi <- c("targetGene", "tf", "tf.rank", "rsid", "chrom", "hg38", "start", "end", "pctRef","pctDelta",
             "tfbs","motif", "betaLasso","spearman","rfScore","rfNorm",
             "fimo_pvalue",
             "phast7", "phast100","gh","oc","tss","motif_id","ampad.eqtl.pval",
             "ampad.eqtl.beta","ampad.eqtl.score", "gtex.eqtl.pval","gtex.eqtl.beta","gtex.eqtl.score")
    setdiff(coi, colnames(tbl.tbt))
    browser()
    tbl.tbt[, coi]

} # combine.tables
#----------------------------------------------------------------------------------------------------
test_combine.tables <- function()
{
    model <- models[[1]]
    targetGene <- "PILRB"
    TF <- "MZF1"
    tbl.trena <- model$trena
    tbl.trena$tf.rank <- seq_len(nrow(tbl.trena))
    tbl.breaks <- model$tbl.breaks
    tbl.tms <- model$tms.filtered
    tbl.tms$targetGene <- targetGene

    tbl.out <- combine.tables(tbl.trena, tbl.tms, tbl.breaks, targetGene="PILRB", TF="MZF1")
    checkEquals(dim(tbl.out), c(4, 28))

} # test_combine.tables
#----------------------------------------------------------------------------------------------------

#------------------------------------------------------------
# "PILRB-GTEx_V8.Brain_Cerebellar_Hemisphere"
#------------------------------------------------------------
model <- models[[1]]
names(models)[1]

targetGene <- "PILRB"
TF <- "MZF1"
tbl.trena <- model$trena
tbl.breaks <- model$tbl.breaks
tbl.tms <- model$tms.filtered
tbl.tms$targetGene <- targetGene

tbl.trena.tf <- subset(tbl.trena, gene==TF)



tbl.tms.strong <- subset(tbl.tms, tf==TF & ampad.eqtl.score < -1)
tbl.breaks.strong <- subset(tbl.breaks, geneSymbol==TF & pctDelta < -0.1)
gr.breaks.strong <- GRanges(tbl.breaks.strong)
gr.tms.strong <- GRanges(tbl.tms.strong)
tbl.ov <- as.data.frame(findOverlaps(gr.breaks.strong, gr.tms.strong))

tbl.breaks.strong.ov <- tbl.breaks.strong[tbl.ov[,1]][, c(1,3,4)]
colnames(tbl.breaks.strong.ov)[2:3] <- c("hg38", "rsid")


tbl.tms.strong.ov <- tbl.tms.strong[tbl.ov[,2],]

tbl.breaks.tms.strong <- cbind(tbl.tms.strong.ov, tbl.breaks.strong.ov)[, c("SNP_id", "

tbl.trena.tms <- merge(tbl.trena.tf, tbl.tms.strong, by.x="gene", by.y="tf")
tbl.trena.tms.breaks <- merge(tbl.trena.tms, tbl.breaks.strong, by.x="gene", by.y="geneSymbol")
coi <- c("targetGene", "gene",  "SNP_id", "chrom.y", "start.x", "end.x", "start.y", "betaLasso", "spearmanCoeff", "rfScore", "pctRef", "pctDelta",  "gtex.eqtl.pval", "ampad.eqtl.pval", "ampad.eqtl.beta", "ampad.eqtl.score", "tss")
setdiff(coi, colnames(tbl.trena.tms.breaks))
tbl.trena.tms.breaks[, coi]



dim(tbl.tms)
tbl.tmp <- merge(subset(tbl.trena, gene==TF),
                 subset(tbl.tms, tf==TF), by.x="gene", by.y="tf")
tbl.mzf1 <- merge(tbl.tmp, subset(tbl.breaks, geneSymbol==goi), by.x="gene", by.y="geneSymbol")
dim(tbl.mzf1)
coi <- c("targetGene", "gene",  "chrom.y", "start.y", "betaLasso", "spearmanCoeff", "rfScore", "pctRef", "pctDelta", "SNP_id", "gtex.eqtl.pval", "ampad.eqtl.pval", "ampad.eqtl.beta", "ampad.eqtl.score", "tss")
setdiff(coi, colnames(tbl.mzf1))
head(tbl.mzf1[, coi])
tbl.mzf1$ampad.eqtl.score <- with(tbl.mzf1, -log10(ampad.eqtl.pval) * ampad.eqtl.beta)


tbl.trenaBreaks.mzf1 <- merge(subset(tbl.trena, gene=="MZF1"), tbl.breaks, by.x="gene", by.y="geneSymbol")
tbl.trenaBreaks.tead1 <- merge(subset(tbl.trena, gene=="TEAD1"), tbl.breaks, by.x="gene", by.y="geneSymbol")
