library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
    #--------------------------
    # hg38
    #--------------------------

tag.snp <- "rs4575098"
chrom <- "chr1"
min.loc.hg38 <- 160684126
max.loc.hg38 <- 161688825
printf("%dk span", round((max.loc.hg38 - min.loc.hg38)/1000))   # 1005k

slice <- GRanges(seqnames="1", IRanges(min.loc.hg38, max.loc.hg38))
snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, slice, genome="GRCh38")
class(snps)
length(snps)
head(snps)
tbl.raw <- as.data.frame(snps)
tbl.hg38 <- tbl.raw[, c(1,2,4,7,8)]
colnames(tbl.hg38) <- c("chrom", "hg38", "rsid", "ref", "alt")
head(tbl.hg38)

    #------------------------------------------------------
    # try liftover of more complete hg38 variants to hg19
    #------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)
dim(tbl.hg38)
gr.hg38 <- GRanges(seqnames="chr1", IRanges(start=tbl.hg38$hg38, end=tbl.hg38$hg38+1))
length(gr.hg38)
gr.hg38$rsid <- tbl.hg38$rsid

seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")["chr1"]
length(gr.hg38)

chain.gz <- "hg38ToHg19.over.chain.gz"
system(sprintf("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/%s", chain.gz))
system(sprintf("gunzip %s", chain.gz))

chain <- import.chain(sub(".gz", "", chain.gz))
x <- liftOver(gr.hg38, chain)
gr.hg19 <- unlist(x)
head(gr.hg19)
tbl.hg19 <- as.data.frame(gr.hg19)
colnames(tbl.hg19)[grep("start", colnames(tbl.hg19))] <- "hg19"
dim(tbl.hg19)

dim(tbl.hg38)
dim(tbl.hg19)
tbl.rsid <- merge(tbl.hg38, tbl.hg19[, c("rsid", "hg19")], all.x=TRUE)
dim(tbl.rsid)
tbl.rsid <- tbl.rsid[, c("chrom", "hg38", "hg19", "rsid", "ref", "alt")]
dim(tbl.rsid)
if(!grepl("chr", tbl.rsid$chrom[1]))
    tbl.rsid$chrom <- paste0("chr", tbl.rsid$chrom)
head(tbl.rsid)
lapply(tbl.rsid, class)
filename <- sprintf("tbl.rsid.%d-%d.hg38.hg19.RData", min.loc.hg38, max.loc.hg38)
full.path <- file.path("~/github/tms-makeAndBreak/studies/rs4575098-ndufs2-region/shared", filename)
save(tbl.rsid, file=full.path)



