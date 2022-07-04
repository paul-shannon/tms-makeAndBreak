library(RPostgreSQL)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
tag.snp <- "rs143080277"
chrom <- "chr2"
tag.snp.hg38 <- 105749599
shoulder <- 1000000
start <- tag.snp.hg38 - shoulder
end   <- tag.snp.hg38 + shoulder
dbGetQuery(geneRegDB, "select * from eqtl2 limit 3")
query <- sprintf("select * from eqtl2 where chrom='%s' and hg38 >= %d and hg38 <= %d", chrom, start, end)

tbl <- dbGetQuery(geneRegDB, query)
dim(tbl) # 45599 12

     # rename C2orf40 to ECRG4 - 

c2orf40.rows <- grep("C2orf40", tbl$genesymbol)
printf("--- rename %d C2orf40 rows to ECRG4", length(c2orf40.rows))
if(length(c2orf40.rows) > 0){
   tbl$genesymbol[c2orf40.rows] <- "ECRG4"
   }

     # fill in empty genesymbol for  ENSG00000272994: C2orf49-DT

ENSG00000272994.rows <- grep("ENSG00000272994", tbl$ensg)
length(ENSG00000272994.rows)  # 3810
tbl$genesymbol[ENSG00000272994.rows] <- "C2orf49-DT"

dim(tbl)
tbl.eqtl.rosmap <- subset(tbl, pvalue < 1e-3)
dim(tbl.eqtl.rosmap)  # 974 at 1e-3 
unique(tbl.eqtl.rosmap$genesymbol)
unique(subset(tbl.eqtl.rosmap, nchar(genesymbol) == 0)$ensg)   # just "ENSG00000269707"
ENSG00000269707.rows <- grep("ENSG00000269707", tbl.eqtl.rosmap$ensg)
length(ENSG00000269707.rows)  # 153
tbl.eqtl.rosmap$genesymbol[ENSG00000269707.rows] <- "Lnc-POU3F3-7"
unique(tbl.eqtl.rosmap$genesymbol)

#"C2orf40" %in% tbl.eqtl.rosmap$genesymbol

as.data.frame(sort(table(tbl.eqtl.rosmap$genesymbol), decreasing=TRUE))
tbl.eqtl.rosmap$score <- -log10(tbl.eqtl.rosmap$pvalue) * tbl.eqtl.rosmap$beta
fivenum(tbl.eqtl.rosmap$score)   #  -71.0453097  -6.7500935  -1.4324213   0.9271386  72.2599049
subset(tbl.eqtl.rosmap, abs(score) > 70)
"C2orf40" %in% tbl.eqtl.rosmap$genesymbol
"C2orf49-DT" %in% tbl.eqtl.rosmap$genesymbol
"Lnc-POU3F3-7"  %in% tbl.eqtl.rosmap$genesymbol
"" %in% tbl.eqtl.rosmap$genesymbol
"ECRG4" %in% tbl.eqtl.rosmap$genesymbol
save(tbl.eqtl.rosmap, file="../shared/tbl.eqtls.rosmap.RData")


