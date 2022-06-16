library(RPostgreSQL)
geneRegDB <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="genereg2021", host="khaleesi")
tag.snp.loc <- 100334426
shoulder <- 1000000
chrom <- "chr7"
start <- tag.snp.loc - shoulder
end   <- tag.snp.loc + shoulder
dbGetQuery(geneRegDB, "select * from eqtl2 limit 3")
dbGetQuery(geneRegDB, "select * from eqtl2 where chrom='chr7'and hg38  limit 3")
query <- sprintf("select * from eqtl2 where chrom='%s' and hg38 >= %d and hg38 <= %d", chrom, start, end)

tbl <- dbGetQuery(geneRegDB, query)
dim(tbl) # 129426 12
head(as.data.frame(sort(table(tbl$genesymbol), decreasing=TRUE)))
save(tbl, file="tbl.eqtls.rosmap.RData")

