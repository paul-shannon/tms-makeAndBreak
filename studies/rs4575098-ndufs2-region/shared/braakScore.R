library(EndophenotypeExplorer)
library(RUnit)
targetGene <- "NDUFS2"
if(!exists("etx")){
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD",
                                    initialize.snpLocs=TRUE)
   tbl.eqtls <- etx$get.ampad.EQTLsForGene()
   dim(tbl.eqtls)
   head(tbl.eqtls)
   table(tbl.eqtls$study)   #   ampad-mayo ampad-rosmap         GTEx
                            #        7532         3774            4
   }
#----------------------------------------------------------------------------------------------------
rosmap.ad.ctl.separation <- function(rsid)
{
   printf("----- rosmap.ad.ctl.separation: %s", rsid)

   result <- list(pval.t=NA,
                  pval.fisher=NA,
                  tbl.geno=data.frame(),
                  tbl.pt=data.frame(),
                  mtx.geno=matrix(),
                  mtx.geno.study=matrix(),
                  pt.ad=c(),
                  pt.ctl=c())

   mtx.geno.1 <- etx$getGenoMatrixByRSID(rsid)
   if(all(is.na(mtx.geno.1)))
       return(result)

   mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")
   if(all(is.na(mtx.geno.pt.rosmap)))
       return(result)

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.rosmap))
   length(rosmap.patients)   # 1143

   tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
   dim(tbl.pt.rosmap)  # 1143 18

   dim(mtx.geno.pt.rosmap)
   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              # 0/0 0/1 1/1
                              # 827 296  28
   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85


   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              #  0/0 0/1 1/1
                              #  827 296  28

   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85

   ad.wt <- 0
   at.het <- 0
   ad.hom <- 0

   ad.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ad])))

   ctl.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ctl])))

   tbl.summary <- data.frame(wt=c(ad.wt, ctl.wt), het=c(ad.het, ctl.het), hom=c(ad.hom, ctl.hom),
                             row.names=c("ad", "ctl"))
   print(tbl.summary)

   ad.vector <- with(tbl.summary["ad",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))
   ctl.vector <- with(tbl.summary["ctl",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))

   #pval.t <- t.test(ad.vector, ctl.vector)$p.value
   #pval.fisher <- fisher.test(tbl.summary)$p.value

   pval.t <- tryCatch({
       t.test(ad.vector, ctl.vector)$p.value
       }, error=function(e){return(1)})

   pval.fisher <- tryCatch({
      fisher.test(tbl.summary)$p.value
      }, error=function(e){return(1)})



   return(list(pval.t=pval.t,
               pval.fisher=pval.fisher,
               tbl.geno=tbl.summary,
               tbl.pt=tbl.pt.rosmap,
               mtx.geno=mtx.geno.1,
               mtx.geno.study=mtx.geno.pt.rosmap,
               pt.ad=pt.ad,
               pt.ctl=pt.ctl))


} # rosmap.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
scores <- list()
rsids <- get(load("rsids.RData"))
deleters <- which(nchar(rsids) == 0)
if(length(deleters) > 0)
    rsids <- rsids[-deleters]
deleters <- which(is.na(rsids))
if(length(deleters) > 0)
    rsids <- rsids[-deleters]

length(rsids)
i <- 0
for(rsid in rsids[2447:2471]){
   i <- i + 1
   printf("--- %d) %s", i, rsid)
   scores[[rsid]] <- rosmap.ad.ctl.separation(rsid)
   }


f <-  sprintf("rsid.braak.scores-%s.RData", gsub(" ", ".", Sys.time(), fixed=TRUE))
save(scores, file=f)


