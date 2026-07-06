args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

require("data.table")

## res.dir <- "/home/nwl15/WORKING_DATA/HRCimpvData/NCL1_2/ResultSummary"
## phenoprefix <- "status"
## dataprefix <- "NCL1_2_finalQCed"

res.dir
phenoprefix
dataprefix

#snptest.cols <- c("rsid", "chromosome", "position", "alleleA", "alleleB", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1")

chroms <- seq(1, 22)
chrom.res <- vector("list", length=22)
for (n in seq_along(chroms)){
    filein <- paste0(res.dir, "/", phenoprefix, "_", dataprefix, "_chrom" , chroms[n], "_cleaned.rds")
    chrom.res[[n]] <- readRDS(filein)
    #[
    #    !grepl("^rs|^esv", rsid), rsid:=paste0("chr", chromosome, ":", position)]

}
chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)

## setnames(chrom.res,
##          c("rsid", "chromosome", "position", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1"),
##          c("SNP", "CHR", "BP", "P", "BETA", "SE")
##          )
## chrom.res[, (c("chr", "snp", "clogP", "alogP")):=NULL]

chrom.res[, (c("chr", "snp", "clogP")):=NULL]
abc <- c("rsid", "chromosome", "position")
efg <- setdiff(names(chrom.res), abc)
#col.orders <-  c("SNP", "CHR", "BP", "alleleA", "alleleB", "info", "all_AA", "all_AB",
#                 "all_BB", "all_maf", "cases_maf", "controls_maf", "BETA", "SE", "P")
col.orders <- c(abc, efg)
col.orders


## disable scitific notation
options(scipen =99999)
write.table(chrom.res[, col.orders, with=FALSE],
            sep="\t",
            file=paste0(phenoprefix, "_", dataprefix, "_assocLZ.res"), quote=FALSE,
            row.names=FALSE, col.names=TRUE)



