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

#res.dir <- "/home/nwl15/WORKING_DATA/Newcastle3/qcedData/ResultSummary"
#phenoprefix <- "status"
#dataprefix <- "NCL1_2_finalQCed"

res.dir
phenoprefix
dataprefix

#snptest.cols <- c("rsid", "chromosome", "position", "alleleA", "alleleB", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1")
snptest.cols <- c("rsid", "chromosome", "position", "alleleA", "alleleB", "minorAllele", "crude.pvalue", "beta", "SE")

filein <- paste0(file.path(res.dir, phenoprefix), "_", dataprefix, "assoc_HRC.gz")
if (!file.exists(filein)){
    stop(paste0(filein, " doesn't exist"))
}


chrom.res <- fread(paste0("zcat ", filein), header = TRUE)
setnames(chrom.res, c("#rsid"), c("rsid"))

chrom.res[, rsid1 := paste0(chromosome,":", position,"_", alleleA, "_", alleleB)][
    grepl("^esv", rsid), rsid1:=rsid]

## chroms <- seq(1, 22)
## chrom.res <- vector("list", length=22)
## for (n in seq_along(chroms)){
##     filein <- paste0(res.dir, "/", phenoprefix, "_", dataprefix, "_chrom" , chroms[n], "_cleaned.rds")
##     ## plink can't cope with extremely long string of SNP names
##     ## nchar(rsid1)>150, all of them are CNV
##     chrom.res[[n]] <- readRDS(filein)[, snptest.cols, with=FALSE][,
##             rsid1:=paste0(chromosome,":",position,"_", alleleA, "_", alleleB)][grepl("^esv", rsid), rsid1:=rsid]
## }
## chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)

chrom.res[alleleB == minorAllele, ':=' (A2=alleleA, A1=alleleB)][
    alleleA == minorAllele, ':=' (A2=alleleB, A1=alleleA) ]


#chrom.res[, ':=' (A2=alleleA, A1=alleleB)]



##chrom.res[alleleA == minorAllele, beta := -1*beta]

## rsid1 should be unique 
exclude.out <- chrom.res[, .N, by=rsid1][N>1, rsid1]
if (length(exclude.out)!=0){
    write.table(chrom.res[rsid1 %chin% exclude.out, ],
                file=paste0(phenoprefix, "_", dataprefix,"_excluded.out"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE )}
## 
## setnames(chrom.res,
##          c("rsid1", "chromosome", "position", "alleleA", "alleleB", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1"),
##          c("SNP", "CHR", "BP", "A2", "A1", "P", "BETA", "SE")
##          )

setnames(chrom.res,
         c("rsid1", "chromosome", "position", "crude.pvalue", "beta"),
         c("SNP", "CHR", "BP", "P", "BETA")
         )


chrom.res[studyMAF < 0.01 | studyMAF > 0.99, ':=' (BETA=NA, P=NA, SE=NA)]
chrom.res[info_score < 0.6, ':=' (BETA=NA, P=NA, SE=NA)]

## disable scitific notation
options(scipen =99999)

fwrite(chrom.res[!(SNP %chin% exclude.out),
                 c("SNP", "CHR", "BP", "A2", "A1", "BETA", "SE", "P"), with=FALSE],
       sep = "\t",
       na = "NA",
       file=paste0(phenoprefix, "_", dataprefix,"_assoc.res"), quote=FALSE,
       row.names=FALSE, col.names=TRUE)


