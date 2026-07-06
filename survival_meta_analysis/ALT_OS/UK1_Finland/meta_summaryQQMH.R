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
#phenoprefix <- "status"
#dataprefix <- "gt"

## read meta res
chrom.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, ".meta"))

## logs of no. variants
chrom.log.id <- paste0(phenoprefix, "_", dataprefix, "_metaNovariants", ".log")
chrom.log     <- file(chrom.log.id, "w")


cat("No. variants in 6 data set are  ", chrom.res[N==6, .N] , "\n", file=chrom.log, sep="")
cat("No. variants in 5 data set are  ", chrom.res[N==5, .N] , "\n", file=chrom.log, sep="")
cat("No. variants in 4 data set are  ", chrom.res[N==4, .N] , "\n", file=chrom.log, sep="")
cat("No. variants in 3 data set are  ", chrom.res[N==3, .N] , "\n", file=chrom.log, sep="")
cat("No. variants in 2 data set are  ", chrom.res[N==2, .N] , "\n", file=chrom.log, sep="")
cat("No. variants in 1 data set are  ", chrom.res[N==1, .N] , "\n", file=chrom.log, sep="")
cat("No. variants with NA p values are ", chrom.res[is.na(P), .N], "\n", file=chrom.log, sep="")
cat("No. variants with NA p values and present in only 1 data set are ", chrom.res[is.na(P) & N==1, .N], "\n", file=chrom.log, sep="")
close(chrom.log)


## qq plots
source("/nobackup/nwl15/CamStuff/codesbackup/genetic_resources/Scripts/r/qq.plot.r")
source("/nobackup/nwl15/CamStuff/codesbackup/genetic_resources/Scripts/r/manhattan.plot.r")

setwd("./MH_QQplots")

png(paste0(phenoprefix, "_meta_", dataprefix, "mh_HRC.png"), type="cairo", width=14, height=9,units="in", res=300)
#highKARN1 <- kkkk[CHR==9 & BP>=731563 & BP<=761563, SNP]
manhattan.plot(chrom.res, chr="CHR", bp="BP", p="P", snp="SNP", lty=2, lcol="firebrick",
               main=paste0(phenoprefix, " meta Results; No variants=", chrom.res[!is.na(P), .N], "/", chrom.res[, .N])
               )
dev.off()

png(paste0(phenoprefix, "_meta_" , dataprefix, "qq_HRC.png"), type="cairo", width=8, height=6, units="in", res=300)
par(oma=c(2,2,3,1))
qq.plot(chrom.res[, P])
b <- chrom.res[!is.na(P), .N]
mtext(paste0("No variants =", b), side=3, cex=1)

dev.off()
