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
#dataprefix <- "NCL4_gtQCed"
## chrom <- 1
## res.dir <- "/home/nwl15/WORKING_DATA/Newcastle4/QCsteps_gtconsole_realPheno/qcedData/associationTest"
phenoprefix
dataprefix
#chrom
#res.dir

chroms <- seq(1, 22)
chrom.res <- vector("list", length=22)
for (n in seq_along(chroms)){
    filein <- paste0(phenoprefix, "_", dataprefix, "_chrom" , chroms[n], "_cleaned.rds")
    #print(filein)
    chrom.res[[n]] <- readRDS(filein)
}

#chrom.res <- rbindlist(lapply(chrom.res, function(x) x[info_score>=0.6 & minorAllele>=0.02, ]), use.names=TRUE, fill=TRUE)
chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)

## qq plots
source("/nobackup/proj/jamgaml/dean_AML/codesbackup/genetic_resources/Scripts/r/qq.plot.r")
source("/nobackup/proj/jamgaml/dean_AML/codesbackup/genetic_resources/Scripts/r/manhattan.plot.r")


if (!dir.exists("MH_QQplots")){
    dir.create(file.path(getwd(), "MH_QQplots"))
    setwd("./MH_QQplots")
} else {
    setwd("./MH_QQplots/")}


## crude

## MH
png(paste0(phenoprefix, dataprefix, "mhcrude95maf3infoHRC.png"), type="cairo", width=14, height=9, units="in", res=300)

manhattan.plot(chrom.res, chr="chromosome", bp="position", p="crude.pvalue",
               snp="rsid", lty=2, lcol="firebrick",
               ##main=paste0(phenoprefix, " (HRC crude, MAFs of >=2% and info of >=0.6); \nN =", chrom.res[!is.na(crude.pvalue), .N], "/", chrom.res[ , .N])
               main=paste0(phenoprefix, " (HRC crude, MAFs of >=0.95% and info of >=0.3); \nN =", chrom.res[!is.na(crude.pvalue), .N], "/", chrom.res[ , .N])
               )
dev.off()

## QQ
png(paste0(phenoprefix, dataprefix, "qqcrude95maf3infoHRC.png"), type="cairo", width=8, height=6, units="in", res=300)

par(oma=c(2,2,3,1))
qq.plot(chrom.res[, crude.pvalue])
b <- chrom.res[!is.na(crude.pvalue), .N]
#mtext(paste0("No variants (MAFs of >= 2% and info of >= 0.6)=", b), side=3, cex=1)
mtext(paste0("No variants ( MAFs of >= 0.95% and info of >=0.3)=", b), side=3, cex=1)

dev.off()

## adjusted

## png(paste0(phenoprefix, dataprefix, "mhadjustedHRC.png"), type="cairo", width=14, height=9,units="in", res=300)

## manhattan.plot(chrom.res, chr="chromosome", bp="position", p="a.pvalue",
##                snp="rsid", lty=2, lcol="firebrick",
##                main=paste0(phenoprefix, " (HRC adjusted); N=", chrom.res[ , .N])
##                )
## dev.off()

## png(paste0(phenoprefix, dataprefix, "qqadjustedHRC.png"), type="cairo", width=8, height=6, units="in", res=300)

## qq.plot(chrom.res[, a.pvalue])

## dev.off()














## MH plot
## png(paste0(phenoprefix, dataprefix, "mhHRCres.png"), type="cairo", width=14, height=9, units="in", res=300)
## manhattan.plot(chrom.res, chr="chromosome", bp="position", p="frequentist_add_pvalue", snp="rsid", lty=2, lcol="firebrick",
##                main=paste0(phenoprefix, " N=", chrom.res[ , .N])
##                )
## dev.off()

## QQ plot
## png(paste0(phenoprefix, dataprefix, "qqHRCres.png"), type="cairo", width=8, height=6, units="in", res=300)
## qq.plot(chrom.res[, frequentist_add_pvalue])
## dev.off()

