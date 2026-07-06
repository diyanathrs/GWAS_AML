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
chrom.res <- chrom.res[N >=6, ]

## read info score
info.data <- fread("zcat /home/nwl15/WORKING_DATA/HRCimpvData/CLL/HRCinfo_CLL_finalQced4impassoc_HRC.gz")[Rsq >=0.3, ]
##
info.data <- info.data[paste0(chromosome, "_", position) %chin% chrom.res[, paste0(CHR, "_", BP), ]]
## appending info score and overall MAF
chrom.res[info.data, info:=i.Rsq, on=c(SNP="#RSID")][info.data, OverMAF:=i.MAF, on=c(SNP="#RSID")]
infokkk <- info.data[paste0(chromosome, "_", position) %chin% chrom.res[is.na(info), paste0(CHR, "_", BP)], .(chromosome, position, Rsq, MAF)]
chrom.res[infokkk, info:=i.Rsq, on=c(CHR="chromosome", BP="position")][infokkk, OverMAF:=i.MAF, on=c(CHR="chromosome", BP="position")]
gc()

## qq plots
source("/home/nwl15/CamStuff/codesbackup/genetic_resources/Scripts/r/qq.plot.r")
#source("/home/nwl15/CamStuff/codesbackup/genetic_resources/Scripts/r/manhattan.plot.r")

setwd("./MH_QQplots")

## png(paste0(phenoprefix, "_meta_", dataprefix, "mh_HRC.png"), type="cairo", width=14, height=9,units="in", res=300)
## #highKARN1 <- kkkk[CHR==9 & BP>=731563 & BP<=761563, SNP]
## manhattan.plot(chrom.res, chr="CHR", bp="BP", p="P", snp="SNP", lty=2, lcol="firebrick",
##                main=paste0(phenoprefix, " N=", chrom.res[!is.na(P), .N], "/", chrom.res[, .N])
##                )
## dev.off()
infothre <- seq(0.3, 0.8, by=0.1)

png(paste0(phenoprefix, "_meta_" , dataprefix, "qq_HRCbyinfo.png"), type="cairo", width=8, height=6, units="in", res=150)
#qq.plot(chrom.res[, P])
par(mfrow=c(2,3), oma=c(2,2,3,1))
for (n in seq_along(infothre)){
    qq.plot(chrom.res[ info >= infothre[n], P])
    b <- chrom.res[!is.na(P) & info >= infothre[n], .N]
    mtext(paste0("Info score of >= ", infothre[n], "\n N=", b), side=3, cex=0.8)
}
mtext(paste0("studyMAF > 1% & present in 6 studies"), side=3, outer=TRUE, cex=1)

dev.off()

mafthre <- seq(0.01, 0.06, by=0.01)

png(paste0(phenoprefix, "_meta_" , dataprefix, "qq_HRCbyMAF.png"), type="cairo", width=8, height=6, units="in", res=150)
#qq.plot(chrom.res[, P])
par(mfrow=c(2,3), oma=c(2,2,3,1))
for (n in seq_along(mafthre)){
    qq.plot(chrom.res[OverMAF >= mafthre[n], P])
    b <- chrom.res[!is.na(P) & OverMAF >= mafthre[n], .N]
    mtext(paste0("overall MAF >= ", mafthre[n], "\n N=", b), side=3, cex=0.8)
}
mtext(paste0("studyMAF > 1% & present in 6 studies, info >=0.3"), side=3, outer=TRUE, cex=1)

dev.off()



