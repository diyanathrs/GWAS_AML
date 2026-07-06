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
inf <- 0.6
maf <- 0.0095
chrom.res <- rbindlist(lapply(chrom.res, function(x) x[info>=inf & all_maf>= maf, ]), use.names=TRUE, fill=TRUE)
#chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)

## qq plots
source("/mnt/storage/nobackup/proj/jamgaml/dean_AML/codesbackup/genetic_resources/Scripts/r/qq.plot.r")
source("/mnt/storage/nobackup/proj/jamgaml/dean_AML/codesbackup/genetic_resources/Scripts/r/manhattan.plot.r")



if (!dir.exists("MH_QQplots_info4_maf01")){
    dir.create(file.path(getwd(), "MH_QQplots_info4_maf01"))
    setwd("./MH_QQplots_info4_maf01")
} else {
    setwd("./MH_QQplots_info4_maf01/")}

## MH plot
png(paste0(phenoprefix, dataprefix, "mhHRCres_info",inf,"maf",maf,".png"), type="cairo", width=14, height=9, units="in", res=300)
manhattan.plot(chrom.res, chr="chromosome", bp="position", p="frequentist_add_pvalue", snp="rsid", lty=2, lcol="firebrick",
               main=paste0(phenoprefix, " N=", chrom.res[ , .N])
               )
dev.off()

## QQ plot
png(paste0(phenoprefix, dataprefix, "qqHRCres_info",inf,"maf",maf,".png"), type="cairo", width=8, height=6, units="in", res=300)
qq.plot(chrom.res[, frequentist_add_pvalue])
dev.off()

