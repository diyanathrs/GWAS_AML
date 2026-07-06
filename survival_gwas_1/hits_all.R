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
## res.dir <- "/home/nwl15/WORKING_DATA/Newcastle4/QCsteps_CRLMMdata_realPheno/qcedData/ResultSummary"
## res.dir1 <- "/home/nwl15/WORKING_DATA/Newcastle4/QCsteps_CRLMMdata_realPheno/qcedData/imputedMergeddata"
## dataprefix <- "NCL4_crlmmQCed"
## phenoprefix <- "status"
#res.dir
#res.dir1

phenoprefix
dataprefix
root.dir

res.dir <- paste0(root.dir, "/", "ResultSummary")
res.dir
info.dir <- paste0(root.dir, "/", "impvRes")
info.dir

## snps of interest
snps <- fread("cat *_hits.lst", header=FALSE, sep="\t")[[1]]

chroms <- seq(1, 22)
chrom.res <- vector("list", length=22)
for (n in seq_along(chroms)){
    snpStats.file <- paste0("zcat ", info.dir , "/", "chr", chroms[n], ".info.gz")
    snpStats <- fread(snpStats.file, header=TRUE, colClasses="character", na.strings=c("-"))[Genotyped!="Imputed", ]
    char2num1 <- c("ALT_Frq", "MAF", "AvgCall", "Rsq", "LooRsq", "EmpR", "EmpRsq", "Dose0", "Dose1")
    snpStats[, (char2num1):=lapply(.SD, as.numeric), .SDcols=char2num1]
    setnames(snpStats, c("REF(0)", "ALT(1)"), c("REF", "ALT"))
    #infodata <- readRDS(paste0(res.dir1, "/" , dataprefix, "_chr", chroms[n], "_variantsInfo.rds"))[type==2, position]
    filein <- paste0(res.dir, "/", phenoprefix, "_", dataprefix, "_chrom" , chroms[n], "_cleaned.rds")
    chrom.res[[n]] <- readRDS(filein)[rsid %chin% snps, ][, typed:="imputed"][
        paste0(chromosome, ":", position, "_", alleleA, "_", alleleB) %chin% snpStats[, paste0(SNP, "_", REF, "_", ALT)], typed:="typed"]
}

chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)[rsid %chin% snps, ]

if (chrom.res[, .N]!=0 ){
   #cat(chrom.res[, rsid], file=paste0(dataprefix, "_", phenoprefix, "_hits.lst"), sep="\n")
   saveRDS(chrom.res, file=paste0(dataprefix, "_", phenoprefix, "_hitsALL.rds" ))
}


