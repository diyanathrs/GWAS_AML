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
## phenoprefix <- "status"
## dataprefix <- "NCL1_2_finalQCed"
## chrom <- 22
## res.dir <- "/home/nwl15/WORKING_DATA/Newcastle4/QCsteps_gtconsole_realPheno/qcedData/associationTest"
#phenoprefix
#dataprefix

#res.dir <- paste0(root.dir, "/", "association")
root.dir

info.dir <- paste0(root.dir, "/", "impvRes")
#root.dir <- "/home/nwl15/WORKING_DATA/HRCimpvData/CLL"
#chrom <- 2

chrom
#dataprefix <- "CLL_finalQced4imp"

dataprefix

chrom.log.id <- paste0(dataprefix, "_chr", chrom, "_HRCMarkers.log")
chrom.log     <- file(chrom.log.id, "w")

## info and maf from michigan imputation server
snpStats.file <- paste0("zcat ", info.dir , "/", "chr", chrom, ".info.gz")
snpStats <- fread(snpStats.file, header=TRUE, colClasses="character", na.strings=c("-"))
char2num1 <- c("ALT_Frq", "MAF", "AvgCall", "Rsq", "LooRsq", "EmpR", "EmpRsq", "Dose0", "Dose1")
snpStats[, (char2num1):=lapply(.SD, as.numeric), .SDcols=char2num1]
setnames(snpStats, c("REF(0)", "ALT(1)"), c("REF", "ALT"))
cat("\nNo. variants on chrom", chrom, " are ", snpStats[, .N], " in the imputation output \n", file=chrom.log, sep="")
cat("  --No. Genotyped variants are ", snpStats[Genotyped=="Genotyped", .N], "\n", file=chrom.log, sep="")
cat("  --No. Imputed variants are ", snpStats[Genotyped=="Imputed", .N], "\n", file=chrom.log, sep="")
cat("  --No. Typed_Only variants are ", snpStats[Genotyped=="Typed_Only", .N], "\n", file=chrom.log, sep="")
snpStats <- snpStats[(MAF >=0.0095 & Rsq >=0.3) | (Genotyped=="Typed_Only" & MAF >= 0.0095), ]
#snpStats[, RSID:= paste0(SNP, "_", REF, "_", ALT)]
cat("No. variants (MAF of >=1% and INFO >=0.3) are ", snpStats[, .N], "\n\n", file=chrom.log, sep="")

close(chrom.log)

cat(snpStats[, SNP], file=paste0(dataprefix, "_chr", chrom, "_MarkersPassed.lst"), sep="\n")

