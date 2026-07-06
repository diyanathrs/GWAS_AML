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
dataprefix
#chrom
root.dir
#res.dir <- paste0(root.dir, "/", "association")
info.dir <- paste0(root.dir, "/", "impvRes")

chrom.res <- vector("list", length=22)
for (chrom in seq(1, 22)){
    snpStats.file <- paste0("zcat ", info.dir , "/", "chr", chrom, ".info.gz")
    chrom.res[[chrom]] <- fread(snpStats.file, header=TRUE, colClasses="character", na.strings=c("-"))
}

char2num1 <- c("ALT_Frq", "MAF", "AvgCall", "Rsq", "LooRsq", "EmpR", "EmpRsq", "Dose0", "Dose1")
lapply(chrom.res, function(x) x[, (char2num1):=lapply(.SD, as.numeric), .SDcols=char2num1])
lapply(chrom.res, function(x) x[, ':=' (chromosome=as.integer(sapply(strsplit(SNP, ":"), function(y) y[1])),
                                        position=as.integer(sapply(strsplit(SNP, ":"), function(y) y[2]))
                                        ) ]
       )

lapply(chrom.res, function(x) x[, RSID:= paste0(SNP, "_", `REF(0)`, "_", `ALT(1)`)])

chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)
setnames(chrom.res, c("REF(0)", "ALT(1)"), c("REF", "ALT"))
col.wanted <- c("RSID", "chromosome", "position", names(chrom.res)[c(1:8)])
setorder(chrom.res, chromosome, position)

options(scipen =99999)
write.table(chrom.res[, col.wanted, with=FALSE],
            sep="\t",
            file=paste0(phenoprefix, "_", dataprefix, "_assocLZ.res"), quote=FALSE,
            row.names=FALSE, col.names=TRUE)
