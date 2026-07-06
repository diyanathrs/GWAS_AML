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
## phenoprefix <- "Hull1_OS_Dx_Death_LFU_Status"
## dataprefix <- "CLL_finalQced4imp"
## chrom <- 22
## res.dir <- "/home/nwl15/WORKING_DATA/Newcastle4/QCsteps_gtconsole_realPheno/qcedData/associationTest"
phenoprefix
dataprefix
#chrom
root.dir
res.dir <- paste0(root.dir, "/", "ResultSummary")
res.dir

#info.dir <- paste0(root.dir, "/", "impvRes")

## info regarding info, typed, and rsid
cleanedrds <- paste0(res.dir, "/", "Hull1_", phenoprefix, "_", dataprefix, "_chrom", seq(1, 22), "_cleaned.rds")
allrds <- lapply(cleanedrds, readRDS)
allrds <- rbindlist(lapply(allrds, function(x) x[alleleA == minorAllele, metaid:=paste0(snp, "_", alleleB, "_", alleleA)][alleleB == minorAllele, metaid:=paste0(snp, "_", alleleA, "_", alleleB)]), use.names=T, fill=T)[, .(rsid, chromosome, position, alleleA, alleleB, minorAllele, metaid, info_score, Genotyped)]

## Plink meta results
plinkMeta <- fread(paste0("../", phenoprefix, "_NCL_", dataprefix, ".meta"), header=T)
plinkMeta[, metaid:=paste0(CHR, ":", BP, "_", A2, "_", A1)]

stopifnot(all(plinkMeta$metaid %chin% allrds$metaid))

## update with rsid, info_score, typed
plinkMeta[allrds, rsid:=i.rsid, on=c(metaid="metaid")][allrds, infoScore:=i.info_score, on=c(metaid="metaid")][allrds,
                                                                                                               Genotyped:=i.Genotyped, on=c(metaid="metaid")]

setnames(plinkMeta, c("P", "P(R)", "OR", "OR(R)"), c("fixedP", "randomP", "fixedHR", "randomHR"))
setorder(plinkMeta, CHR, BP)

## disable scitific notation
options(scipen =999)
write.table(plinkMeta, sep="\t",
            file=paste0(phenoprefix, "_", dataprefix, "_assocLZ.res"), quote=FALSE, row.names=FALSE, col.names=TRUE)


