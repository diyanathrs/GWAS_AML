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
phenoprefix
dataprefix
res.dir <- "../"
#chrom

## metaRes

metaRes <- fread(paste0(res.dir, phenoprefix, "_NCL_", dataprefix, ".meta"), header = TRUE)
metaRes[, rsid := SNP]


## HRC rsid
hrcsnpid <- fread("cut -f1-5 /nobackup/proj/jamgaml/dean_AML/HRCstuff/HRC.r1-1.GRCh37.wgs.mac5.sites.tab",
                  header = TRUE
                  )[ID != ".", ]
##hrcsnpid[POS==693731,]
hrcsnpid[, SNP := paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]

metaRes[hrcsnpid, rsid := i.ID, on=c(SNP="SNP")]

setorder(metaRes, CHR, BP)


options(scipen =999999)
fwrite(x=metaRes,
       sep="\t",
       file=paste0(phenoprefix, "_", dataprefix, "_assocLZ.res"), quote=FALSE,
       row.names=FALSE, col.names=TRUE)
