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
#root.dir
#res.dir <- paste0(root.dir, "/", "ResultSummary")
#res.dir

#info.dir <- paste0(root.dir, "/", "impvRes")

## info regarding info, typed, and rsid
## infoscore and typed are no longer valid as more studies are pulled together
## get rsid from HRC variant list 

## cleanedrds <- paste0(res.dir, "/", "Hull1_", phenoprefix, "_", dataprefix, "_chrom", seq(1, 22), "_cleaned.rds")
## allrds <- lapply(cleanedrds, readRDS)
## allrds <- rbindlist(lapply(allrds, function(x) x[alleleA == minorAllele, metaid:=paste0(snp, "_", alleleB, "_", alleleA)][alleleB == minorAllele, metaid:=paste0(snp, "_", alleleA, "_", alleleB)]), use.names=T, fill=T)[, .(rsid, chromosome, position, alleleA, alleleB, minorAllele, metaid, info_score, Genotyped)]

hrcfilein <- paste0("cut -f1-8 /nobackup/proj/jamgaml/HRCstuff/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)[ID != ".", ]
hrc.rsid[, ':=' (metaid=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT),
                 metaid1=paste0(`#CHROM`, ":", POS, "_", ALT, "_", REF)
                 )
         ]

## hrc.rsid[AF <=0.5, metaid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)][
##    AF > 0.5, metaid:=paste0(`#CHROM`, ":", POS, "_", ALT, "_", REF)] 
## hrc.rsid <- hrc.rsid[ID != ".", ]

## Plink meta results
plinkMeta <- fread(paste0("../", phenoprefix, "_NCL_", dataprefix, ".meta"), header=T)
plinkMeta[, metaid:=paste0(CHR, ":", BP, "_", A2, "_", A1)]

## stopifnot(all(plinkMeta$metaid %chin% allrds$metaid))
## stopifnot(all(plinkMeta$metaid %chin% hrc.rsid$metaid))

## update with rsid, info_score, typed
plinkMeta[hrc.rsid, rsid:=i.ID, on=c(metaid="metaid")][hrc.rsid, rsid:=i.ID, on=c(metaid="metaid1")]
plinkMeta[is.na(rsid), rsid := metaid]

#[allrds, infoScore:=i.info_score, on=c(metaid="metaid")][allrds,
                                                      #                                                         Genotyped:=i.Genotyped, on=c(metaid="metaid")]

setnames(plinkMeta, c("P", "P(R)", "OR", "OR(R)"), c("fixedP", "randomP", "fixedHR", "randomHR"))
setorder(plinkMeta, CHR, BP)

## disable scitific notation
options(scipen =99999)
write.table(plinkMeta, sep="\t",
            file=paste0(phenoprefix, "_", dataprefix, "_assocLZ.res"), quote=FALSE, row.names=FALSE, col.names=TRUE)


