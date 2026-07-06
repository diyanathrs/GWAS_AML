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
chrom
root.dir
assoc.fol 
res.dir <- paste0(root.dir, "/", assoc.fol)
info.dir <- paste0(root.dir, "/", "impvRes")


chrom.log.id <- paste0(phenoprefix, "_", dataprefix, "_chr", chrom, "_HRCsur.log")
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
snpStats[, RSID:= paste0(SNP, "_", REF, "_", ALT)]
cat("No. variants (MAF of >=1% and INFO >=0.3) are ", snpStats[, .N], "\n\n", file=chrom.log, sep="")

######################
## SNPTEST results  ##
######################
chrom.files <- paste0(res.dir, "/", phenoprefix, "_chr", chrom, "_HRCSurRes.out")
#filein <- paste0("zcat ", chrom.files, " | grep -v ^#  | cut -d' ' -f 2-6,9,14-16,29-31,45,47,48")
#filein <- paste0("zcat ", chrom.files, " | grep -v ^# ")
#chrom.res <-  fread(filein, header=TRUE, colClasses="character")
chrom.res <-  fread(chrom.files, header=TRUE, colClasses="character")
cat("No. variants in SNPTEST results are ", chrom.res[, .N], "\n", file=chrom.log, sep="")
#chrom.res1 <- rbindlist(lapply(chrom.res, function(x) x[info>=0.3 & all_maf>=0.01, ]), use.names=TRUE, fill=TRUE)
#chrom.res1 <- rbindlist(lapply(chrom.res, function(x) x[rsid==".",rsid:=paste0(chrom,":", position, "_", alleleA, "_", alleleB)][rsid %chin% snpStats[,RSID], ]), use.names=TRUE, fill=TRUE)
chrom.res1 <- chrom.res[paste0(snp, "_", alleleA, "_", alleleB) %chin% snpStats[, RSID], ][, rsid:=paste0(snp, "_", alleleA, "_", alleleB)]
#[,':=' (a=sapply(strsplit(frequentist_add_pvalue, split="e|E"), function(x) x[1]), b=sapply(strsplit(frequentist_add_pvalue, split="e|E"), function(x) x[2]))][is.na(b), logP:=-log10(as.numeric(a))][!is.na(b), logP:=-(log10(as.numeric(a)) + as.numeric(b))][,(c("a","b")):=NULL]
cat("No. variants with MAFs of >=1% and info >= 0.3 are ", chrom.res1[,.N], "\n\n", file=chrom.log, sep="")
close(chrom.log)

stopifnot(snpStats[,.N]==chrom.res1[,.N])

## crude p values
chrom.res1[,':=' (a=sapply(strsplit(crude.pvalue, split="e|E"), function(x) x[1]), b=sapply(strsplit(crude.pvalue, split="e|E"), function(x) x[2]))][is.na(b), clogP:=-log10(as.numeric(a))][!is.na(b), clogP:=-(log10(as.numeric(a)) + as.numeric(b))][,(c("a","b")):=NULL]
## adjusted p values
#chrom.res1[,':=' (a=sapply(strsplit(a.pvalue, split="e|E"), function(x) x[1]), b=sapply(strsplit(a.pvalue, split="e|E"), function(x) x[2]))][is.na(b), alogP:=-log10(as.numeric(a))][!is.na(b), alogP:=-(log10(as.numeric(a)) + as.numeric(b))][,(c("a","b")):=NULL]


#chrom.res1[,':=' (a=sapply(strsplit(frequentist_add_pvalue, split="e|E"), function(x) x[1]), b=sapply(strsplit(frequentist_add_pvalue, split="e|E"), function(x) x[2]))][is.na(b), logP:=-log10(as.numeric(a))][!is.na(b), logP:=-(log10(as.numeric(a)) + as.numeric(b))][,(c("a","b")):=NULL]
char2num <- setdiff(names(chrom.res1),
                    c("rsid", "chr", "snp", "alleleA", "alleleB", "minorAllele", "clogP"))

chrom.res1[, (char2num):=lapply(.SD, as.numeric), .SDcols=char2num][, chromosome:=as.numeric(chrom)]

## appending info score and type of variants
chrom.res1[snpStats, info_score:=i.Rsq, on=c(rsid="RSID")][snpStats, Genotyped:=i.Genotyped, on=c(rsid="RSID")]

stopifnot(chrom.res1[is.na(info_score), .N] ==0 )

setnames(chrom.res1, c("posb37"), c("position"))
setorder(chrom.res1, chromosome, position)

#####################
## rsid assignment ##
#####################
##hrcfilein <- paste0("cut -f1-5 /home/nwl15/WORKING_DATA/HRCstuff/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | awk '$1==\"#CHROM\" || $1==", chrom, "'")

hrcfilein <- paste0("cut -f1-5 /nobackup/proj/jamgaml/dean_AML/HRCstuff/HRC.r1-1.GRCh37.wgs.mac5.sites.tab | awk '$1==\"#CHROM\" || $1==", chrom, "'")

hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]

## this not holds for chromosome 9
##stopifnot(all(hrc.rsid[, rsid] %chin% chrom.res[, paste0(rsid, "_", alleleA, "_", alleleB)]))

## update rsid if they have one
chrom.res1[hrc.rsid, SNP:=i.ID, on=c(rsid="rsid")]
chrom.res1[grepl("^rs", SNP), rsid:=SNP][, (c("SNP")):=NULL]

######################################
## studyMAF of 2.5% and info of 0.9 ##
######################################

## chrom.res1[studyMAF < 0.025, ':=' (crudeHR=NA,
##                                    crudeHR_L95=NA,
##                                    crudeHR_H95=NA,
##                                    crude.pvalue=NA,
##                                    beta=NA,
##                                    SE=NA,
##                                    clogP=NA)]
## chrom.res1[info_score < 0.9, ':=' (crudeHR=NA,
##                                    crudeHR_L95=NA,
##                                    crudeHR_H95=NA,
##                                    crude.pvalue=NA,
##                                    beta=NA,
##                                    SE=NA,
##                                    clogP=NA)]


setorder(chrom.res1, chromosome, position)


saveRDS(chrom.res1, file=paste0(phenoprefix, "_", dataprefix, "_chrom" , chrom, "_cleaned.rds"))



