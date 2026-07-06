## args=(commandArgs(TRUE))
## if(length(args)==0){
##     print("No arguments supplied.")
##     ##supply default values
## } else {
##     for(i in 1:length(args)){
##         eval(parse(text=args[[i]]))
##     }
## }

## generate IID based on IID_IID to match up with ID in VCF
## pheno.mod <- fread("./association/CLL_allphenos.txt", header=T)
## pheno.mod[, IID:=paste0(IID, "_", IID)]
## write.table(pheno.mod, file="./association/CLL_allphenosMOD.txt",
##             quote=F, sep="\t", row.names=F, col.names=T)

require("data.table")
impv.path <- "/nobackup/proj/jamgaml/HRCimpvData/OxfCLL/impvRes"
out.path <- "/nobackup/proj/jamgaml/HRCimpvData/OxfCLL/impvSubset"

prefix <- "OxfCLL4imp"


#chrom <- 2
#pheno
#cov
#pc.num <- as.numeric(pcnum)
#pc.num

for (chrom in seq(1,22)){

    vcf <- paste0(impv.path, "/chr", chrom, ".dose.vcf.gz")
    vcf
    sampleinVCF <- fread(paste0("zcat ", vcf, " | grep -v ^## | head -n 1 | awk '{for (i=1;i<=NF;i++) print $i}' "), header=FALSE)[[1]]
    check.cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    stopifnot(sampleinVCF[1:9]==check.cols)

    sample.temple <- data.table(ID_1=c("0", sampleinVCF[10:length(sampleinVCF)]),
                            ID_2=c("0", sampleinVCF[10:length(sampleinVCF)]),
                            missing=0
                            )[, oriOrder:=1:.N]
    setorder(sample.temple, oriOrder)
    
    write.table(sample.temple[, (c("oriOrder")):=NULL], file=paste0(out.path, "/", prefix, "_chr", chrom, "subsetMOD.sample"),
                quote=F, col.names=T, row.names=F
                )

}
