args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}


## sample by phenoname and chrom

impv.path
prefix
##chrom
pheno
cov

phenoname

if (phenoname=="t15_17"){
    phenoname <- c("t.15.17.")
}

if (phenoname=="del5_7"){
    phenoname <- c("del5.del7")
}


##pc.num <- as.numeric(pcnum)
##pc.num


#q(save="no")

require("data.table")

## plink status 2 for cases, 1 for controls, 0 for missing
pheno.data <- fread(pheno, header=TRUE)
## accounting for discrepant column names
##if (sum(names(pheno.data) %chin% c("Monosomal.karyotype", "Any.monosomy"))!=0){
##    setnames(pheno.data, c("Monosomal.karyotype", "Any.monosomy"), c("monosomal.karyotype", "Any.Monosomy"))
##}
## cov.data <- fread(cov, header=TRUE)

cov.data <- fread(paste0(phenoname, cov), header=TRUE)

if ("Monosomal.karyotype" %in% names(pheno.data)){
    setnames(pheno.data, c("Monosomal.karyotype"), c("monosomal.karyotype"))
}

if ("Any.monosomy" %in% names(pheno.data)){
    setnames(pheno.data, c("Any.monosomy"), c("Any.Monosomy"))
}



##pheno.cols <- c("status", "Normal", "CBF", "Trans", "Complex", "t.15.17.", "del5.del7", "Trisomies", "Any.Monosomy", "monosomal.karyotype")
##pcs.cols <- paste0("PC",seq(1, pc.num))

pheno.cols <- phenoname

## read no.pcs by pheno
pcsbypheno <- fread("noPCsbyphenoTW.txt", header = TRUE)
pc.num <- pcsbypheno$nopcs[pcsbypheno$pheno == pheno.cols]
##pc.num1 <- unlist(strsplit(pc.num, ","))
##pcs.cols <- paste0("PC", seq(1, pc.num))
pcs.cols <- paste0("PC", pc.num)

if (pcs.cols=="PC0"){
    pcs.cols <- NULL
}

stopifnot(all(pheno.cols %chin% names(pheno.data)), all(pcs.cols %chin% names(cov.data)))



pheno.cov <- merge(pheno.data[, c("FID", "IID", pheno.cols), with=FALSE],
                   cov.data[, c("FID", "IID", pcs.cols), with=FALSE],
                   by=c("FID","IID")
                    )
## converting 2 to 1, 1 to 0, and 0 to missing for SNPTEST binary trait
pheno.cov[, (pheno.cols):=lapply(.SD, function(x) x-1), .SDcols=pheno.cols]
stopifnot(all(pheno.cov[[pheno.cols]]>=0))

if (phenoname=="t.15.17."){
    setnames(pheno.cov, c("t.15.17."), c("t15_17"))
    pheno.cols <- "t15_17"
}

if (phenoname=="del5.del7"){
    setnames(pheno.cov, c("del5.del7"), c("del5_7"))
    pheno.cols <- "del5_7"
}

##row.list <- apply(pheno.cov[, pheno.cols, with=FALSE], 2, function(x) which(x==-1))
##cols <- as.integer(which(names(pheno.cov) %in% pheno.cols))
##stopifnot( names(row.list)==names(pheno.cov)[cols])
##for (i in seq_along(row.list)){
##    set(pheno.cov, i=row.list[[i]], j=cols[i], value=NA)}


##
pheno.cov[, ':=' (ID_1=paste(FID, IID, sep="_"), ID_2=paste(FID, IID, sep="_"))][, (c("FID","IID")):=NULL]


for (chrom in seq(1, 22)){

    vcf <- paste0(impv.path, "/chr", chrom, ".dose.vcf.gz")
    vcf

    sampleinVCF <- fread(paste0("zcat ", vcf, " | grep -v ^## | head -n 1 | awk '{for (i=1;i<=NF;i++) print $i}' "), header=FALSE)[[1]]
                                          
    check.cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    stopifnot(sampleinVCF[1:9]==check.cols)
    sample.temple <- data.table(ID_1=c("0", sampleinVCF[10:length(sampleinVCF)]),
                            ID_2=c("0", sampleinVCF[10:length(sampleinVCF)]),
                            missing=0
                            )[, oriOrder:=1:.N]

    sample.file.pheno <- merge(sample.temple, pheno.cov, by=c("ID_1", "ID_2"), all.x=TRUE)

#firrow <- data.table(status=0, Normal=0, CBF=0, Trans=0, Complex=0, t.15.17.=0,
#                     del5.del7=0,Trisomies=0, Any.Monosomy=0, monosomal.karyotype=0, PC1=0, PC2=0, PC3=0, ID_1="0", ID_2="0", missing=0)
#sample.file.pheno <- rbindlist(list(firrow, pheno.cov), use.names=TRUE, fill=TRUE)


    ## convert to character type
    sample.file.pheno[,(c(pcs.cols, pheno.cols)):=lapply(.SD, as.character),
                      .SDcols=c(pcs.cols, pheno.cols)]


    ## pc columns
    ## cols <- as.integer(which(names(sample.file.pheno) %in% pcs.cols))
    cols <- as.integer(match(pcs.cols, names(sample.file.pheno)))
    names(sample.file.pheno)[cols]
    for (i in seq_along(cols)){set(sample.file.pheno, i=1L, j=cols[i],value="C") }

    ## pheno columns
    ## cols <- as.integer(which(names(sample.file.pheno) %in% pheno.cols))
    cols <- as.integer(match(pheno.cols, names(sample.file.pheno)))
    names(sample.file.pheno)[cols]
    for (i in seq_along(cols)){set(sample.file.pheno, i=1L,j=cols[i],value="B") }


    ### setnames(sample.file.pheno, c("t.15.17.", "del5.del7"), c("t15_17", "del5_7"))
    setorder(sample.file.pheno, oriOrder)

    if (pheno.cols=="t.15.17."){
        setnames(sample.file.pheno, c("t.15.17."), c("t15_17"))
        pheno.cols <- "t15_17"
    }
    
    if (pheno.cols=="del5.del7"){
        setnames(sample.file.pheno, c("del5.del7"), c("del5_7"))
        pheno.cols <- "del5_7"
    }
    
    ##cols.output <- c("ID_1", "ID_2", "missing", pcs.cols,  c("status", "Normal", "CBF", "Trans", "Complex", "t15_17", "del5_7", "Trisomies", "Any.Monosomy", "monosomal.karyotype"))
    cols.output <- c("ID_1", "ID_2", "missing", pcs.cols, pheno.cols)

    write.table(sample.file.pheno[, cols.output, with=FALSE], file=paste0(pheno.cols, prefix, "_chr", chrom, "_VCFphenocov.sample"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
}
