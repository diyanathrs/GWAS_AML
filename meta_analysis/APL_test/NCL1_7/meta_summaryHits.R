args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}


folderlist

require("data.table")
#phenoprefix <- "status"
#dataprefix <- "gt"
## read meta res
chrom.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, ".meta"))
if (chrom.res[P<10^-6, .N]!=0){

    folder.lst <- fread(paste0("../", folderlist), header=FALSE)
    hits <- chrom.res[P<10^-6, ][, pheno:=phenoprefix]

    setwd("./Hitsmeta")
    write.table(hits, file=paste0(phenoprefix, "_NCL_", dataprefix, "_hits.lst"), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
    ## folder list
    ## folder.lst <- fread("/nobackup/proj/jamgaml/HRCimpvData/metaNCL1_4/NCLdataset_ukb500k.lst", header=FALSE)
    

    ## if (dataprefix=="gt"){
    ##     folder.lst <- folder.lst[1:3, ]
    ## } else {
    ##     folder.lst <- folder.lst[c(1:2, 4),]
    ## }
    
    ## extract function
    extract.data <- function(x){
        x <- as.numeric(x)
        hits.chr <- hits[, unique(CHR)]
        out.res <- vector("list", length=length(hits.chr))
        for (n in seq_along(hits.chr)){
            out.res[[n]] <- readRDS(paste0(folder.lst[x, V1], "/", phenoprefix, "_", folder.lst[x, V2], "_chrom", hits.chr[n], "_cleaned.rds"))[
                paste0(chromosome,":",position,"_", alleleA, "_", alleleB) %chin% hits[CHR==hits.chr[n], SNP] | rsid %chin% hits[CHR==hits.chr[n], SNP], ][,
                   rsid1:=paste0(chromosome,":",position,"_", alleleA, "_", alleleB)][grepl("^esv", rsid), rsid1:=rsid]

        }                               #loop for hits.chr
        out.res <- rbindlist(out.res, use.names=TRUE, fill=TRUE)
    }
    
    study.res <- vector("list", length=nrow(folder.lst))
    names(study.res) <- folder.lst[, V2]
    for (nn in seq_along(folder.lst[,V2])){
       study.res[[nn]] <- extract.data(nn)
    }
    study.res <- rbindlist(study.res, use.names=TRUE, fill=TRUE, idcol=TRUE)
    setnames(study.res, c(".id"), c("study"))
    write.table(study.res, file=paste0(phenoprefix, "_NCL_", dataprefix, "_hitsByStudy.txt"), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
}                                       # 

