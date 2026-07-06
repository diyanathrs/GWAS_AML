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
#dataprefix <- "gt"
## read meta res
chrom.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, ".meta"))
if (chrom.res[P < 10^-7, .N]!=0){
    hits <- chrom.res[P < 10^-7, ][, pheno:=phenoprefix]
    setwd("./Hitsmeta")
    write.table(hits, file=paste0(phenoprefix, "_NCL_", dataprefix, "_hits.lst"), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)

    folder.lst <- data.table(
        V1=c(rep("/nobackup/proj/jamgaml/HRCimpvData/CLLbystudy/ResultSummary",6),
             "/nobackup/proj/jamgaml/HRCimpvData/CLL_OEE/ResultSummary",
             "/nobackup/proj/jamgaml/HRCimpvData/OxfCLL/ResultSummary"),
        V2=c("Bournemouth", "Cardiff", "Hull1", "Hull2", "Newcastle1", "Newcastle2", "Hull3", "Oxford"),
        V3=c(rep("CLL_finalQced4imp", 6), "CLL_OEEfinalQced4imp", "OxfCLL4imp")
                             )
    ## extract function
    extract.data <- function(x){
        x <- as.numeric(x)
        hits.chr <- hits[, unique(CHR)]
        out.res <- vector("list", length=length(hits.chr))
        for (n in seq_along(hits.chr)){
            out.res[[n]] <- readRDS(paste0(folder.lst[x, V1], "/", folder.lst[x, V2], "_", phenoprefix, "_", folder.lst[x, V3], "_chrom", hits.chr[n], "_cleaned.rds"))[
                paste0(chromosome,":",position,"_", alleleA, "_", alleleB) %chin% hits[CHR==hits.chr[n], SNP] | rsid %chin% hits[CHR==hits.chr[n], SNP], ][,
                   rsid1:=paste0(chromosome,":",position,"_", alleleA, "_", alleleB)][grepl("^esv", rsid), rsid1:=rsid]

        }                               #loop for hits.chr
        out.res <- rbindlist(out.res, use.names=TRUE, fill=TRUE)
    }

    
    study.res <- vector("list", length=folder.lst[, .N])
    names(study.res) <- folder.lst[, V2]
    for (nn in seq_along(study.res)){
       study.res[[nn]] <- extract.data(nn)
    }
    study.res <- rbindlist(study.res, use.names=TRUE, fill=TRUE, idcol=TRUE)
    setnames(study.res, c(".id"), c("study"))
    write.table(study.res, file=paste0(phenoprefix, "_NCL_", dataprefix, "_hitsByStudy.txt"), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=TRUE)
}                                       # 

