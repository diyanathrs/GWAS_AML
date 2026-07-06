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

## read hits from AML eatiology meta
hits <- fread("AML_survHits.txt", header=T)

folder.lst <- fread("awk '{print $1, $2, $3}' mystudy.lst", header= FALSE)
    setnames(folder.lst, c("V1", "V2"), c("V2", "V1"))
    folder.lst[, V1:=sub("/GWAStabix", "/ResultSummary_noAPL", V1)]
    
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
    write.table(study.res, file="AML_hitsByStudy.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
 

