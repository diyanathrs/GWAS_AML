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

res.dir
phenoprefix
dataprefix

snptest.cols <- c("rsid", "chromosome", "position", "alleleA", "alleleB", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1", "info_score", "cases_maf", "controls_maf")

chroms <- seq(1, 22)
chrom.res <- vector("list", length=22)
for (n in seq_along(chroms)){
    filein <- paste0(res.dir, "/", phenoprefix, "_", dataprefix, "_chrom" , chroms[n], "_cleaned.rds")
    ## plink can't cope with extremely long string of SNP names
    ## nchar(rsid1)>150, all of them are CNV
    chrom.res[[n]] <- readRDS(filein)[, snptest.cols, with=FALSE][,
            rsid1:=paste0(chromosome,":",position,"_", alleleA, "_", alleleB)][grepl("^esv", rsid), rsid1:=rsid]
}

chrom.res <- rbindlist(chrom.res, use.names=TRUE, fill=TRUE)


chrom.res[info_score <= 0.6,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

chrom.res[cases_maf <= 0.0095,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

chrom.res[controls_maf <= 0.0095,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

## rsid1 should be unique 
exclude.out <- chrom.res[, .N, by=rsid1][N>1, rsid1]
if (length(exclude.out)!=0){
    write.table(chrom.res[rsid1 %chin% exclude.out, ],
                file=paste0(phenoprefix, "_", dataprefix,"_excluded.out"),
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE )}
## 
setnames(chrom.res,
         c("rsid1", "chromosome", "position", "alleleA", "alleleB", "frequentist_add_pvalue", "frequentist_add_beta_1", "frequentist_add_se_1"),
         c("SNP", "CHR", "BP", "A2", "A1", "P", "BETA", "SE")
         )

write.table(chrom.res[!(SNP %chin% exclude.out),  c("SNP", "CHR", "BP", "A2", "A1", "BETA", "SE", "P"), with=FALSE],
            sep="\t",
            file=paste0(phenoprefix, "_", dataprefix,"_assoc.res"), quote=FALSE,
            row.names=FALSE, col.names=TRUE)


