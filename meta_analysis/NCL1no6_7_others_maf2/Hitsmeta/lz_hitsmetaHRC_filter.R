args=(commandArgs(TRUE))

require("data.table")

##pheno <- "status"
##dataprefix <- "AML1_5tidyup"

pheno <- args[1]
pheno[]

dataprefix <- 'PCspeAMLHRC'
dataprefix[]

##resdir <- "/nobackup/proj/jamgaml/AML_meta/metaNCL1_5tidy/Hitsmeta"

resdir <- getwd()

hits <- fread(paste0(resdir, "/", pheno, "_NCL_", dataprefix, "_hits.lst"), header = TRUE)
hits <- hits[N>=3, ]

fwrite(hits, file=paste0(pheno, "_NCL_", dataprefix, "filtered_hits.lst"),
                quote = FALSE, sep="\t",
                row.names = FALSE, col.names = TRUE )

quit(save='no')

if (length(hits)==0){quit(save="no")}

for (n in seq_along(1:nrow(hits))){
    region.chr <- hits[n, CHR]
    region.start <- max(1, hits[n, BP] - 250000)
    region.end <- hits[n, BP] + 250000
    tabix.res <- paste0(pheno, hits[n, SNP], "query.res")
    tabix.cmd <- paste0("tabix -hf ../metaTabix/", pheno, "_", dataprefix, "assoc_HRC.gz ",
                        paste0(region.chr, ":", region.start, "-", region.end, " > ", tabix.res)
                        )
    system(tabix.cmd)
    lzdsn1 <- fread(tabix.res, header = TRUE)#[N==4, ]

    lzdsn1[, rsid1 := SNP]
    lzdsn1[, SNP:= rsid]
    lzdsn1[!grepl("^rs", SNP), SNP:=paste0("chr", `#CHR`, ":", BP)]
    snp <- lzdsn1[rsid1 %in% hits[n, SNP], SNP]

    fwrite(lzdsn1, file=paste0(pheno, "_", dataprefix, "_", snp, "_lz.txt"),
                quote = FALSE, sep="\t",
                row.names = FALSE, col.names = TRUE
            )
    fwrite(data.table(pheno, dataprefix, snp, region.chr, hits[n, BP]),
                file=paste0(pheno, "_", dataprefix, "_", snp, "_lzname.txt"),
                sep="\t", quote= FALSE, row.names=FALSE, col.names = FALSE
                
                )

    system(paste0("rm ", paste0(c(tabix.res), collapse=" ")))

}


