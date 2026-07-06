meta.path <- '/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new/'
meta.files <- list.files(pattern = '.meta.gz',path =meta.path ,full.names = T)
GWAS_SIG_THR <- 5e-8


#require(data.table)
hits.lst <- list()
pheno.no <- 0
for (meta in meta.files) {
        pheno.no <- pheno.no + 1
        GWASData <- read.table(meta, header=T )
        print(paste('Reading meta file',meta))
        pheno <- gsub('/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new//','',meta)
        pheno <- gsub('_NCL_PCspeAMLHRC.meta.gz','',pheno)
        print(paste('Phenotype is',pheno))
        print(head(GWASData))

        n.max <- max(GWASData$N)
        print(paste('n.max is',n.max))
        idx <- which(as.numeric(GWASData$P) < GWAS_SIG_THR & GWASData$N == n.max)
        if (length(idx) > 0) {
                GWASSigData <- GWASData[idx, ]
                GWASSigData$pheno <- paste(pheno)
                hits.lst[[pheno.no]] <- GWASSigData
                print(paste(length(idx),'GWAS hits above the threashold detected..'))
        } else {
                cat(sprintf("\n\n *** No GWAS significant SNPs (p-value < 5e-8) - exit !!! \n\n"))
        }

}
hits <- do.call(rbind, hits.lst)
print(hits)
write.table(hits,'AMLmeta_hits.txt',sep='\t',quote=F,row.names=F)

#subset 1mb boundary from the snps
require(GenomicRanges)

hits <- read.table('AMLmeta_hits.txt',header=T)

