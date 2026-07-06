meta.path <- '/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new/'
meta.files <- list.files(pattern = '.meta.gz',path =meta.path ,full.names = T)
GWAS_SIG_THR <- 5e-8

if (!file.exists('AMLmeta_hits.txt')) {
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

require(dplyr)
hits <- do.call(rbind, hits.lst)
hits <- hits %>% arrange(CHR,BP)
print(hits)
write.table(hits,'AMLmeta_hits.txt',sep='\t',quote=F,row.names=F)

} else {

print('AML hits file already exists..')
# start from the hits file
require(GenomicRanges)
hits <- read.table('AMLmeta_hits.txt',header=T)

# start from the hits file
require(dplyr)
#hits <- read.table('AMLmeta_hits.txt',header=T)

hits.grn <- GRanges(seqnames=hits$CHR,ranges=IRanges(hits$BP-1000000,hits$BP+1000000))
hits.grn <- reduce(hits.grn)

#for (i in 1:length(hits.grn)) {
#       chrom <- data.frame(hits.grn[i])$seqnames
#       header.line <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2')
#        write.table(header.line, file=paste0('chr',chrom,'_sub.z'),row.names=F,sep=' ',quotes=F)
#}

for (i in 1:length(hits.grn)) {
        chrom <- data.frame(hits.grn[i])$seqnames
        print(paste('Hit no',i,'chromosome',chrom))
        info <- read.table(paste0('../merged_bgen_new/chr',chrom,'_merged.info'),header=T)
        #subset info snps based on hits.grn
        start <- data.frame(hits.grn[i])$start
        end <- data.frame(hits.grn[i])$end
        print(paste('Filtering SNPs between',start,'and',end))
        info.sub <- filter(info, between(position,start,end))
        print(paste(nrow(info.sub),'SNPs filtered..'))
        #header.line <-c(rsid', 'chromosome','position', 'allele1', 'allele2')
        #write(header.line, file=paste0('chr',chrom,'_sub.z'),row.names=F,sep=' ',quotes=F)

        info.sub <- info.sub[c(2,3,4,5,6,14)]
        data.table::setnames(info.sub,c("alleleA","alleleB"),c("allele1","allele2"))
        print(head(info.sub))
        z.file <- paste0('chr',chrom,'_sub.z')
        write.table(info.sub, z.file, quote=F, row.names=F, col.names=!file.exists(z.file), sep=' ',append=T)

}

}
