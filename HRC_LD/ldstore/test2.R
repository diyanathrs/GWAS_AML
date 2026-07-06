
# start from the hits file
require(GenomicRanges)
require(dplyr)
hits <- read.table('AMLmeta_hits.txt',header=T)

hits.grn <- GRanges(seqnames=hits$CHR,ranges=IRanges(hits$BP-1500000,hits$BP+1500000))
hits.grn <- reduce(hits.grn)

#for (i in 1:length(hits.grn)) {
#	chrom <- data.frame(hits.grn[i])$seqnames
#	header.line <- c('rsid', 'chromosome', 'position', 'allele1', 'allele2')
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



