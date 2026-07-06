z.files <- list.files(pattern = '\\.z$',path ='.')
finemapdir <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/finemap/'
finemapexec <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64'

#require(GenomicRanges, quietly=T)
#require(dplyr, quietly=T)

# write master file for finemap and its first line
master.first.line <- "z;ld;snp;config;cred;log;n_samples"
write(master.first.line, file= paste0('HRC_finemap','.master'), append=FALSE)


## for each z file and relevent ld files, write one line in master
for (z in z.files) {
        prefix <- gsub('.z','',z)
        print(paste('Writing finemap master for', prefix))
	
    	ld.file <- paste0('HRCout_',prefix,'.ld')
        snpfile <- paste0(finemapdir,prefix,'.snp')
        config <- paste0(finemapdir,prefix,'.config')
        cred <- paste0(finemapdir,prefix, '.cred')
        log <- paste0(finemapdir, prefix, '.log')
	samplecount <- 12938
        currstr <- paste0(z,";",ld.file,";",snpfile,";",config,";",cred,";",log,";",samplecount)
              #  print(paste('Writing master file..'))
                write(currstr, file= paste0('HRC_finemap','.master'), append=TRUE)

       }

print('Running Finemap on z files')
system(paste(finemapexec, "--sss --in-files", "HRC_finemap.master", "--log --n-causal-snps", 10, "--n-threads", 8))
print(paste('Finemap output written to','/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/finemap'))



## for each pheno file - filter 1.5mb both directions of the hits
#for (i in 1:length(hits.grn)) {
#        chrom <- data.frame(hits.grn[i])$seqnames
#        print(paste('Hit no',i,'chromosome',chrom))
#        info <- read.table(paste0('../merged_bgen_new/chr',chrom,'_merged.info'),header=T)
#        #subset info snps based on hits.grn
#        start <- data.frame(hits.grn[i])$start
#        end <- data.frame(hits.grn[i])$end
#        print(paste('Filtering SNPs between',start,'and',end))
#        info.sub <- filter(info, between(position,start,end))
#        print(paste(nrow(info.sub),'SNPs filtered..'))
#        #header.line <-c(rsid', 'chromosome','position', 'allele1', 'allele2')
#        #write(header.line, file=paste0('chr',chrom,'_sub.z'),row.names=F,sep=' ',quotes=F)
	
#        info.sub <- info.sub[c(2,3,4,5,6,14)]
#       data.table::setnames(info.sub,c("alleleA","alleleB"),c("allele1","allele2"))
#        print(head(info.sub))
#        z.file <- paste0('chr',chrom,'_sub.z')
#        write.table(info.sub, z.file, quote=F, row.names=F, col.names=!file.exists(z.file), sep=' ',append=T)

#}

#}
