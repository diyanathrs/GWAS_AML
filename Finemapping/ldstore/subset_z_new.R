meta.path <- '/nobackup/proj/jamgaml/dean_AML/HRCimpvData/AMLaetioMeta_forFinemap/NCL1_7/'
meta.files <- list.files(pattern = '.meta.gz',path =meta.path ,full.names = T)
GENOTYPEDIR <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/merged_bgen_new/'
GWAS_SIG_THR <- 5e-8
offset <- 1e6


require(GenomicRanges, quietly=T)
require(dplyr, quietly=T)
hits.lst <- list()
pheno.no <- 0

# write master file and its first line
master.first.line <- "z;bgen;bgi;bcor;ld;n_samples;bdose"
write(master.first.line, file= paste0('HRC_all','.master'), append=FALSE)


## for each pheno - find hits with sig and N=nmax. Filter snps 1.5mb both ways from each hit
## match the pos with bgen.info file and write those as z file. (change snp names to bgen.info names) 
## also add SE, OR etc. from the sum.stats file to the z file

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

	pheno.range <- GRanges(seqnames=GWASData[idx,]$CHR, ranges=IRanges(GWASData[idx,]$BP-offset, GWASData[idx,]$BP+offset))
        pheno.range <- reduce(pheno.range)

	# filter snps based on ranges and write to a z file. chr wise? yes
	chr.lst <- unique(seqnames(pheno.range))
	# second for loop - do for each range in pheno.range, write chr based z file
	for (i in 1:length(pheno.range)) {
		start <- data.frame(pheno.range)$start[i]
		end <- data.frame(pheno.range)$end[i]	
		chrom <- data.frame(pheno.range)$seqname[i]		
		gwas.data.chr <- filter(GWASData,CHR==chrom)
		gwas.data.chr.fil <- filter(gwas.data.chr,between(BP,start,end))
		print(paste('Filtering SNPs between',start,'and',end, 'of chromosome', chrom))
		print(paste(nrow(gwas.data.chr.fil),'SNPs filtered..'))
		
		names(gwas.data.chr.fil) <- c('chromosome','position','rsid','allele2','allele1','N','P','P_R','beta','beta_R','Q','I') # changes alleles
	        #calculate SE(std error) from plink output
        	gwas.data.chr.fil$se <- abs(log(gwas.data.chr.fil$beta)/qnorm(gwas.data.chr.fil$P/2))
		# if Standard error is 0 add a small floating positive value
		gwas.data.chr.fil$se[gwas.data.chr.fil$se==0] <- 0.1 
       		# add dummy maf
        	gwas.data.chr.fil$maf <- 0.1 #can add MAF from bgen.info file?
		#change rsid - 13:50246074,13:50246074:C:T (should be A2:A1)
		#make rsid names according to the info
		gwas.data.chr.fil$rsid <- paste0(gwas.data.chr.fil$chromosome,':',gwas.data.chr.fil$position,',',gwas.data.chr.fil$chromosome,':',gwas.data.chr.fil$position,':',gwas.data.chr.fil$allele1,':',gwas.data.chr.fil$allele2)

		# load corresponding info file and match snps
	#	info.file <- read.table(paste0('../merged_bgen_new/chr',chrom,'_merged.info'),header=T)	
	#	gwas.data.chr.fil$position

		finemap_z_file <- paste0(pheno,'_hit',i,'_chr',chrom,'.z')
        	print(paste('Writing z file for LDstore2 in',finemap_z_file))
        	gwas.data.chr.fil <- gwas.data.chr.fil[, c(3,1,2,5,4,14,9,13,6,7,8,10,11,12,13)]
        	#print(head(gwas.data.chr.fil))
        	write.table(gwas.data.chr.fil, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F) #isn't working
	       
		 #match with the info.chr positions and use rsids of those - last resort

		#Write .master file too
		# line per each .z file but one master file
				
		outbcorfile <- paste0('HRCout_', pheno,'_hit',i,'_chr',chrom,'.bcor')
		finemap_ld_file <- paste0('HRCout_',pheno,'_hit',i,'_chr',chrom, '.ld')
		bgenfile <- paste0(GENOTYPEDIR,'chr',chrom,'_merged.bgen')
		bgenbgifile <- paste0(bgenfile,'.bgi')
		bdosefile <- paste0('HRCout_', pheno,'_hit',i,'_chr',chrom, '.bdose')
		samplecount <- 12938
		currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount,";",bdosefile)
		print(paste('Writing master file..'))
		write(currstr, file= paste0('HRC_all','.master'), append=TRUE)


		# match with the info.chr positions and use rsids of those - last resort
		}  
	
        }  else {
                cat(sprintf("\n\n *** No GWAS significant SNPs (p-value < 5e-8) - exit !!! \n\n"))
        } }


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
