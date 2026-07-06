# prep ldstore masterfile and run

##===============
## run LDStore
##===============
ldstore.exec <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64'
GENOTYPEDIR <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/merged_bgen/'
FINEMAPInDir <- 'finemap'

system(paste("mkdir -p ", FINEMAPInDir))

# for Z file, write SNPin each region file in GWAS_Regions. Z file should name with the chromosome and pheno so
# it can be run with matching masterfile line (chr)
# ex - status_region1_chr2.z;HRC_merged_chr2.bgen;same.bgi;out
#for loop

# pull GWASregions list
gwas.regions <- list.files(path="LDstore2_out/GWAS_Regions", pattern="*.txt", full.names = T)
print(paste("Detected gwas region files.."))
print(gwas.regions)

master.first.line <- "z;bgen;bgi;bcor;ld;n_samples;bdose"
write(master.first.line, file= paste0('LDstore2_out', '/HRC_LD.master'), append=FALSE)

for (gwasfile in gwas.regions) {
	currgwasdata <- data.table::fread(gwasfile, header=T)
	cat(sprintf("\n\n==>> processing gwasfile : %s Number of GWAS entries : %s ", gwasfile, nrow(currgwasdata)))
 	pheno <- gsub('LDstore2_out/GWAS_Regions/Region_','',gwasfile)
	pheno <- gsub('.txt','',pheno)
	print(paste("\n",'Phenotype is',pheno))
	
	finemap_z_file <- paste0('LDstore2_out/', 'Region', pheno, '.z') # going to be the same for Ldstore2 z file
	finemapdf <- currgwasdata
	names(finemapdf) <- c('chromosome','position','rsid','allele1','allele2','N','P','P_R','beta','beta_R','Q','I')
	#calculate SE(std error) from plink output
	finemapdf$se <- abs(log(finemapdf$beta)/qnorm(finemapdf$P/2))
	
	# add dummy maf
	finemapdf$maf <- 0.1
	print(paste('Writing z file for LDstore2 in',finemap_z_file))
	finemapdf <- finemapdf[, c(3,1,2,4,5,14,9,13,6,7,8,10,11,12)]
	print(head(finemapdf))
	write.table(finemapdf, finemap_z_file, row.names=F, col.names=T, sep=" ", quote=F, append=F)	
	
	#get the chr of currgwasdata so we can write the masterfile line
	chrom <- unique(finemapdf$chromosome)
	bgenfile <- paste0(GENOTYPEDIR,'HRC_merged_chr',chrom,'.bgen')
	bgenbgifile <- paste0(GENOTYPEDIR, 'HRC_merged_chr', chrom, '.bgen.bgi')
	bimfile <- paste0(GENOTYPEDIR, '../HRC_chr', chrom, '.bim')

	## Each line of the masterfile should be appended into a list
	Masterfile <- paste0('LDstore2_out', '/HRC_LD.master')
	
	##========= write the master file entries
	outbcorfile <- paste0(FINEMAPInDir, '/Region', pheno, '.bcor')
	finemap_ld_file <- paste0(FINEMAPInDir, '/Region', pheno, '.ld')
	bdosefile <- paste0(FINEMAPInDir, '/Region', pheno, '.bdose')
	samplecount <- 18140
	currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount,";",bdosefile)
	write(currstr, file=Masterfile, append=TRUE)
 	}

# HRC insample hardcalls
#bgenfile <- paste0(GENOTYPEDIR, '/HRC_merged_chr', chrom, '.bgen')
#bgenbgifile <- paste0(GENOTYPEDIR, '/HRC_merged_chr', chrom, '.bgen.bgi')
#bimfile <- paste0(GENOTYPEDIR, '../HRC_chr', chrom, '.bim')

##========= write the master file entries
#outbcorfile <- paste0(FINEMAPInDir, '/Region', i, '.bcor')
#finemap_ld_file <- paste0(FINEMAPInDir, '/Region', i, '.ld')
#bdosefile <- paste0(FINEMAPInDir, '/Region', i, '.bdose')
#currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount, ";", bdosefile)
#write(currstr, file=Masterfile, append=TRUE)


## now execute ldstore using the generared master file
system(paste(ldstore.exec, "--in-files", Masterfile, "--read-only-bgen --write-text --write-bdose --n-threads", 12, "--memory 40"))
