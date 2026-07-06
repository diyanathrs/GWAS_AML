args=(commandArgs(TRUE))
if(length(args)==0){
     print("No arguments supplied.")

 } else {
     for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
     }
 }


chrom

##===============
## run LDStore
##===============
ldstore.exec <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64'
GENOTYPEDIR <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/merged_bgen_new/'
GWAS_SIG_THR <- 5e-8

# for Z file, write SNPin each region file in GWAS_Regions. Z file and master should name with the chromosome
master.first.line <- "z;bgen;bgi;bcor;ld;n_samples;bdose"
write(master.first.line, file= paste0('HRC_',chrom,'.master'), append=FALSE)


## have to filter snps on z file to onlyh include 10000 surrounding snps for each hit
#meta.path <- '/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new/'

#meta.files <- list.files(pattern = '.meta.gz',path =meta.path ,full.names = T)



# load snp info
#info.file <- read.table(paste0(GENOTYPEDIR,'chr',chrom,'_merged.info'), header=T)
#info.file <- info.file[c(2,3,4,5,6,14)]
#data.table::setnames(info.file,c("alleleA","alleleB"),c("allele1","allele2"))
#head(info.file)
#write.table(info.file, paste0('chr',chrom,'.z'),quote=F, row.names=F, sep=' ')

##========= write the master file entries
outbcorfile <- paste0('HRCout_chr',chrom, '.bcor')
finemap_z_file <- paste0('chr',chrom,'_sub.z')
finemap_ld_file <- paste0('HRCout_chr', chrom, '.ld')
bgenfile <- paste0(GENOTYPEDIR,'chr',chrom,'_merged.bgen')
bgenbgifile <- paste0(bgenfile,'.bgi')
bdosefile <- paste0('HRCout_chr', chrom, '.bdose')
samplecount <- 12938
currstr <- paste0(finemap_z_file, ";", bgenfile, ";", bgenbgifile, ";", outbcorfile, ";", finemap_ld_file, ";", samplecount,";",bdosefile)
print(currstr)	
write(currstr, file= paste0('HRC_',chrom,'.master'), append=TRUE)
 

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
system(paste(ldstore.exec, "--in-files", paste0('HRC_',chrom,'.master'), "--read-only-bgen --write-text --write-bdose --n-threads", 2, "--memory 32"))
