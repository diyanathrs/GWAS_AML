##=============
## Run polyfin within R 
##=============
# use z files generated earlier to get boundries and chr for each polyfun run
# suppressMessages(require(GenomicRanges))

z.file.path <- '/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/ldstore'
z_files <- list.files(pattern= '\\.z', path= z.file.path, full.names=F)
print(length(z_files))

for (z in z_files) {
	z.data <- read.table(paste0("/mnt/storage/nobackup/proj/jamgaml/dean_AML/HRCimpvData/HRC_LD/ldstore/",z), header=T )
	cat("\n")
	#print(paste('Reading z file',z))
	start <- min(z.data$position)
	end <- max(z.data$position)
	chr <- unique(z.data$chromosome)
	name <- gsub('.z','',z)
	pheno <- gsub('_hit.*z','',z)
	print("##################################################")
	print(paste("#### Running Polyfun for",name,"######"))
	print("##################################################")
	#print(paste("Running polyfun for",name ))



susie.cmd <- paste0("python ../finemapper.py", 
   " --geno ../../merged_bgen_new/HRC_chr.",chr, 
   " --sumstats ",pheno,"_munged.parquet", 
   " --n 13000", 
   " --chr ",chr, 
   " --start ",start,
   " --end ",end, 
   " --method susie",
   " --non-funct",   
   " --max-num-causal 5", 
   " --cache-dir LD_cache", 
   " --out results_nofun/",name,"_susie_nofun.out")
#   " --no-sort-pip",
#   " --susie-outfile results_out/susie_files")

finemap.cmd <- paste0("python ../finemapper.py", 
  "  --geno ../../merged_bgen_new/HRC_chr.",chr, 
  "  --sumstats ",pheno,"_munged.parquet", 
  "  --n 13000",
  "  --chr ",chr,
  "  --start ",start,
  "  --end ",end,
  "  --method finemap", 
  "  --max-num-causal 5",
  " --non-funct",
  "  --cache-dir LD_cache",
  "  --ldstore2 /mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64",
  "  --out results_nofun/",name,"_finemap_nofun.out",
  "  --finemap-exe /mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64")
#  "  --finemap-dir results_out/finemap_files/",name)


print(paste("Running SUSIE for",pheno, "-" ,"chr",chr,":",start,"-",end))
print(susie.cmd)
system(susie.cmd)
cat("\n")
print(paste("Running FINEMAP for",pheno, "-" ,"chr",chr,":",start,"-",end))
print(finemap.cmd)
system(finemap.cmd)
cat("\n")
 
}



print('Done...')
## Step 2: run LDStore and FINEMAP
#Run_Finemap_GWAS(FINEMAPInpDir, OutDir, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, GENOTYPEDIR, ldstoreexec, finemapexec, samplecount, NUMCAUSALSNP, NUMTHREAD)

## Step 3: Summarize finemap output
#regiondata <- read.table(Input_GWAS_Region_File, header=F, sep="\t", stringsAsFactors=F)
#cat(sprintf("\n Number of GWAS regions: %s ", nrow(regiondata)))	
#Summary_Finemap(OutDir, regiondata, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, "GWAS")
