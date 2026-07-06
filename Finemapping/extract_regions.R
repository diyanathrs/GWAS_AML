##=============
## function to extract GWAS loci
##=============
suppressMessages(require(GenomicRanges))
## significance threshold for GWAS
GWAS_SIG_THR <- 5e-8


Extract_GWAS_Regions <- function(GWASData,  regionfile) {

	OFFSET <- 1000000	

	tempregionfile <- paste0(dirname(regionfile), '/temp_GWAS_Regions.txt')

	## extract GWAS significant SNPs
	n.max <- max(GWASData$N)
	print(paste('n.max is',n.max))	
	idx <- which(as.numeric(GWASData$P) < GWAS_SIG_THR & GWASData$N == n.max)
	if (length(idx) > 0) {
		GWASSigData <- GWASData[idx, ]
	} else {
		cat(sprintf("\n\n *** No GWAS significant SNPs (p-value < 5e-8) - exit !!! \n\n"))
		return()
	}

	## for each significant variant, extract 3 Mb region surrounding it
	## and then merge those regions using bedtools
	## each line of the generated file denotes individual non-overlapping regions
	tempRegionData <- GWASSigData[, c(1, 2, 2)]
	colnames(tempRegionData) <- c('chr', 'start', 'end')
	tempRegionData[, 2] <- tempRegionData[, 2] - (OFFSET / 2)
	tempRegionData[, 3] <- tempRegionData[, 3] + (OFFSET / 2)
	idx <- which(tempRegionData[, 2] < 0)
	if (length(idx) > 0) {
		tempRegionData[idx, 2] <- 0
	}
	tempRegionData <- tempRegionData[order(tempRegionData[,1], tempRegionData[,2]), ]
	write.table(tempRegionData, tempregionfile, row.names=F, col.names=F, sep="\t", quote=F, append=F)

	## bedtools merge
	system(paste0("bedtools merge -i ", tempregionfile, " > ", regionfile))

	## now divide the non-overlapping regions and extract the summary statistics for individual regions
	regiondata <- data.table::fread(regionfile, header=F)
	GWAS_RegionDir <- paste0(dirname(regionfile), '/GWAS_Regions')
	system(paste("mkdir -p", GWAS_RegionDir))
	for (i in 1:nrow(regiondata)) {
		currregiondata <- as.data.frame(regiondata[i, ])
		ov <- Overlap1D(GWASData[, c(1, 2, 2)], currregiondata, boundary=0)
		curroutfile <- paste0(GWAS_RegionDir, '/Region_',pheno,i, '.txt')
		write.table(GWASData[ov$A_AND_B, ], curroutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}

	## remove temporary files
	if (file.exists(tempregionfile)) {
		system(paste("rm", tempregionfile))
	}
}


##=================
## overlapping 1D genomic regions
##=================
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {
	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(ov1[,1])
		ov_idx_file2 <- unique(ov1[,2])		
	} else {
		ov_idx_file1 <- ov1[,1]
		ov_idx_file2 <- ov1[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)
	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)
}



	
## Step 1: extracting GWAS regions
OutDir <- 'LDstore2_out' 
system(paste("mkdir -p", OutDir))

meta.path <- '/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new/'

meta.files <- list.files(pattern = '.meta.gz',path =meta.path ,full.names = T)

Input_GWAS_Region_File <- paste0(OutDir, '/GWAS_Regions.txt')	
#require(data.table)
for (meta in meta.files) {
GWASData <- read.table(meta, header=T )
print(paste('Reading meta file',meta))
pheno <- gsub('/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new//','',meta)
pheno <- gsub('_NCL_PCspeAMLHRC.meta.gz','',pheno)
print(paste('Phenotype is',pheno))
print(head(GWASData))
Extract_GWAS_Regions(GWASData,  Input_GWAS_Region_File)

}

print('Done...')
## Step 2: run LDStore and FINEMAP
#Run_Finemap_GWAS(FINEMAPInpDir, OutDir, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, GENOTYPEDIR, ldstoreexec, finemapexec, samplecount, NUMCAUSALSNP, NUMTHREAD)

## Step 3: Summarize finemap output
#regiondata <- read.table(Input_GWAS_Region_File, header=F, sep="\t", stringsAsFactors=F)
#cat(sprintf("\n Number of GWAS regions: %s ", nrow(regiondata)))	
#Summary_Finemap(OutDir, regiondata, BaseOutDir_FINEMAP_cond, BaseOutDir_FINEMAP_sss, NUMCAUSALSNP, "GWAS")
