args=(commandArgs(TRUE))
## if(length(args)==0){
##     stop("No arguments supplied.")
## } else {
##     for(i in 1:length(args)){
##         eval(parse(text=args[[i]]))
##     }
## }

#genfileName <- "flipped_chr22_shapteit_phased_Imputed_chunk_stitched_filtered0.8info.gen"

#args

chrom <- args[1]
##chrom <- 22

ncpu <- as.numeric(Sys.getenv("NSLOTS"))
ncpu


genfileName <- paste0("CLL_finalQced4imp_chr", chrom, "subset.gen")
gensample <- paste0("CLL_finalQced4imp_chr", chrom, "subsetMOD.sample")
#phenofile <- "HullB1_OS_Dx_Death_LFU.txt"
#cen.col <- "OS_Dx_Death_LFU_Status"
#time.col <- "OS_Dx_Death_LFU"
phenofile <- args[2]
cen.col <- args[3]
cen.col
time.col <- args[4]
time.col
data.dir <- args[5]
data.dir


cat("chrom", chrom, "data\n")
cat("gen file is", genfileName, "\n")
cat("sample file is", gensample, "\n" )
cat("pheno file is ", phenofile, "\n")
cat("censored variable used is", cen.col, "\n")
cat("time variable used is", time.col, "\n")


#subtype
kkk <- c("survival", "data.table", "doParallel")
lapply(kkk, require, character.only=TRUE)


index <- 0
chunkSize <- 20000

## checking no. columns in gen are the same as 3*no.samples
## samples <- fread(paste0(args[5], "/", gensample), header=TRUE)[-1, ]
## genfile <- file(paste0(args[5], "/", genfileName), open="r")
samples <- fread(paste0(data.dir, "/", gensample), header=TRUE)[-1, ]
genfile <- file(paste0(data.dir, "/", genfileName), open="r")

sapply(strsplit(readLines(genfile, n=1)," "), length)
stopifnot(sapply(strsplit(readLines(genfile, n=1)," "), length)==(samples[,.N]*3+6))
close(genfile)

## the type of columns
cols <- rep("numeric", samples[, .N]*3+6 )
cols[c(1,2,3,5,6)] <- "character"

## read the data for a chunk
con <- file(paste0(data.dir, "/", genfileName), open="r")


dataChunk <- read.table(con,  nrows=chunkSize, skip=0, header=FALSE, colClasses=cols)

source("../alleleDosageSurvMOD.R")
summary.out <- gen2dosageSurv(dataChunk, phenofile, cen.col, time.col, ncpu, chrom)

write.table(summary.out, paste0(cen.col,"_chr", chrom, "_HRCSurRes.out"), append=FALSE, sep="\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

repeat {
        index <- index + 1
        print(paste('Processing rows:', index * chunkSize))
 
        if (nrow(dataChunk) != chunkSize){
                print('Processed all rows!')
                break}

        dataChunk <- read.table(con,  nrows=chunkSize, skip=0, header=FALSE, colClasses=cols)
        summary.out <- gen2dosageSurv(dataChunk, phenofile, cen.col, time.col, ncpu, chrom)
        
        write.table(summary.out, file=paste0(cen.col,"_chr", chrom, "_HRCSurRes.out"),
                    append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
close(con)
       

