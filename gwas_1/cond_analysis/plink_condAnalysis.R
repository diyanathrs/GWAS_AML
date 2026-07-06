args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
} else {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

print("Input Cov file and Plink file..")
input 
bfile 
proj <- gsub("_.*","",input)
proj <- gsub("...*/","",proj)
#snp <- '11:67931761_G_A' 
snp <- '11:67820335_T_A'

pheno <- read.table(input, header=T)
head(pheno)
#remove first col and convert 0>1 and 1>2 an NA > 0
pheno$status <- as.numeric(pheno$status)
table(pheno$status)
pheno$status <- pheno$status + 1 
pheno$status[is.na(pheno$status)] <- '0'
table(pheno$status)
pheno <- pheno[-1]
names(pheno)[1] <- 'IID'
head(pheno)

write.table(pheno, file = 'phenoCov.txt', quote = F, row.names = F, col.names = T)
# run plink
# input params
cov.idx <- paste(3:(length(pheno)-1),collapse = ',')

print(paste("Running Plink conditional analysis on", snp))
plink.cmd <- paste0("plink2 --bfile ", bfile, " --pheno phenoCov.txt ","--pheno-name status ", "--condition ", snp, 
                    " --covar-col-nums ",cov.idx, " --glm hide-covar --out ", proj,'_',snp,"_cond.out" )

system(plink.cmd)
