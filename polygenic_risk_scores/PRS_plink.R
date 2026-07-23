# Simple internal PRS calculation using PLINK
# get cojo selected snps for panAML
library(dplyr)
library(data.table)

setwd('cond_out_nofiltr/')

# load panAML cojo
pheno <- 'status'
# load cond results
# only use cond_c
cond.lst <- list.files('cojo_slct/', pattern=paste0(pheno,'.*jma'), full.names = T)
cojo.comb <- rbindlist(lapply(cond.lst, fread))
range(cojo.comb$p)
length(cojo.comb$p)
head(cojo.comb)

#correst SNP string
cojo.comb <- cojo.comb %>% tidyr::separate(SNP, into = c("discard", "snp"), sep = ",") %>%
  tidyr::separate(snp, into = c("chr", "pos","A1","A2"), sep = ":") 

cojo.comb$SNP <- paste0(cojo.comb$chr,':',cojo.comb$pos,'_',cojo.comb$A1,'_',cojo.comb$A2)
names(cojo.comb)

score <- cojo.comb[, c("SNP","refA","b")]

colnames(score) <- c("ID","A1","BETA")

write.table(score, "PanAML_score.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

# to PLINK in HPC

# get pheno from linux station 
pheno.lst <- list.files('../pheno/', pattern = 'txt', full.names = T)
pheno <- rbindlist(lapply(pheno.lst, function(x){
  dat <- fread(x)
  dat <- dat %>% select('FID', 'IID', 'status')
}))

dat <- fread(pheno.lst[2]) 
names(pheno)

##################################
### get polygenic scores from HPC
###################################
files <- list.files(path = 'sscores/',pattern="*.sscore", full.names = T)
prs <- rbindlist(lapply(files, fread))

# standardize PRS
prs$PRS <- scale(prs$SCORE1_SUM)
length(prs$IID)
