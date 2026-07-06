library(tidyverse)
#library(vroom)
#library(TwoSampleMR)
library(readxl)
library(data.table)

# for kar study use snps from 
# get the snps of int list from papers
snps.of.int_kar <- read_xlsx('kar_snps.xlsx', sheet = 1) # this is hg19/37, sheet4 is grch38
snps.of.int_kar$snp.id <- paste0(snps.of.int_kar$CHR,':',snps.of.int_kar$POS)
#remove large, small
snps.of.int_kar <- snps.of.int_kar[snps.of.int_kar$GWAS %in% c('CH', 'DNMT3A', 'TET2'),]
head(snps.of.int_kar)

snps.of.int_kar.cond <- read_xlsx('kar_snps.xlsx', sheet = 2)
snps.of.int_kar.cond$snp.id <- paste0(snps.of.int_kar.cond$Chr,':',snps.of.int_kar.cond$bp)
snps.of.int_kar.cond <- snps.of.int_kar.cond[snps.of.int_kar.cond$GWAS %in% c('CH', 'DNMT3A', 'TET2'),]
head(snps.of.int_kar.cond)

snps.of.int_Weinstock <- read_xlsx('kar_snps.xlsx', sheet = 3) # Joshua S Weinstock - for rs2887399
snps.of.int_Weinstock$snp.id <- '14:96180695' #14:95714358' #95714358 this is grch38 too
snps.of.int_Weinstock$ref <- 'G'
snps.of.int_Weinstock$rsid <- 'rs2887399'
head(snps.of.int_Weinstock)

#get top prob from finemap results of kessler
snps.of.int_kessler <- read_xlsx('kar_snps.xlsx', sheet = 6) # this is grch38
# take highest pip snp from each locus ?
#snps.of.int_kessler <- snps.of.int_kessler[snps.of.int_kessler$Prob > 0.95,]
length(unique(snps.of.int_kessler$Locus))
#snps.of.int_kessler$snp.id <- sub(":\\w+:\\w$",'',snps.of.int_kessler$SNP, perl = T)
head(snps.of.int_kessler)
#write.table(snps.of.int_kessler$snp.id, file = 'kes_for_liftover.txt', quote = F, row.names = F, col.names = F)
#get 5th sheet as liftover results
#kes.liftover <- read_xlsx('kar_snps.xlsx', sheet = 5) # liftover output
#snps.of.int_kessler <- merge(snps.of.int_kessler, kes.liftover, by.x = 'snp.id', by.y='Query')
head(snps.of.int_kessler)
#new.snp id
snps.of.int_kessler$hg19 <- sub('chr','',snps.of.int_kessler$hg19)
snps.of.int_kessler$snp.id <- sub('-',':',snps.of.int_kessler$hg19)

snp.list <- c(snps.of.int_kar$snp.id, snps.of.int_kar.cond$snp.id, snps.of.int_kessler$snp.id, snps.of.int_Weinstock$snp.id)
snp.list <- unique(snp.list)

# load aml sumstats and check
pan.aml <- read.table('../../AMLmeta_results/unfiltered/status_NCL_PCspeAMLHRC.meta', header = T)
# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

pan.aml$rsid <- paste0(pan.aml$CHR,':',pan.aml$BP,'_',pan.aml$A2,'_', pan.aml$A1)
pan.aml$snp.id <- paste0(pan.aml$CHR,':',pan.aml$BP)
head(pan.aml)
#merge
pan.aml <- merge(pan.aml, hrc.rsid[,c(3,6)], by='rsid')

# check against snp
table(snp.list %in% pan.aml$snp.id)
pan.aml[pan.aml$snp.id %in% snp.list,]
saveRDS(pan.aml[pan.aml$snp.id %in% snp.list,], 'panAML_snplist.Rds')


#test
#calculate ORs for 1-5
cal.or <- function(beta, se){
  odds_ratio <- exp(beta)
  ci_lower <- exp(beta - 1.96 * se) 
  ci_upper <-  exp(beta + 1.96 * se)
  out <- paste0(round(odds_ratio,2),'(',round(ci_lower,2),'-',round(ci_upper,2),')')
  return(out)
}
beta0 <- 0.00286715
beta <- beta0/((2041/(2041+173918))*(1-(2041/(2041+173918))))
exp(beta)
se <- se0/((2041/(2041+173918))*(1-(2041/(2041+173918))))

cal.or(0.00449409, 7.64954e-4)
