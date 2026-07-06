# download sumstat from kar_etal
library(tidyverse)
library(vroom)
library(TwoSampleMR)
library(readxl)
library(data.table)
library(topr)

#################
### initial work
#################
# load aml sumstats and check
pan.aml <- fread('../../AMLmeta_results/unfiltered/status_NCL_PCspeAMLHRC.meta', header = T) %>% mutate(study='Pan-AML')
ch.aml <- fread('../../AMLmeta_results/unfiltered/Normal_NCL_PCspeAMLHRC.meta', header = T) %>% mutate(study='CN-AML')
del57.aml <- fread('../../AMLmeta_results/unfiltered/del5_7_NCL_PCspeAMLHRC.meta', header = T) %>% mutate(study='del57-AML')
com.aml <- fread('../../AMLmeta_results/unfiltered/Complex_NCL_PCspeAMLHRC.meta', header = T) %>% mutate(study='complex-AML')
head(ch.aml)

# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# need to keep aml lead snps for retrospective lookup
aml.lead <- c('rs4665765', 'rs11481', 'rs3916765', 'rs79918355', 'rs12988876', 'rs12078864')

#comb aml
aml.comb <- list(pan.aml, ch.aml, del57.aml, com.aml)
rm(pan.aml, ch.aml, del57.aml, com.aml)

# isolate top hits
library(parallel)
aml.top <- mclapply(aml.comb, function(x){
  n.max <- max(x$N)
  x.sig <- x %>% filter(P < 5e-8 & N == n.max )
  x.sig$rsid <- paste0(x.sig$CHR,':',x.sig$BP,'_',x.sig$A2,'_', x.sig$A1)
  x.sig$snp.id <- paste0(x.sig$CHR,':',x.sig$BP)
  #merge
  x <- merge(x.sig, hrc.rsid[,c(3,6)], by='rsid')
  return(x)
}, mc.cores = 4)

aml.comb <- aml.comb[-c(3,4)] # remove other 2 phenos
#aml comb annotate all
aml.comb <- mclapply(aml.comb, function(x){
  n.max <- max(x$N)
  x <-  x %>% filter(N == n.max)
  x$rsid <- paste0(x$CHR,':',x$BP,'_',x$A2,'_', x$A1)
  x$snp.id <- paste0(x$CHR,':',x$BP)
  #merge
  x <- merge(x, hrc.rsid[,c(3,6)], by='rsid')
  return(x)
}, mc.cores = 2)
head(aml.comb)

##########################
##### load snps of int ###
##########################
# get the snps of int list from papers
# include only CH sig snps from kar (i have included ch,dnm3a and tet 2 sig for kar)and kessler studies 
# then test them in their respective dnmt3a and tet2 trait studies + panaml and cnAML
snps.of.int_kar <- read_xlsx('kar_snps.xlsx', sheet = 1) # this is hg19/37, sheet4 is grch38
snps.of.int_kar$snp.id <- paste0(snps.of.int_kar$CHR,':',snps.of.int_kar$POS)
#remove large, small
snps.of.int_kar <- snps.of.int_kar[snps.of.int_kar$GWAS %in% c('CH', 'DNMT3A', 'TET2'),]
head(snps.of.int_kar)

snps.of.int_kar.cond <- read_xlsx('kar_snps.xlsx', sheet = 2) # only keep 1st snp of this
#snps.of.int_kar.cond <- snps.of.int_kar.cond[1,]
snps.of.int_kar.cond$snp.id <- paste0(snps.of.int_kar.cond$Chr,':',snps.of.int_kar.cond$bp)
snps.of.int_kar.cond <- snps.of.int_kar.cond[snps.of.int_kar.cond$GWAS %in% c('CH', 'DNMT3A', 'TET2'),]
head(snps.of.int_kar.cond)

snps.of.int_Weinstock <- read_xlsx('kar_snps.xlsx', sheet = 3) # Joshua S Weinstock - for rs2887399
snps.of.int_Weinstock$snp.id <- '14:96180695' #14:95714358' #95714358 this is grch38 too
snps.of.int_Weinstock$ref <- 'G'
snps.of.int_Weinstock$rsid <- 'rs2887399'
head(snps.of.int_Weinstock)

#kessler subtypes - keep snps only replicated in ukb and GHS
snps.of.int_kessler <- read_xlsx('kar_snps.xlsx', sheet = 6) # this is grch38
snps.of.int_kessler$hg19 <- sub('chr','', snps.of.int_kessler$hg19)
snps.of.int_kessler$snp.id <- sub('-',':', snps.of.int_kessler$hg19)
#remove snps that were not replicated
snps.of.int_kessler <- snps.of.int_kessler[snps.of.int_kessler$Pval_GHS < 0.05,]
range(snps.of.int_kessler$Pval_GHS)
head(snps.of.int_kessler)
length(unique(snps.of.int_kessler$snp.id))
table(snps.of.int_kessler$Trait)

# #load Bick data
# snps.of.int_bick <- read_xlsx('kar_snps.xlsx', sheet = 9) # this is grch38
# snps.of.int_bick$hg19 <- sub('chr','', snps.of.int_bick$hg19)
# snps.of.int_bick$snp.id <- sub('-',':', snps.of.int_bick$hg19)
# head(snps.of.int_bick)
# #add bick snps manually
# bick.snps <- c('rs34002450', '')

#get hits and their source
snps.source <- data.frame(CHROM=c(snps.of.int_kar$CHR, snps.of.int_kar.cond$Chr), snp.id=c(snps.of.int_kar$snp.id,snps.of.int_kar.cond$snp.id), source='kar et al' )
snps.source <- rbind(snps.source, data.frame(CHROM=snps.of.int_kessler$Chr, snp.id=snps.of.int_kessler$snp.id, source='kesslar et al'))
snps.source <- rbind(snps.source, data.frame(CHROM=14, snp.id=snps.of.int_Weinstock$snp.id, source='Weinstock et al'))
#snps.source <- rbind(snps.source, data.frame(snp.id=snps.of.int_bick$snp.id, source='Bick et al'))
snps.source <- unique(snps.source)
table(snps.source$source)
length(unique(snps.source$snp.id))
saveRDS(snps.source, 'ch_snpstoCheck.Rds')

# for snp.lists, better to have both hg19 and hg38 ids
aml.top <- do.call(rbind, aml.top)

ch.top.list <- c(snps.of.int_kar$snp.id, snps.of.int_kar.cond$snp.id, snps.of.int_kessler$snp.id, snps.of.int_Weinstock$snp.id)
ch.top.list <- unique(ch.top.list)

# check against snp
table(ch.top.list %in% aml.comb[[2]]$snp.id)
#check snps for each input
check.snps <- function(x) {
  table(unique(x$snp.id) %in% aml.comb[[2]]$snp.id)
}

check.snps(snps.of.int_kar)
check.snps(snps.of.int_kar.cond)
check.snps(snps.of.int_kessler)
check.snps(snps.of.int_Weinstock)

#filter and keep those 72 snps in aml
aml.comb <- lapply(aml.comb, function(x){
  x <- x[x$snp.id %in% ch.top.list,]
  return(x)
})

saveRDS(aml.top, 'AML_all_topforCH_snplist.Rds')
saveRDS(aml.comb, 'AML_combforCH_snplist.Rds')

##################
##start from here <<<<<
#################
aml.comb <- readRDS('AML_combforCH_snplist.Rds')
aml.top <- readRDS('AML_all_topforCH_snplist.Rds')
#add hg38 to both aml.top and aml.comb
all.hits <- do.call(rbind, aml.comb)
all.hits <- rbind(all.hits, aml.top)
write.table(unique(all.hits$snp.id), file = 'snps_for_liftover.txt', quote = F, row.names = F, col.names = F)
#add hg38 coords too after doing this on https://genebe.net/tools/liftover
ch38.coords <- unique(read.table('commonSNPS_hg38.txt', header = T))
ch38.coords$hg38 <- sub('chr','', ch38.coords$hg38)
ch38.coords$hg38 <- sub('-',':', ch38.coords$hg38)
head(ch38.coords)
#merge them back to aml.comb and aml.top
aml.comb <- lapply(aml.comb, function(x){
  x <- merge(x, ch38.coords, by.x ='snp.id', by.y = 'Query', all.x = T)
  #x$snp.id.hg38 <- paste0(x$hg38,':',x$A2,':', x$A1) don't need alleles as this could make us lose
  #opposite direction hits. Grep anyway matches the first
  return(x)
  })
aml.top <- merge(aml.top, ch38.coords, by.x ='snp.id', by.y = 'Query', all.x = T)
#aml.top$snp.id.hg38 <- paste0(aml.top$hg38,':',aml.top$A2,':', aml.top$A1)

######################################
#process kessler separately - use hg38
#######################################
#get both aml top and aml comb snps 
all.hits <- do.call(rbind, aml.comb) # update all hits again
all.hits <- rbind(all.hits, aml.top)
all.hits <- all.hits[order(all.hits$CHR)]
head(all.hits)
tmpfile <- tempfile()
writeLines(unique(all.hits$hg38), tmpfile) # use hg38 snps
kes.files <- list.files('sumstats', pattern = 'GRCh38',full.names = T)
cmd <- sprintf("(zcat %s | head -n 1 && zgrep -Fwf %s %s)", kes.files, tmpfile, kes.files)
chip.full <- parallel::mclapply(cmd, fread, mc.cores = 3)
chip.full <- do.call(rbind, chip.full)
# add rsid to chip.full
chip.full$snp.id.hg38 <- paste0(chip.full$chromosome,':',chip.full$base_pair_location)
chip.full <- merge(chip.full, unique(all.hits[,c(1,16:18)]), by.x = 'snp.id.hg38', by.y ='hg38', all.x = T) # need to amend <<
chip.full$case_con <- paste0(chip.full$num_cases,'/',chip.full$num_controls)
setnames(chip.full, 'ID', 'variant_id')
names(chip.full)

#######################
#### add Kar studies
#######################
# process other sumstats
tmpfile <- tempfile()
writeLines(unique(all.hits$ID), tmpfile) # this time its rsids
kar.lst <- list.files('sumstats', pattern = 'GRCh37',full.names = T)
cmd <- sprintf("(zcat %s | head -n 1 && zgrep -Fwf %s %s)", kar.lst, tmpfile, kar.lst)
kar.full <- parallel::mclapply(cmd, fread, mc.cores = 3)
kar.full[[1]]$trait <- 'CH_inclusive'
kar.full[[2]]$trait <- 'DNMT3A'
kar.full[[3]]$trait <- 'TET2'

#add case con nums and nearest gene
kar.full[[1]]$case_con <- '9386/211899'
kar.full[[2]]$case_con <- '6829/211899'
kar.full[[3]]$case_con <- '2818/211899'
#gwas.full[[4]]$case_con <- '25657/342869'
#gwas.full[[5]]$case_con <- '4710/12938'

gc()

#change beta0 to beta
kar.full[[1]] <- kar.full[[1]] %>% mutate(BETA = beta/((10203/(10203+173918))*(1-(10203/(10203+173918)))))
kar.full[[1]] <- kar.full[[1]] %>% mutate(SE = standard_error/((10203/(10203+173918))*(1-(10203/(10203+173918)))))
kar.full[[2]] <- kar.full[[2]] %>% mutate(BETA = beta/((5185/(5185+173918))*(1-(5185/(5185+173918)))))
kar.full[[2]] <- kar.full[[2]] %>% mutate(SE = standard_error/((5185/(5185+173918))*(1-(5185/(5185+173918)))))
kar.full[[3]] <- kar.full[[3]] %>% mutate(BETA = beta/((2041/(2041+173918))*(1-(2041/(2041+173918)))))
kar.full[[3]] <- kar.full[[3]] %>% mutate(SE = standard_error/((2041/(2041+173918))*(1-(2041/(2041+173918)))))
head(kar.full[3])

#calculate ORs for 1-5 - this needs to be corrected beta0 is in the sumstats
kar.full <- lapply(kar.full, function(x){
  odds_ratio <- exp(x[,'BETA'])
  ci_lower <- exp(x[,'BETA'] - 1.96 * x[,'SE']) 
  ci_upper <-  exp(x[,'BETA'] + 1.96 * x[,'SE'])
  x[,'odds_ratio'] <- odds_ratio
  x[,'ci_lower'] <- ci_lower
  x[,'ci_upper'] <- ci_upper
  return(x)
})

# comb all kar
kar.full <- do.call(rbind, kar.full)
names(kar.full)
head(kar.full)
length(unique(kar.full$variant_id))
length(unique(chip.full$name))

#final checks before saving
table(unique(all.hits$ID) %in% chip.full$variant_id)
table(unique(all.hits$ID) %in% kar.full$variant_id)

saveRDS(chip.full, 'Kessler_chip_full.Rds')
saveRDS(kar.full, 'Kar_ch_full.Rds')

### continue in 2.aml_ch script <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

