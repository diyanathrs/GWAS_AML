###################################
## Two sample MR for AML vs CHIP
##################################
library(data.table)
library(TwoSampleMR)

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

####################
# start from there
####################
setwd('../CHIP/ch_new/')
# load pre-processed data
kessler.full <- readRDS('Kessler_chip_full.Rds')
kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
kar.full <- readRDS('Kar_ch_full.Rds')
kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
AML.full <- readRDS('AML_combforCH_snplist.Rds')
#ch.snps.chk <- read.csv('nearest_gene.csv', header = T)

# declare exp and effect
exp <- kessler.full
names(exp)

exposure <- format_data(
  exp, type = "exposure",
  snp_col = "name",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  samplesize_col = "case_con"
)

# clump  with boarderline significant 
exposure_sig <- subset(exposure, pval.exposure < 5e-8)

exposure_sig <- clump_data(
  exposure_sig,
  clump_kb = 10000,
  clump_r2 = 0.001,
  pop = "EUR"
)

# load outcome data
outcome <- format_data(
  AML.full,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "p",
  samplesize_col = "N"
)

#harmonize
outcome <- outcome[outcome$SNP %in% exposure_sig$SNP, ]
dat <- harmonise_data(
  exposure_sig,
  outcome)

mr_ch_aml <- mr(dat)

mr_heterogeneity(dat)
mr_pleiotropy_test(dat) # pleiotropy
loo <- mr_leaveoneout(dat)

mr_leaveoneout_plot(loo)
#single snp effects
single <- mr_singlesnp(dat)

mr_forest_plot(single)

# plot
res <- mr(dat)
plots <- mr_scatter_plot(res, dat)
plots[[1]]

