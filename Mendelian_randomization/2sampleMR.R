###################################
## Two sample MR for AML vs CHIP
##################################
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

kar.ch <- vroom('../CHIP/ch_new/sumstats/GCST90102618_buildGRCh37.tsv.gz')
names(kar.ch)
#correct beta
kar.ch <- kar.ch %>% mutate(BETA = beta/((10203/(10203+173918))*(1-(10203/(10203+173918)))))
kar.ch <- kar.ch %>% mutate(SE = standard_error/((10203/(10203+173918))*(1-(10203/(10203+173918)))))

kar.dnmt3a <- vroom('../CHIP/ch_new/sumstats/GCST90102619_buildGRCh37.tsv.gz')
kar.dnmt3a <- kar.dnmt3a %>% mutate(BETA = beta/((5185/(5185+173918))*(1-(5185/(5185+173918)))))
kar.dnmt3a <- kar.dnmt3a %>% mutate(SE = standard_error/((5185/(5185+173918))*(1-(5185/(5185+173918)))))

kar.tet2 <- vroom('../CHIP/ch_new/sumstats/GCST90102620_buildGRCh37.tsv.gz')
kar.tet2 <- kar.tet2 %>% mutate(BETA = beta/((2041/(2041+173918))*(1-(2041/(2041+173918)))))
kar.tet2 <- kar.tet2 %>% mutate(SE = standard_error/((2041/(2041+173918))*(1-(2041/(2041+173918)))))

kes.CH <- vroom('../CHIP/ch_new/sumstats/GCST90165267_buildGRCh38.tsv.gz')
kes.dnmt3a <- vroom('../CHIP/ch_new/sumstats/GCST90165271_buildGRCh38.tsv.gz')
kes.tet2 <- vroom('../CHIP/ch_new/sumstats/GCST90165281_buildGRCh38.tsv.gz')

names(kar.dnmt3a)

# load panAML and annotate SNPs
pan.aml <- vroom('../AMLmeta_results/unfiltered/status_NCL_PCspeAMLHRC.meta')
head(pan.aml)
# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

pan.aml$rsid <- paste0(pan.aml$CHR,':',pan.aml$BP,'_',pan.aml$A2,'_', pan.aml$A1)
pan.aml$snp.id <- paste0(pan.aml$CHR,':',pan.aml$BP)
head(pan.aml)
#merge
pan.aml <- merge(pan.aml, hrc.rsid[,c(3,6)], by='rsid')

# check snps across gwas
table(pan.aml$ID %in% kar.ch$variant_id)

# get beta from OR
pan.aml$Beta <- log(pan.aml$OR)
pan.aml$Z_score <- abs(qnorm(pan.aml$P / 2))
pan.aml$SE_Beta <- abs(pan.aml$Beta / pan.aml$Z_score)

table(duplicated(pan.aml$rsid))

aml.out <- format_data(
  data.frame(pan.aml), type="outcome",
  snp_col="ID", 
  beta_col="Beta",
  se_col="SE_Beta",
  effect_allele_col="A1",
  other_allele_col="A2",
  pval_col="P")

##############
# start MR
##############
run_mr <- function(exposure,
                   outcome,
                   trait_name){
  
  ###########################
  # load exposure sumstats
  ###########################
  # process exp
  chip.exp <- TwoSampleMR::format_data(
    exposure,
    type="exposure",
    snp_col="variant_id",
    beta_col="BETA",
    se_col="SE",
    effect_allele_col="effect_allele",
    other_allele_col="other_allele",
    eaf_col="effect_allele_frequency",
    pval_col="p_value")
  
  #independent instruments
  chip.exp <- subset(chip.exp, pval.exposure < 5e-7)
  # prune ld
  #check API token first
  ieugwasr::user()
  ieugwasr::get_opengwas_jwt()
  
  chip.exp <- TwoSampleMR::clump_data( chip.exp, clump_r2=0.001,
                                       clump_kb=10000, pop = 'EUR')
  
  # harmonize
  dat <- harmonise_data(chip.exp, aml.out)
  
  # run MR
  mr.results <- mr(dat)
  
  mr.results
  mr_heterogeneity(dat)
  mr_pleiotropy_test(dat)
  loo <- mr_leaveoneout(dat)
  
  # plots
  p1 <- mr_scatter_plot(mr.results, dat)
  p1 
  # forest
  single <- mr_singlesnp(dat)
  forest <- mr_forest_plot(single)
  
  #loo
  mr_leaveoneout_plot(loo)
  funnel <- mr_funnel_plot(single)
  
  # odds ratio table
  generate_odds_ratios(mr.results)
  directionality_test(dat)
  
  library(MRPRESSO)
  dat$F <- (dat$beta.exposure^2) / (dat$se.exposure^2)
  mean(dat$F)
  
  
  return(list(mr.results,loo,p1, forest, funnel))
}




