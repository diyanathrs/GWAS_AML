# kar and kessler hits are saved
# also top hits and aml hits are saved
library(data.table)
library(dplyr)
library(forestploter)
library(topr)
library(grid)
library(forestploter)

##################
## get genes first
##################
#sig.snps <- c('rs2736100', 'rs2853677', 'rs7705526', 'rs228606','rs11212666', 'rs10890839', 'rs188761458')
#Gene <- c('TERT', 'TERT', 'TERT','NPAT', 'ATM', 'C11orf65',  'MSI2')
#nearest gene - correct the kessler bp location first
ch.snps.chk <- readRDS('ch_snpstoCheck.Rds')
ch.snps.chk$POS <- as.numeric(sub('.*:','',ch.snps.chk$snp.id))
ch.snps.chk <- annotate_with_nearest_gene(ch.snps.chk, protein_coding_only = T, build = 37)
# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS)]
head(hrc.rsid)
# add rsid to ch.snps.chk
ch.snps.chk <- merge(ch.snps.chk, hrc.rsid[,c(3,6)], by.y='rsid', by.x ='snp.id', all.x = T)
rm(hrc.rsid)
write.csv(ch.snps.chk, 'nearest_gene.csv', quote = F, row.names = F)


####################
# start from there
####################
ch_get_dir <- function() {
  kessler.full <- readRDS('Kessler_chip_full.Rds')
  kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
  kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
  kar.full <- readRDS('Kar_ch_full.Rds')
  kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
  kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
  AML.full <- readRDS('AML_combforCH_snplist.Rds')
  ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
  #AML.top <- readRDS('AML_all_topforCH_snplist.Rds')
  
  #add weinstock snp too
  weinstock.snp <- read.delim('sumstats/phs001974.pha005266.txt', comment.char = '#')
  names(weinstock.snp)
  weinstock.snp <- weinstock.snp[weinstock.snp$Rank <= 1,]
  weinstock.snp$study <- 'Weinstock et al 2023'
  weinstock.snp$trait <- 'CHIP-DNMT3A'
  weinstock.snp$case_con <- '3931/70277'
  weinstock.snp$snp.id <- '14:96180695'
  setnames(weinstock.snp, c('SNP.ID', 'P.value', 'Chr.ID', 'Chr.Position', 'Allele1', 'Allele2'),
           c('variant_id','p_value','chromosome','base_pair_location','other_allele','effect_allele'))
  # add OR
  weinstock.snp$odds_ratio <- exp(weinstock.snp$X..beta..)
  weinstock.snp$ci_lower <- exp(weinstock.snp$X..beta.. - 1.96 * weinstock.snp$SE)
  weinstock.snp$ci_upper <- exp(weinstock.snp$X..beta.. + 1.96 * weinstock.snp$SE)
  setDT(weinstock.snp)
  
  #add study
  kessler.full$study <- 'Kessler et al 2022'
  kar.full$study <- 'Kar et al 2022'
  #get snp.id for kar
  kar.full$snp.id <- paste0(kar.full$chromosome,':', kar.full$base_pair_location)
  table(ch.snps.chk$source)
  # divide into AML>ch dataset and CH>AML
  common.cols <- intersect(names(kessler.full), names(kar.full))
  common.cols <- intersect(common.cols, names(weinstock.snp))
  
  #from kar study, remove kes snps and wise versa
  kar.full <- kar.full[kar.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source %in% c('kar et al', 'Weinstock et al')]]
  kessler.full <- kessler.full[kessler.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source%in% c('kesslar et al', 'Weinstock et al')]]
  #combine
  
  kar.kess <- lapply(list(kessler.full, kar.full, weinstock.snp), function(x){
    x[,common.cols, with=F]
  })
  
  kar.kess <- do.call(rbind,kar.kess)
  sp <- unique(kar.kess[,c(11,13)])
  table(table(sp$variant_id))
  
  #prune kar.kess again
  common.snps <- intersect(AML.full[[2]]$snp.id, kar.kess$snp.id)
  grep('14:96180695', x = common.snps ,value = T)
  #prune
  kar.kess <- kar.kess[kar.kess$snp.id %in% common.snps]
  kar.kess$p_adj <- kar.kess$p_value
  
  #get common snps again for panaml and cnaml then adjust p
  table(AML.full[[2]]$snp.id %in% common.snps)#check
  names(AML.full[[1]])
  names(kar.kess)
  
  #get 95CIs for panAML and cnaml
  AML.full <- lapply(AML.full, function(x){
    x <- x[x$snp.id %in% common.snps]
    beta <- log(x[,OR])
    Z <- abs(qnorm(x[,P]/2, lower.tail = FALSE))
    SE <- abs(beta) / Z
    x$ci_lower <- exp(beta - 1.96 * SE)
    x$ci_upper <- exp(beta + 1.96 * SE)
    p.adj <- p.adjust(x[,P])
    x[,'p_adj'] <- p.adj
    setnames(x, c('CHR', 'A2', 'A1', 'OR', 'study', 'ID', 'P'),c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'trait', 'variant_id', 'p_value'))
    x$study <- 'AML GWAS'
    return(x)})
  
  AML.full[[1]]$case_con <- '4710/12938'
  AML.full[[2]]$case_con <- '1583/12938'
  #check
  common.col <- intersect(names(AML.full[[2]]), names(kar.kess))
  ch_aml.final <- lapply(list(AML.full[[1]], AML.full[[2]], kar.kess), function(x){
    x[,common.col, with=F]
  })
  
  ch_aml.final <- do.call(rbind, ch_aml.final)
  sp <- unique(ch_aml.final[,c(8,12)])
  table(table(sp$variant_id))
  length(unique(sp$variant_id))
  table(table(sp$variant_id))
  
  ####################
  ## harmonize alleles
  ####################
  # harmonize ORs based on EA
  # Define your reference effect allele per SNP
  #ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  ref <- ch_aml.final[trait ==  "Pan-AML", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  # Merge reference info to all rows
  ch_aml.final <- merge(ch_aml.final, ref, by = "variant_id", no.dups = T, all.x = T)
  # Flip odds ratios if study's effect allele != reference effect allele
  #correct this to only flip if EA=NEA
  ch_aml.final[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
  ch_aml.final[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_lower)]
  ch_aml.final[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_upper)]
  #remove NA OR_harm rows
  ch_aml.final <- ch_aml.final[!is.na(ch_aml.final$OR_harmonized)]
  
  #remove rs2887399 from comb data as this additional snp is both sig in Weinstock and kessler dnmt3a
  rs2887399_sub <- ch_aml.final[ch_aml.final$variant_id=='rs2887399']
  #ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id=='rs2887399'] # and remove
  
  ## add genes
  ch_aml.final <- merge(ch_aml.final, unique(ch.snps.chk[,c(1,5)]), by = 'snp.id', all.x = T)
  
  # order by rsid and p-value
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  table(ch_aml.final$trait)
  
  #sort snps based on trait first
  ch_aml.final$trait <- factor(ch_aml.final$trait, 
                               levels = c('Pan-AML', 'CN-AML', 'CH_inclusive', 'CH-DNMT3A', 'CH-TET2','CHIP_inclusive', 
                                          'CHIP-DNMT3A', 'CHIP-TET2', 'CHIP'))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$trait)]
  
  

  snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome) # use chr order
  ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  #ch_aml.final <- ch_aml.final %>% group_by(chromosome) %>% arrange(variant_id, .by_group = T)
  #drop 7 and 9 rows - multi-allelic from AML
  ch_aml.final <- ch_aml.final[-c(7,9)] # use this if not want multiallelic
  
  ################
  ### check direction
  #################
  head(ch_aml.final)
  # get gwas
  gwas.ch.aml <- ch_aml.final[ch_aml.final$trait %in% c('CHIP_inclusive', 'CH_inclusive')]
  gwas.ch.aml$dir <- ifelse(gwas.ch.aml$OR_harmonized >= 1, 1, 0)
  gwas.ch.aml.snp <- unique(gwas.ch.aml[,c(2,20)])
  table(gwas.ch.aml.snp$variant_id)
  #check against pan-aml
  ch_aml.final <- merge(ch_aml.final, gwas.ch.aml.snp, by='variant_id')
  ch_aml.final.gwas <- ch_aml.final[ch_aml.final$study =='AML GWAS']
  ch_aml.final.gwas$same.dir <- ifelse(ch_aml.final.gwas$OR_harmonized >= ch_aml.final.gwas$dir, 'same', 'opposite')
  table(unique(ch_aml.final.gwas[,c(1,21)])$same.dir)
  

}

ch_forest <- function(p.val) {
kessler.full <- readRDS('Kessler_chip_full.Rds')
kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
kar.full <- readRDS('Kar_ch_full.Rds')
kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
AML.full <- readRDS('AML_combforCH_snplist.Rds')
ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
#AML.top <- readRDS('AML_all_topforCH_snplist.Rds')

#add weinstock snp too
weinstock.snp <- read.delim('sumstats/phs001974.pha005266.txt', comment.char = '#')
names(weinstock.snp)
weinstock.snp <- weinstock.snp[weinstock.snp$Rank <= 1,]
weinstock.snp$study <- 'Weinstock et al 2023'
weinstock.snp$trait <- 'CHIP'
weinstock.snp$case_con <- '3931/70277'
weinstock.snp$snp.id <- '14:96180695'
setnames(weinstock.snp, c('SNP.ID', 'P.value', 'Chr.ID', 'Chr.Position', 'Allele1', 'Allele2'),
         c('variant_id','p_value','chromosome','base_pair_location','other_allele','effect_allele'))
# add OR
weinstock.snp$odds_ratio <- exp(weinstock.snp$X..beta..)
weinstock.snp$ci_lower <- exp(weinstock.snp$X..beta.. - 1.96 * weinstock.snp$SE)
weinstock.snp$ci_upper <- exp(weinstock.snp$X..beta.. + 1.96 * weinstock.snp$SE)
setDT(weinstock.snp)

#add study
kessler.full$study <- 'Kessler et al 2022'
kar.full$study <- 'Kar et al 2022'
#get snp.id for kar
kar.full$snp.id <- paste0(kar.full$chromosome,':', kar.full$base_pair_location)
table(ch.snps.chk$source)
# divide into AML>ch dataset and CH>AML
common.cols <- intersect(names(kessler.full), names(kar.full))
common.cols <- intersect(common.cols, names(weinstock.snp))

#from kar study, remove kes snps and wise versa
kar.full <- kar.full[kar.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source %in% c('kar et al', 'Weinstock et al')]]
kessler.full <- kessler.full[kessler.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source%in% c('kesslar et al', 'Weinstock et al')]]
#combine

kar.kess <- lapply(list(kessler.full, kar.full, weinstock.snp), function(x){
  x[,common.cols, with=F]
})

kar.kess <- do.call(rbind,kar.kess)
# add padj.col
kar.kess$p_adj <- ' '
kar.kess$p_adj2 <- kar.kess$p_value
sp <- unique(kar.kess[,c(11,13)])
table(table(sp$variant_id))

#prune kar.kess again
common.snps <- intersect(AML.full[[2]]$snp.id, kar.kess$snp.id)
grep('14:96180695', x = common.snps ,value = T)
#prune
kar.kess <- kar.kess[kar.kess$snp.id %in% common.snps]
#kar.kess$p_adj <- kar.kess$p_value
 
#get common snps again for panaml and cnaml then adjust p
table(AML.full[[2]]$snp.id %in% common.snps)#check
names(AML.full[[1]])
names(kar.kess)

#get 95CIs for panAML and cnaml
AML.full <- lapply(AML.full, function(x){
  x <- x[x$snp.id %in% common.snps]
  beta <- log(x[,OR])
  Z <- abs(qnorm(x[,P]/2, lower.tail = FALSE))
  SE <- abs(beta) / Z
  x$ci_lower <- exp(beta - 1.96 * SE)
  x$ci_upper <- exp(beta + 1.96 * SE)
  p.adj.std <- p.adjust(x[,P], 'bonferroni') 
  #p.adj <- p.adj.std * (51 / length(p.adj.std))
  #p.adj <- x[,P] * 51 / rank(x[,P])
  #p.adj <- x[,P] * 51
  #p.adj[p.adj > 1] <- 1
  x[,'p_adj'] <- formatC(p.adj.std, format = 'e', digits = 1)
  x[,'p_adj2'] <- p.adj.std
  setnames(x, c('CHR', 'A2', 'A1', 'OR', 'study', 'ID', 'P'),c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'trait', 'variant_id', 'p_value'))
  x$study <- 'AML GWAS'
  return(x)})

AML.full[[1]]$case_con <- '4710/12938'
AML.full[[2]]$case_con <- '1583/12938'
#check
common.col <- intersect(names(AML.full[[2]]), names(kar.kess))
ch_aml.final <- lapply(list(AML.full[[1]], AML.full[[2]], kar.kess), function(x){
  x[,common.col, with=F]
})

ch_aml.final <- do.call(rbind, ch_aml.final)
sp <- unique(ch_aml.final[,c(8,12)])
table(table(sp$variant_id))
length(unique(sp$variant_id))
table(table(sp$variant_id))

####################
## harmonize alleles
####################
# harmonize ORs based on EA
# Define your reference effect allele per SNP
#ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
ref <- ch_aml.final[trait ==  "Pan-AML", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
# Merge reference info to all rows
ch_aml.final <- merge(ch_aml.final, ref, by = "variant_id", no.dups = T, all.x = T)
# Flip odds ratios if study's effect allele != reference effect allele
#correct this to only flip if EA=NEA
ch_aml.final[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
ch_aml.final[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_upper)]
ch_aml.final[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_lower)]

#remove NA OR_harm rows
ch_aml.final <- ch_aml.final[!is.na(ch_aml.final$OR_harmonized)]

#remove rs2887399 from comb data as this additional snp is both sig in Weinstock and kessler dnmt3a
rs2887399_sub <- ch_aml.final[ch_aml.final$variant_id=='rs2887399']
#ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id=='rs2887399'] # and remove

## add genes
ch_aml.final <- merge(ch_aml.final, unique(ch.snps.chk[,c(1,5)]), by = 'snp.id', all.x = T)

# order by rsid and p-value
ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
table(ch_aml.final$trait)

#sort snps based on trait first
trait.order <- c('Pan-AML', 'CN-AML', 'CH_inclusive', 'CH-DNMT3A', 'CH-TET2','CHIP_inclusive', 
                                        'CHIP-DNMT3A', 'CHIP-TET2', 'CHIP')
ch_aml.final <- ch_aml.final %>% arrange(match(trait, trait.order))
snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome) # use chr order
ch_aml.final <- ch_aml.final[-c(37,97)] # use this if not want multiallelic of rs3219104
#ch_aml.final <- ch_aml.final %>% mutate(variant_id=factor(variant_id, levels = unique(snp.order$variant_id))) %>%  arrange(variant_id)

###############
##try plotting based on p cutoff
##############
names(ch_aml.final)
#plot each study separately to reduce size
p.cutoff <- p.val
snps.2.plt <- unique(ch_aml.final$variant_id[ch_aml.final$study=='AML GWAS' & ch_aml.final$p_adj2 <= p.cutoff])
test <- ch_aml.final[ch_aml.final$variant_id %in% snps.2.plt,]
snp.order <- test[test$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome)
test <- test %>% mutate(variant_id=factor(variant_id, levels = unique(snp.order$variant_id))) %>% arrange(variant_id)

test$` ` <- paste(rep(" ", 28), collapse = " ")
# Create a confidence interval column to display
test$`OR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", test$OR_harmonized, test$CI_low_harmonized, test$CI_upper_harmonized)
#test$p_adj <- format(test$p_adj, scientific=T, digits = 2)
test$p_value <- formatC(test$p_value, format = "e", digits = 1)
test$test.allele <- paste0(test$NEA_ref,'/',test$EA_ref)
names(test)

setnames(test, c('variant_id', 'chromosome', 'case_con', 'trait', 'study', 'test.allele', 'p_adj', 'p_value','Gene_Symbol'),
         c('SNP', 'CHR', 'Case/Con', 'Trait', 'Study', 'REF/EFF', 'P value \n(adjusted)','P value \n(unadjusted)', 'Gene'))

#add bg color
test <- test %>%  mutate(bg_group = as.numeric(as.factor(SNP)) %% 2)

bg_colors <- ifelse(test$bg_group == 0, "#F9FAF9", "#DEE4DD") #b6d1b9 #edfced

#add group labels
#test$snp.ea <- paste0(test$snp.id,'_',test$EA_ref)
test <- test %>% group_by(SNP) %>% mutate(`REF/EFF` = paste0("  ",`REF/EFF`))
test <- test %>% group_by(SNP) %>% mutate(CHR = ifelse(row_number() == 1, paste0("  ",CHR), ""))
test <- test %>% group_by(SNP) %>% mutate(`REF/EFF` = ifelse(row_number() == 1, paste0(" ",`REF/EFF`), "")) #remove this for type2 all plt
test <- test %>% group_by(SNP) %>% mutate(Gene = ifelse(row_number() == 1, Gene, ""))
test <- test %>% mutate(SNP = ifelse(duplicated(SNP), "", as.character(SNP)))

#separate AML ors and others
test <- test %>% mutate(OR_AML = ifelse(Study=='AML GWAS', OR_harmonized, NA))
test <- test %>% mutate(OR_other = ifelse(!Study=='AML GWAS', OR_harmonized, NA))
test <- test %>% mutate(CI_low_AML = ifelse(Study=='AML GWAS', CI_low_harmonized, NA))
test <- test %>% mutate(CI_low_other = ifelse(!Study=='AML GWAS', CI_low_harmonized, NA))
test <- test %>% mutate(CI_hi_AML = ifelse(Study=='AML GWAS', CI_upper_harmonized, NA))
test <- test %>% mutate(CI_hi_other = ifelse(!Study=='AML GWAS', CI_upper_harmonized, NA))

#change ch_invlusive to ch
test$Trait <- sub('_inclusive',' ',test$Trait)
#write.csv(test, 'plot_manuall_all.csv', quote = F, row.names = F)
# Define theme
tm <- forest_theme(base_size =  8,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = c(21),
                   ci_col = c("#519474", "#945171"),#9e0e13, #0d4012
                  # ci_fill = c("#9F608C", "#609F73"),
                   ci_alpha = 1,
                   ci_lty = 1,
                   ci_lwd = 1.5, 
                   legend_value = c('AML','CH'),
                   ci_Theight = 0, # Set a T end at the end of CI 
                   # Reference line width/type/color
                   refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey30"), # grey20
                   # Vertical line width/type/color
                   vertline_lwd = 1,legend_position = 'none',
                   vertline_lty = "dashed", 
                   # Change summary color for filling and borders
                   core = list(bg_params = list(fill=bg_colors)))

names(test)
p <- forest(test[,c(2,3,20,23,14,8,13,21,22,6,11)],
       est = list(test$OR_AML, test$OR_other),
       lower = list(test$CI_low_AML, test$CI_low_other), 
       upper = list(test$CI_hi_AML, test$CI_hi_other),
       ref_line = 1, sizes = 0.6, xlab = 'Odds Ratio',
       ci_column = 8, xlim = c(0.35,1.65), ticks_at = c(0.5,1,1.5),
       theme = tm)

p_wh <- get_wh(plot = p, unit = "in")
if (p.cutoff <= 0.05) {
pdf("forest_ch_aml_sig_new.pdf", width = p_wh[1], height = p_wh[2]) # 57 for all  # adjust dimensions as needed
grid.draw(p)
dev.off()
} else {
pdf("forest_ch_aml_all_V2.pdf", width = p_wh[1], height = p_wh[2]) # 57 for all  # adjust dimensions as needed
grid.draw(p)
dev.off()
}
}

ch_forest(0.05)

ch_forest_all <- function(snp.start, snp.end, prefix) {
  kessler.full <- readRDS('Kessler_chip_full.Rds')
  kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
  kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
  kar.full <- readRDS('Kar_ch_full.Rds')
  kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
  kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
  AML.full <- readRDS('AML_combforCH_snplist.Rds')
  ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
  #AML.top <- readRDS('AML_all_topforCH_snplist.Rds')
  
  #add weinstock snp too
  weinstock.snp <- read.delim('sumstats/phs001974.pha005266.txt', comment.char = '#')
  names(weinstock.snp)
  weinstock.snp <- weinstock.snp[weinstock.snp$Rank <= 1,]
  weinstock.snp$study <- 'Weinstock et al 2023'
  weinstock.snp$trait <- 'CHIP'
  weinstock.snp$case_con <- '3931/70277'
  weinstock.snp$snp.id <- '14:96180695'
  setnames(weinstock.snp, c('SNP.ID', 'P.value', 'Chr.ID', 'Chr.Position', 'Allele1', 'Allele2'),
           c('variant_id','p_value','chromosome','base_pair_location','other_allele','effect_allele'))
  # add OR
  weinstock.snp$odds_ratio <- exp(weinstock.snp$X..beta..)
  weinstock.snp$ci_lower <- exp(weinstock.snp$X..beta.. - 1.96 * weinstock.snp$SE)
  weinstock.snp$ci_upper <- exp(weinstock.snp$X..beta.. + 1.96 * weinstock.snp$SE)
  setDT(weinstock.snp)
  
  #add study
  kessler.full$study <- 'Kessler et al 2022'
  kar.full$study <- 'Kar et al 2022'
  #get snp.id for kar
  kar.full$snp.id <- paste0(kar.full$chromosome,':', kar.full$base_pair_location)
  table(ch.snps.chk$source)
  # divide into AML>ch dataset and CH>AML
  common.cols <- intersect(names(kessler.full), names(kar.full))
  common.cols <- intersect(common.cols, names(weinstock.snp))
  
  #from kar study, remove kes snps and wise versa
  kar.full <- kar.full[kar.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source %in% c('kar et al', 'Weinstock et al')]]
  kessler.full <- kessler.full[kessler.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source%in% c('kesslar et al', 'Weinstock et al')]]
  #combine
  
  kar.kess <- lapply(list(kessler.full, kar.full, weinstock.snp), function(x){
    x[,common.cols, with=F]
  })
  
  kar.kess <- do.call(rbind,kar.kess)
  # add padj.col
  kar.kess$p_adj <- ' '
  kar.kess$p_adj2 <- kar.kess$p_value
  sp <- unique(kar.kess[,c(11,13)])
  table(table(sp$variant_id))
  
  #prune kar.kess again
  common.snps <- intersect(AML.full[[2]]$snp.id, kar.kess$snp.id)
  grep('14:96180695', x = common.snps ,value = T)
  #prune
  kar.kess <- kar.kess[kar.kess$snp.id %in% common.snps]
  #kar.kess$p_adj <- kar.kess$p_value
  
  #get common snps again for panaml and cnaml then adjust p
  table(AML.full[[2]]$snp.id %in% common.snps)#check
  names(AML.full[[1]])
  names(kar.kess)
  
  #get 95CIs for panAML and cnaml
  AML.full <- lapply(AML.full, function(x){
    x <- x[x$snp.id %in% common.snps]
    beta <- log(x[,OR])
    Z <- abs(qnorm(x[,P]/2, lower.tail = FALSE))
    SE <- abs(beta) / Z
    x$ci_lower <- exp(beta - 1.96 * SE)
    x$ci_upper <- exp(beta + 1.96 * SE)
    p.adj.std <- p.adjust(x[,P], 'bonferroni') 
    #p.adj <- p.adj.std * (51 / length(p.adj.std))
    #p.adj <- x[,P] * 51 / rank(x[,P])
    #p.adj <- x[,P] * 51
    #p.adj[p.adj > 1] <- 1
    x[,'p_adj'] <- formatC(p.adj.std, format = "e", digits = 2)
    x[,'p_adj2'] <- p.adj.std
    setnames(x, c('CHR', 'A2', 'A1', 'OR', 'study', 'ID', 'P'),c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'trait', 'variant_id', 'p_value'))
    x$study <- 'AML GWAS'
    return(x)})
  
  AML.full[[1]]$case_con <- '4710/12938'
  AML.full[[2]]$case_con <- '1583/12938'
  #check
  common.col <- intersect(names(AML.full[[2]]), names(kar.kess))
  ch_aml.final <- lapply(list(AML.full[[1]], AML.full[[2]], kar.kess), function(x){
    x[,common.col, with=F]
  })
  
  ch_aml.final <- do.call(rbind, ch_aml.final)
  sp <- unique(ch_aml.final[,c(8,12)])
  table(table(sp$variant_id))
  length(unique(sp$variant_id))
  table(table(sp$variant_id))
  
  ####################
  ## harmonize alleles
  ####################
  # harmonize ORs based on EA
  # Define your reference effect allele per SNP
  #ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  ref <- ch_aml.final[trait ==  "Pan-AML", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  # Merge reference info to all rows
  ch_aml.final <- merge(ch_aml.final, ref, by = "variant_id", no.dups = T, all.x = T)
  # Flip odds ratios if study's effect allele != reference effect allele
  #correct this to only flip if EA=NEA
  ch_aml.final[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
  ch_aml.final[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_upper)]
  ch_aml.final[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_lower)]
  
  #remove NA OR_harm rows
  ch_aml.final <- ch_aml.final[!is.na(ch_aml.final$OR_harmonized)]
  
  #remove rs2887399 from comb data as this additional snp is both sig in Weinstock and kessler dnmt3a
  rs2887399_sub <- ch_aml.final[ch_aml.final$variant_id=='rs2887399']
  #ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id=='rs2887399'] # and remove
  
  ## add genes
  ch_aml.final <- merge(ch_aml.final, unique(ch.snps.chk[,c(1,5)]), by = 'snp.id', all.x = T)
  
  # order by rsid and p-value
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  table(ch_aml.final$trait)
  
  #sort snps based on trait first
  trait.order <- c('Pan-AML', 'CN-AML', 'CH_inclusive', 'CH-DNMT3A', 'CH-TET2','CHIP_inclusive', 
                   'CHIP-DNMT3A', 'CHIP-TET2', 'CHIP')
  ch_aml.final <- ch_aml.final %>% arrange(match(trait, trait.order))
  
  snos.2.exc <- unique(ch_aml.final$variant_id[ch_aml.final$study=='AML GWAS' & ch_aml.final$p_adj2 <= 0.05])
  ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id %in% snos.2.exc]
  snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome) # use chr order
  ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  #ch_aml.final <- ch_aml.final %>% group_by(chromosome) %>% arrange(variant_id, .by_group = T)
  #drop 7 and 9 rows - multi-allelic from AML
  ch_aml.final <- ch_aml.final[-c(7,9)] # use this if not want multiallelic
  #exclude sig snps in aml
  
  ###############
  ##try plotting based on p cutoff
  ##############
  names(ch_aml.final)
  #plot each study separately to reduce size
  snp.slice <- unique(snp.order$variant_id)[snp.start:snp.end]
  test <- ch_aml.final[ch_aml.final$variant_id %in% snp.slice]
  snps.2.plt <- unique(test$variant_id[test$study=='AML GWAS' & test$p_adj2 <= 1])
  test <- ch_aml.final[ch_aml.final$variant_id %in% snps.2.plt,]
  snp.order <- test[test$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome)
  test <- test %>% mutate(variant_id=factor(variant_id, levels = unique(snp.order$variant_id))) %>% arrange(variant_id)
  
  test$` ` <- paste(rep(" ", 24), collapse = " ")
  # Create a confidence interval column to display
  test$`OR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", test$OR_harmonized, test$CI_low_harmonized, test$CI_upper_harmonized)
  #test$p_adj <- format(test$p_adj, scientific=T, digits = 2)
  test$p_value <- formatC(test$p_value, format = "e", digits = 2)
  test$test.allele <- paste0(test$NEA_ref,'/',test$EA_ref)
  names(test)
  
  setnames(test, c('variant_id', 'chromosome', 'case_con', 'trait', 'study', 'test.allele', 'p_adj', 'p_value','Gene_Symbol'),
           c('SNP', 'CHR', 'Case/Con', 'Trait', 'Study', 'REF/EFF', 'P value \n(adjusted)','P value \n(unadjusted)', 'Gene'))
  
  #add bg color
  test <- test %>%  mutate(bg_group = as.numeric(as.factor(SNP)) %% 2)
  
  bg_colors <- ifelse(test$bg_group == 0, "#F9FAF9", "#DEE4DD") #b6d1b9 #edfced
  
  #add group labels
  #test$snp.ea <- paste0(test$snp.id,'_',test$EA_ref)
  test <- test %>% group_by(SNP) %>% mutate(`REF/EFF` = paste0("  ",`REF/EFF`))
  test <- test %>% group_by(SNP) %>% mutate(CHR = ifelse(row_number() == 1, paste0("  ",CHR), ""))
  test <- test %>% group_by(SNP) %>% mutate(`REF/EFF` = ifelse(row_number() == 1, paste0(" ",`REF/EFF`), "")) #remove this for type2 all plt
  test <- test %>% group_by(SNP) %>% mutate(Gene = ifelse(row_number() == 1, Gene, ""))
  test <- test %>% mutate(SNP = ifelse(duplicated(SNP), "", as.character(SNP)))
  
  #separate AML ors and others
  test <- test %>% mutate(OR_AML = ifelse(Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(OR_other = ifelse(!Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(CI_low_AML = ifelse(Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_low_other = ifelse(!Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_hi_AML = ifelse(Study=='AML GWAS', CI_upper_harmonized, NA))
  test <- test %>% mutate(CI_hi_other = ifelse(!Study=='AML GWAS', CI_upper_harmonized, NA))
  
  #change ch_invlusive to ch
  test$Trait <- sub('_inclusive',' ',test$Trait)
  #write.csv(test, 'plot_manuall_all.csv', quote = F, row.names = F)
  # Define theme
  tm <- forest_theme(base_size =  7,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = c(21),
                     ci_col = c("#519474", "#945171"),#9e0e13, #0d4012
                     # ci_fill = c("#9F608C", "#609F73"),
                     ci_alpha = 1,
                     ci_lty = 1,
                     ci_lwd = 1.5, 
                     legend_value = c('AML','CH'),
                     ci_Theight = 0, # Set a T end at the end of CI 
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey30"), # grey20
                     # Vertical line width/type/color
                     vertline_lwd = 1,legend_position = 'none',
                     vertline_lty = "dashed", 
                     # Change summary color for filling and borders
                     core = list(bg_params = list(fill=bg_colors)))
  
  names(test)
  p <- forest(test[,c(2,3,20,23,14,8,13,21,22,6,11)],
              est = list(test$OR_AML, test$OR_other),
              lower = list(test$CI_low_AML, test$CI_low_other), 
              upper = list(test$CI_hi_AML, test$CI_hi_other),
              ref_line = 1, sizes = 0.6, xlab = 'Odds Ratio',
              ci_column = 8, xlim = c(0.35,1.65), ticks_at = c(0.5,1,1.5),
              theme = tm)
  
  p_wh <- get_wh(plot = p, unit = "in")
  print(p_wh)
  
  pdf(paste0("ch_aml_split//forest_ch_aml_all_split",prefix,".pdf"), width = p_wh[1], height = p_wh[2]) #width = 8.3, height = 11.7
  grid.draw(p)
  dev.off()
  
  } # use this for split

ch_forest_all(43,56, '5')

##################################
## plot all forest in separate A4
#################################
# a4 size- 8.3 x 11.7in
ch_forest.A4 <- function(snp.start, snp.end, prefix) {
  kessler.full <- readRDS('Kessler_chip_full.Rds')
  kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
  kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
  kar.full <- readRDS('Kar_ch_full.Rds')
  kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
  kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
  AML.full <- readRDS('AML_combforCH_snplist.Rds')
  ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
  #AML.top <- readRDS('AML_all_topforCH_snplist.Rds')
  
  #add weinstock snp too
  weinstock.snp <- read.delim('sumstats/phs001974.pha005266.txt', comment.char = '#')
  names(weinstock.snp)
  weinstock.snp <- weinstock.snp[weinstock.snp$Rank <= 1,]
  weinstock.snp$study <- 'Weinstock et al 2023'
  weinstock.snp$trait <- 'CHIP'
  weinstock.snp$case_con <- '3931/70277'
  weinstock.snp$snp.id <- '14:96180695'
  setnames(weinstock.snp, c('SNP.ID', 'P.value', 'Chr.ID', 'Chr.Position', 'Allele1', 'Allele2'),
           c('variant_id','p_value','chromosome','base_pair_location','other_allele','effect_allele'))
  # add OR
  weinstock.snp$odds_ratio <- exp(weinstock.snp$X..beta..)
  weinstock.snp$ci_lower <- exp(weinstock.snp$X..beta.. - 1.96 * weinstock.snp$SE)
  weinstock.snp$ci_upper <- exp(weinstock.snp$X..beta.. + 1.96 * weinstock.snp$SE)
  setDT(weinstock.snp)
  
  #add study
  kessler.full$study <- 'Kessler et al 2022'
  kar.full$study <- 'Kar et al 2022'
  #get snp.id for kar
  kar.full$snp.id <- paste0(kar.full$chromosome,':', kar.full$base_pair_location)
  table(ch.snps.chk$source)
  # divide into AML>ch dataset and CH>AML
  common.cols <- intersect(names(kessler.full), names(kar.full))
  common.cols <- intersect(common.cols, names(weinstock.snp))
  
  #from kar study, remove kes snps and wise versa
  kar.full <- kar.full[kar.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source %in% c('kar et al', 'Weinstock et al')]]
  kessler.full <- kessler.full[kessler.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source%in% c('kesslar et al', 'Weinstock et al')]]
  #combine
  
  kar.kess <- lapply(list(kessler.full, kar.full, weinstock.snp), function(x){
    x[,common.cols, with=F]
  })
  
  kar.kess <- do.call(rbind,kar.kess)
  
  #prune kar.kess again
  common.snps <- intersect(AML.full[[2]]$snp.id, kar.kess$snp.id)
  grep('14:96180695', x = common.snps ,value = T)
  #prune
  kar.kess <- kar.kess[kar.kess$snp.id %in% common.snps]
  kar.kess$p_adj <- kar.kess$p_value
  
  #get common snps again for panaml and cnaml then adjust p
  table(AML.full[[2]]$snp.id %in% common.snps)#check
  names(AML.full[[1]])
  names(kar.kess)
  
  #get 95CIs for panAML and cnaml
  AML.full <- lapply(AML.full, function(x){
    x <- x[x$snp.id %in% common.snps]
    beta <- log(x[,OR])
    Z <- abs(qnorm(x[,P]/2, lower.tail = FALSE))
    SE <- abs(beta) / Z
    x$ci_lower <- exp(beta - 1.96 * SE)
    x$ci_upper <- exp(beta + 1.96 * SE)
    p.adj <- p.adjust(x[,P])
    x[,'p_adj'] <- p.adj
    setnames(x, c('CHR', 'A2', 'A1', 'OR', 'study', 'ID', 'P'),c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'trait', 'variant_id', 'p_value'))
    x$study <- 'AML GWAS'
    return(x)})
  
  AML.full[[1]]$case_con <- '4710/12938'
  AML.full[[2]]$case_con <- '1583/12938'
  #check
  common.col <- intersect(names(AML.full[[2]]), names(kar.kess))
  ch_aml.final <- lapply(list(AML.full[[1]], AML.full[[2]], kar.kess), function(x){
    x[,common.col, with=F]
  })
  
  ch_aml.final <- do.call(rbind, ch_aml.final)
  
  ####################
  ## harmonize alleles
  ####################
  # harmonize ORs based on EA
  # Define your reference effect allele per SNP
  #ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  ref <- ch_aml.final[trait ==  "Pan-AML", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  # Merge reference info to all rows
  ch_aml.final <- merge(ch_aml.final, ref, by = "variant_id", no.dups = T, all.x = T)
  # Flip odds ratios if study's effect allele != reference effect allele
  #correct this to only flip if EA=NEA
  ch_aml.final[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
  ch_aml.final[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_lower)]
  ch_aml.final[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_upper)]
  #remove NA OR_harm rows
  ch_aml.final <- ch_aml.final[!is.na(ch_aml.final$OR_harmonized)]
  
  #remove rs2887399 from comb data as this additional snp is both sig in Weinstock and kessler dnmt3a
  rs2887399_sub <- ch_aml.final[ch_aml.final$variant_id=='rs2887399']
  #ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id=='rs2887399'] # and remove
  
  ## add genes
  ch_aml.final <- merge(ch_aml.final, unique(ch.snps.chk[,c(1,5)]), by = 'snp.id', all.x = T)
  
  # order by rsid and p-value
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  table(ch_aml.final$trait)
  
  #sort snps based on trait first
  ch_aml.final$trait <- factor(ch_aml.final$trait, 
                               levels = c('Pan-AML', 'CN-AML', 'CH_inclusive', 'CH-DNMT3A', 'CH-TET2','CHIP_inclusive', 
                                          'CHIP-DNMT3A', 'CHIP-TET2', 'CHIP'))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$trait)]
  
  
  #create a snp list to sort based on high p - don't use
  # snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,11)] %>% arrange(p_adj)
  # ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  # ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome) # use chr order
  ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  #ch_aml.final <- ch_aml.final %>% group_by(chromosome) %>% arrange(variant_id, .by_group = T)
  #drop 7 and 9 rows - multi-allelic from AML
  ch_aml.final <- ch_aml.final[-c(7,9)] # use this if not want multiallelic
  
  ###############
  ##try plotting
  ##############
  names(ch_aml.final)
# slice for A4
  snp.slice <- unique(snp.order$variant_id)[snp.start:snp.end]
  test <- ch_aml.final[ch_aml.final$variant_id %in% snp.slice]
  snos.2.plt <- unique(test$variant_id[test$study=='AML GWAS' & test$p_adj <= 1])
  test <- test[test$variant_id %in% snos.2.plt,] # only take aml or ch sig hits
  test$` ` <- paste(rep(" ", 25), collapse = " ")
  # Create a confidence interval column to display
  test$`OR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", test$OR_harmonized, test$CI_low_harmonized, test$CI_upper_harmonized)
  test$p_adj <- formatC(test$p_adj, format = "e", digits = 2)
  test$p_value <- formatC(test$p_value, format = "e", digits = 2)
  test$test.allele <- paste0(test$NEA_ref,'/',test$EA_ref)
  names(test)
  
  setnames(test, c('variant_id', 'chromosome', 'case_con', 'trait', 'study', 'test.allele', 'p_adj', 'Gene_Symbol'),
           c('SNP', 'CHR', 'Case/Con', 'Trait', 'Study', 'OA/EA', 'P value', 'Gene'))
  
  #add bg color
  test <- test %>%
    mutate(bg_group = as.numeric(as.factor(SNP)) %% 2)
  
  bg_colors <- ifelse(test$bg_group == 0, "#F9FAF9", "#DEE4DD") #b6d1b9 #edfced
  
  #add group labels
  #test$snp.ea <- paste0(test$snp.id,'_',test$EA_ref)
  test <- test %>% group_by(SNP) %>% mutate(`OA/EA` = paste0("  ",`OA/EA`))
  test <- test %>% group_by(SNP) %>% mutate(CHR = ifelse(row_number() == 1, paste0("  ",CHR), ""))
  test <- test %>% group_by(SNP) %>% mutate(`OA/EA` = ifelse(row_number() == 1, paste0(" ",`OA/EA`), "")) #remove this for type2 all plt
  test <- test %>% group_by(SNP) %>% mutate(Gene = ifelse(row_number() == 1, Gene, ""))
  test <- test %>% mutate(SNP = ifelse(duplicated(SNP), "", as.character(SNP)))
  
  #separate AML ors and others
  test <- test %>% mutate(OR_AML = ifelse(Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(OR_other = ifelse(!Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(CI_low_AML = ifelse(Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_low_other = ifelse(!Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_hi_AML = ifelse(Study=='AML GWAS', CI_upper_harmonized, NA))
  test <- test %>% mutate(CI_hi_other = ifelse(!Study=='AML GWAS', CI_upper_harmonized, NA))
  
  #change ch_invlusive to ch
  test$Trait <- sub('_inclusive',' ',test$Trait)
  #write.csv(test, 'plot_manuall_all.csv', quote = F, row.names = F)
  # Define theme
  tm <- forest_theme(base_size =  7,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = c(21),
                     ci_col = c("#519474", "#945171"),#9e0e13, #0d4012
                     # ci_fill = c("#9F608C", "#609F73"),
                     ci_alpha = 1,
                     ci_lty = 1,
                     ci_lwd = 1.5, 
                     legend_value = c('AML','CH'),
                     ci_Theight = 0, # Set a T end at the end of CI 
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey30"), # grey20
                     # Vertical line width/type/color
                     vertline_lwd = 1,legend_position = 'none',
                     vertline_lty = "dashed", 
                     # Change summary color for filling and borders
                     core = list(bg_params = list(fill=bg_colors)))
  
  names(test)
  p <- forest(test[,c(2,3,19,22,13,8,12,20,21,11)],
              est = list(test$OR_AML, test$OR_other),
              lower = list(test$CI_low_AML, test$CI_low_other), 
              upper = list(test$CI_hi_AML, test$CI_hi_other),
              ref_line = 1, sizes = 0.5, xlab = 'Odds Ratio',
              ci_column = 8, xlim = c(0.35,1.65), ticks_at = c(0.5,1,1.5),
              theme = tm)
  
p_wh <- get_wh(plot = p, unit = "in")
print(p_wh)

pdf(paste0("forest_ch_aml_all_A4_",prefix,".pdf"), width = 8.3, height = 11.7) #width = 8.3, height = 11.7
grid.draw(p)
dev.off()

}

ch_forest.split <- function(snp.start, snp.end, prefix) {
  kessler.full <- readRDS('Kessler_chip_full.Rds')
  kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
  kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
  kar.full <- readRDS('Kar_ch_full.Rds')
  kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
  kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
  AML.full <- readRDS('AML_combforCH_snplist.Rds')
  ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
  #AML.top <- readRDS('AML_all_topforCH_snplist.Rds')
  
  #add weinstock snp too
  weinstock.snp <- read.delim('sumstats/phs001974.pha005266.txt', comment.char = '#')
  names(weinstock.snp)
  weinstock.snp <- weinstock.snp[weinstock.snp$Rank <= 1,]
  weinstock.snp$study <- 'Weinstock et al 2023'
  weinstock.snp$trait <- 'CHIP'
  weinstock.snp$case_con <- '3931/70277'
  weinstock.snp$snp.id <- '14:96180695'
  setnames(weinstock.snp, c('SNP.ID', 'P.value', 'Chr.ID', 'Chr.Position', 'Allele1', 'Allele2'),
           c('variant_id','p_value','chromosome','base_pair_location','other_allele','effect_allele'))
  # add OR
  weinstock.snp$odds_ratio <- exp(weinstock.snp$X..beta..)
  weinstock.snp$ci_lower <- exp(weinstock.snp$X..beta.. - 1.96 * weinstock.snp$SE)
  weinstock.snp$ci_upper <- exp(weinstock.snp$X..beta.. + 1.96 * weinstock.snp$SE)
  setDT(weinstock.snp)
  
  #add study
  kessler.full$study <- 'Kessler et al 2022'
  kar.full$study <- 'Kar et al 2022'
  #get snp.id for kar
  kar.full$snp.id <- paste0(kar.full$chromosome,':', kar.full$base_pair_location)
  table(ch.snps.chk$source)
  # divide into AML>ch dataset and CH>AML
  common.cols <- intersect(names(kessler.full), names(kar.full))
  common.cols <- intersect(common.cols, names(weinstock.snp))
  
  #from kar study, remove kes snps and wise versa
  kar.full <- kar.full[kar.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source %in% c('kar et al', 'Weinstock et al')]]
  kessler.full <- kessler.full[kessler.full$snp.id %in% ch.snps.chk$snp.id[ch.snps.chk$source%in% c('kesslar et al', 'Weinstock et al')]]
  #combine
  
  kar.kess <- lapply(list(kessler.full, kar.full, weinstock.snp), function(x){
    x[,common.cols, with=F]
  })
  
  kar.kess <- do.call(rbind,kar.kess)
  
  #prune kar.kess again
  common.snps <- intersect(AML.full[[2]]$snp.id, kar.kess$snp.id)
  grep('14:96180695', x = common.snps ,value = T)
  #prune
  kar.kess <- kar.kess[kar.kess$snp.id %in% common.snps]
  kar.kess$p_adj <- kar.kess$p_value
  
  #get common snps again for panaml and cnaml then adjust p
  table(AML.full[[2]]$snp.id %in% common.snps)#check
  names(AML.full[[1]])
  names(kar.kess)
  
  #get 95CIs for panAML and cnaml
  AML.full <- lapply(AML.full, function(x){
    x <- x[x$snp.id %in% common.snps]
    beta <- log(x[,OR])
    Z <- abs(qnorm(x[,P]/2, lower.tail = FALSE))
    SE <- abs(beta) / Z
    x$ci_lower <- exp(beta - 1.96 * SE)
    x$ci_upper <- exp(beta + 1.96 * SE)
    p.adj.std <- p.adjust(x[,P], 'bonferroni') 
    #p.adj <- p.adj.std * (51 / length(p.adj.std))
    #p.adj <- x[,P] * 51 / rank(x[,P])
    #p.adj <- x[,P] * 51
    #p.adj[p.adj > 1] <- 1
    x[,'p_adj'] <- p.adj.std
    setnames(x, c('CHR', 'A2', 'A1', 'OR', 'study', 'ID', 'P'),c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'trait', 'variant_id', 'p_value'))
    x$study <- 'AML GWAS'
    return(x)})
  
  AML.full[[1]]$case_con <- '4710/12938'
  AML.full[[2]]$case_con <- '1583/12938'
  #check
  common.col <- intersect(names(AML.full[[2]]), names(kar.kess))
  ch_aml.final <- lapply(list(AML.full[[1]], AML.full[[2]], kar.kess), function(x){
    x[,common.col, with=F]
  })
  
  ch_aml.final <- do.call(rbind, ch_aml.final)
  
  ####################
  ## harmonize alleles
  ####################
  # harmonize ORs based on EA
  # Define your reference effect allele per SNP
  #ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  ref <- ch_aml.final[trait ==  "Pan-AML", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
  # Merge reference info to all rows
  ch_aml.final <- merge(ch_aml.final, ref, by = "variant_id", no.dups = T, all.x = T)
  # Flip odds ratios if study's effect allele != reference effect allele
  #correct this to only flip if EA=NEA
  ch_aml.final[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
  ch_aml.final[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_upper)]
  ch_aml.final[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_lower)]
  #remove NA OR_harm rows
  ch_aml.final <- ch_aml.final[!is.na(ch_aml.final$OR_harmonized)]
  
  #remove rs2887399 from comb data as this additional snp is both sig in Weinstock and kessler dnmt3a
  rs2887399_sub <- ch_aml.final[ch_aml.final$variant_id=='rs2887399']
  #ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id=='rs2887399'] # and remove
  
  ## add genes
  ch_aml.final <- merge(ch_aml.final, unique(ch.snps.chk[,c(1,5)]), by = 'snp.id', all.x = T)
  
  # order by rsid and p-value
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  table(ch_aml.final$trait)
  
  #sort snps based on trait first
  ch_aml.final$trait <- factor(ch_aml.final$trait, 
                               levels = c('Pan-AML', 'CN-AML', 'CH_inclusive', 'CH-DNMT3A', 'CH-TET2','CHIP_inclusive', 
                                          'CHIP-DNMT3A', 'CHIP-TET2', 'CHIP'))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$trait)]
  
  
  #create a snp list to sort based on high p - don't use
  # snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,11)] %>% arrange(p_adj)
  # ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  # ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  snos.2.exc <- unique(ch_aml.final$variant_id[ch_aml.final$study=='AML GWAS' & ch_aml.final$p_adj <= 0.05])
  ch_aml.final <- ch_aml.final[!ch_aml.final$variant_id %in% snos.2.exc]
  snp.order <- ch_aml.final[ch_aml.final$study=='AML GWAS'][,c(2,3)] %>% arrange(chromosome) # use chr order
  ch_aml.final$variant_id <- factor(ch_aml.final$variant_id, levels = unique(snp.order$variant_id))
  ch_aml.final <- ch_aml.final[order(ch_aml.final$variant_id)]
  #ch_aml.final <- ch_aml.final %>% group_by(chromosome) %>% arrange(variant_id, .by_group = T)
  #drop 7 and 9 rows - multi-allelic from AML
  ch_aml.final <- ch_aml.final[-c(7,9)] # use this if not want multiallelic
  #exclude sig snps in aml

  length(unique(snp.order$variant_id))
  
  
  ###############
  ##try plotting
  ##############
  names(ch_aml.final)
  # slice for A4
  snp.slice <- unique(snp.order$variant_id)[snp.start:snp.end]
  test <- ch_aml.final[ch_aml.final$variant_id %in% snp.slice]
  snos.2.plt <- unique(test$variant_id[test$study=='AML GWAS' & test$p_adj <= 1])
  #test <- test[test$variant_id %in% snos.2.plt,] # only take aml or ch sig hits
  test$` ` <- paste(rep(" ", 25), collapse = " ")
  # Create a confidence interval column to display
  test$`OR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", test$OR_harmonized, test$CI_low_harmonized, test$CI_upper_harmonized)
  test$p_adj <- formatC(test$p_adj, format = "e", digits = 2)
  test$p_value <- formatC(test$p_value, format = "e", digits = 2)
  test$test.allele <- paste0(test$NEA_ref,'/',test$EA_ref)
  names(test)
  
  setnames(test, c('variant_id', 'chromosome', 'case_con', 'trait', 'study', 'test.allele', 'p_adj', 'Gene_Symbol'),
           c('SNP', 'CHR', 'Case/Con', 'Trait', 'Study', 'OA/EA', 'P value', 'Gene'))
  
  #add bg color
  test <- test %>%
    mutate(bg_group = as.numeric(as.factor(SNP)) %% 2)
  
  bg_colors <- ifelse(test$bg_group == 0, "#F9FAF9", "#DEE4DD") #b6d1b9 #edfced
  
  #add group labels
  #test$snp.ea <- paste0(test$snp.id,'_',test$EA_ref)
  test <- test %>% group_by(SNP) %>% mutate(`OA/EA` = paste0("  ",`OA/EA`))
  test <- test %>% group_by(SNP) %>% mutate(CHR = ifelse(row_number() == 1, paste0("  ",CHR), ""))
  test <- test %>% group_by(SNP) %>% mutate(`OA/EA` = ifelse(row_number() == 1, paste0(" ",`OA/EA`), "")) #remove this for type2 all plt
  test <- test %>% group_by(SNP) %>% mutate(Gene = ifelse(row_number() == 1, Gene, ""))
  test <- test %>% mutate(SNP = ifelse(duplicated(SNP), "", as.character(SNP)))
  
  #separate AML ors and others
  test <- test %>% mutate(OR_AML = ifelse(Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(OR_other = ifelse(!Study=='AML GWAS', OR_harmonized, NA))
  test <- test %>% mutate(CI_low_AML = ifelse(Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_low_other = ifelse(!Study=='AML GWAS', CI_low_harmonized, NA))
  test <- test %>% mutate(CI_hi_AML = ifelse(Study=='AML GWAS', CI_upper_harmonized, NA))
  test <- test %>% mutate(CI_hi_other = ifelse(!Study=='AML GWAS', CI_upper_harmonized, NA))
  
  #change ch_invlusive to ch
  test$Trait <- sub('_inclusive',' ',test$Trait)
  #write.csv(test, 'plot_manuall_all.csv', quote = F, row.names = F)
  # Define theme
  tm <- forest_theme(base_size =  7,
                     # Confidence interval point shape, line type/color/width
                     ci_pch = c(21),
                     ci_col = c("#519474", "#945171"),#9e0e13, #0d4012
                     # ci_fill = c("#9F608C", "#609F73"),
                     ci_alpha = 1,
                     ci_lty = 1,
                     ci_lwd = 1.5, 
                     legend_value = c('AML','CH'),
                     ci_Theight = 0, # Set a T end at the end of CI 
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = "dashed", col = "grey30"), # grey20
                     # Vertical line width/type/color
                     vertline_lwd = 1,legend_position = 'none',
                     vertline_lty = "dashed", 
                     # Change summary color for filling and borders
                     core = list(bg_params = list(fill=bg_colors)))
  
  names(test)
  p <- forest(test[,c(2,3,19,22,13,8,12,20,21,11)],
              est = list(test$OR_AML, test$OR_other),
              lower = list(test$CI_low_AML, test$CI_low_other), 
              upper = list(test$CI_hi_AML, test$CI_hi_other),
              ref_line = 1, sizes = 0.5, xlab = 'Odds Ratio',
              ci_column = 8, xlim = c(0.35,1.65), ticks_at = c(0.5,1,1.5),
              theme = tm)
  
  p_wh <- get_wh(plot = p, unit = "in")
  print(p_wh)
  
  pdf(paste0("ch_aml_split//forest_ch_aml_all_split",prefix,".pdf"), width = p_wh[1], height = p_wh[2]) #width = 8.3, height = 11.7
  grid.draw(p)
  dev.off()
  
}

ch_forest.split(46,55, '6')

#############################
## check LD within final snps
#############################
library(LDlinkR) # ld thresh 0.2
ch.snps.chk <- read.csv('nearest_gene.csv', header = T)
ch.snps.chk <- ch.snps.chk[!is.na(ch.snps.chk$ID),]
get.ld <- function(input, chr) {
  input <- input[!is.na(input$ID),]
  input <- input[input$CHROM==chr,]
  out <- LDmatrix(snps=input$ID, pop = 'EUR', token = 'fc368a631dc0')
  return(out) 
}

get.ld(ch.snps.chk, 11) # chr1,5, 11, 12, 17
# total snps for multiple correction - 51

########################
##delete - check snps #
#######################
pan.aml <- fread('../../AMLmeta_results/unfiltered/status_NCL_PCspeAMLHRC.meta', header = T) %>% mutate(study='Pan-AML')
head(pan.aml)
# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

pan.aml <- merge(pan.aml, hrc.rsid[,c(3,6)], by.x='SNP', by.y='rsid')
pan.aml[grep('rs4930561', pan.aml$ID)]
