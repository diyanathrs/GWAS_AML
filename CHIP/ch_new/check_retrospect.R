# load libs
library(tidyverse)
library(TwoSampleMR)
library(readxl)
library(data.table)

case.nos <- data_frame(trait= c('Pan-AML','del57-AML','complex-AML', 'CN-AML'), 
                                case_con= c('4710/12938','319/12938','358/12938', '1583/12938' ))
                       

# add N case no
#aml.sumstat$n <-case.nos[pheno]

plot_retro <- function() {
# start from there
kessler.full <- readRDS('Kessler_chip_full.Rds')
kessler.full$trait <- sub('DNMT3A','CHIP-DNMT3A',kessler.full$trait)
kessler.full$trait <- sub('TET2','CHIP-TET2',kessler.full$trait)
kar.full <- readRDS('Kar_ch_full.Rds')
kar.full$trait <- sub('TET2','CH-TET2',kar.full$trait)
kar.full$trait <- sub('DNMT3A','CH-DNMT3A',kar.full$trait)
#AML.full <- readRDS('AML_combforCH_snplist.Rds')
AML.top <- readRDS('AML_all_topforCH_snplist.Rds')
head(AML.top)

# finemapped snps
finemap.lead <- c('rs4665765', 'rs11481', 'rs3916765', 'rs79918355', 'rs12988876', 'rs12078864')
# correct direction of AML top hits
AML.top[, EA := ifelse(OR < 1, A2, A1)]
AML.top[, A2 := ifelse(OR < 1, A1, A2)]
AML.top[, OR := ifelse(OR < 1, 1 / OR, OR)]

#add study
kessler.full$study <- 'Kessler et al 2022'
kar.full$study <- 'Kar et al 2022'
AML.top$trait <- AML.top$study
AML.top$study <- 'AML GWAS'
# add case con nums
AML.top <- merge(AML.top, case.nos, by='trait')

#check against rsid
AML.top$ID %in% kar.full$variant_id  
AML.top$ID %in% kessler.full$variant_id 

#calculate CIs for AML
beta <- log(AML.top[,OR])
Z <- abs(qnorm(AML.top[,P]/2, lower.tail = FALSE))
SE <- abs(beta) / Z
AML.top$ci_lower <- exp(beta - 1.96 * SE)
AML.top$ci_upper <- exp(beta + 1.96 * SE)

setnames(AML.top, c('CHR', 'A2', 'EA', 'OR',  'ID', 'P'),
         c('chromosome', 'other_allele', 'effect_allele', 'odds_ratio', 'variant_id', 'p_value'))

common.cols <- intersect(names(AML.top), names(kar.full))
common.cols <- c(common.cols, 'p_adj')

AML.to.CH.all <- lapply(list(AML.top, kar.full, kessler.full), function(x){
  setDT(x)
  x <- x[x$variant_id %in% unique(AML.top$variant_id)]
  x <- x[x$variant_id %in% finemap.lead] #to check only finemap index
  # added new to multiple test correction of kar kess
  p.adj.std <- formatC(p.adjust(x[,p_value], 'bonferroni'),  format = "f", digits = 2)
  #p.adj.std <- format(p.adj.std, digits = 2)
  x[,'p_adj'] <- p.adj.std
  x[,common.cols, with=F]
})

#keep aml padj blank
AML.to.CH.all[[1]]$p_adj <- ' '

AML.to.CH.all <- do.call(rbind, AML.to.CH.all)
length(unique(AML.to.CH.all$variant_id))#check how many

####################
## harmonize alleles
####################
# harmonize ORs based on EA
# Define your reference effect allele per SNP
#ref <- ch_aml.final[trait %in%  c("CH_inclusive", "CHIP_inclusive"), .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
ref <- AML.to.CH.all[study ==  "AML GWAS", .(variant_id, EA_ref = effect_allele, NEA_ref = other_allele)]
# Merge reference info to all rows
AML.to.CH.all <- merge(AML.to.CH.all, ref, by = "variant_id", no.dups = T, all.x = T)
# Flip odds ratios if study's effect allele != reference effect allele
#correct this to only flip if EA=NEA
AML.to.CH.all[, OR_harmonized := ifelse(effect_allele == NEA_ref, 1 / odds_ratio, ifelse(effect_allele == EA_ref, odds_ratio, NA))]
AML.to.CH.all[, CI_low_harmonized := ifelse(effect_allele == EA_ref, ci_lower, 1 / ci_lower)]
AML.to.CH.all[, CI_upper_harmonized := ifelse(effect_allele == EA_ref, ci_upper, 1 / ci_upper)]
#remove NA OR_harm rows
AML.to.CH.all <- AML.to.CH.all[!is.na(AML.to.CH.all$OR_harmonized)]

## plot forest
# shorten p values
#AML.to.CH.all <- AML.to.CH.all[AML.to.CH.all$variant_id %in% finemap.lead,] # only take finemapped
AML.to.CH.all$` ` <- paste(rep(" ", 28), collapse = " ")
# Create a confidence interval column to display
AML.to.CH.all$`OR (95% CI)` <- sprintf("%.2f (%.2f to %.2f)", AML.to.CH.all$OR_harmonized, AML.to.CH.all$CI_low_harmonized,
                                       AML.to.CH.all$CI_upper_harmonized)
AML.to.CH.all$p_value <- formatC(AML.to.CH.all$p_value, format = "e", digits = 2)
AML.to.CH.all$test.allele <- paste0(AML.to.CH.all$NEA_ref,'/',AML.to.CH.all$EA_ref)
names(AML.to.CH.all)

setnames(AML.to.CH.all, c('variant_id', 'chromosome', 'case_con','EA_ref', 'trait', 'study', 'test.allele' ,'p_adj'),
         c('SNP', 'CHR', 'Case/Con','EA', 'Trait', 'Study', 'REF/EFF', 'P value \n(adjusted)')) 
#sort snps based on chr and add genes
AML.to.CH.all$SNP <- factor(AML.to.CH.all$SNP, levels = c('rs12078864', 'rs4665765','rs79918355', 'rs12988876', 'rs3916765', 'rs11481'))
unique(AML.to.CH.all$Trait)
# correct trait names
AML.to.CH.all$Trait <- sub('del57-AML', 'Del(5/7)-AML', AML.to.CH.all$Trait)
AML.to.CH.all$Trait <- sub('complex-AML', 'Complex-AML', AML.to.CH.all$Trait)

AML.to.CH.all$Trait <- factor(AML.to.CH.all$Trait, 
                             levels = c("Pan-AML", "CN-AML", "Complex-AML" , "Del(5/7)-AML", "CH_inclusive" ,"CH-DNMT3A",
                                        "CH-TET2", "CHIP_inclusive", "CHIP-DNMT3A", "CHIP-TET2"))
AML.to.CH.all <- AML.to.CH.all[order(AML.to.CH.all$Trait)]
AML.to.CH.all <- AML.to.CH.all[order(AML.to.CH.all$SNP)]
# add genes
genes <- data.frame(SNP= c('rs12078864', 'rs4665765','rs79918355', 'rs12988876', 'rs3916765', 'rs11481'), 
                    Gene=c('DUSP23', 'DNMT3A', 'EPCAM', 'PARD3B', 'HLA-DQA2', 'CHKA'))

#merge genes
head(genes)
AML.to.CH.all <- merge(AML.to.CH.all, genes, by = 'SNP')
#AML.to.CH.all <- AML.to.CH.all %>% arrange(SNP, Study, p_value, .by_group = F)

#add bg color
AML.to.CH.all <- AML.to.CH.all %>%
  mutate(bg_group = as.numeric(as.factor(SNP)) %% 2)

bg_colors <- ifelse(AML.to.CH.all$bg_group == 0, "#F9FAF9", "#DEE4DD") #b6d1b9

#add group labels
AML.to.CH.all <- AML.to.CH.all %>% group_by(SNP) %>% mutate(CHR = ifelse(row_number() == 1, paste0("  ",CHR), ""))
AML.to.CH.all <- AML.to.CH.all %>% group_by(SNP) %>% mutate(`REF/EFF` = ifelse(row_number() == 1, paste0("  ",`REF/EFF`), ""))
AML.to.CH.all <- AML.to.CH.all %>% group_by(SNP) %>% mutate(Gene = ifelse(row_number() == 1, Gene, ""))
AML.to.CH.all <- AML.to.CH.all %>% mutate(SNP = ifelse(duplicated(SNP), "", as.character(SNP)))

#separate AML ors and others
AML.to.CH.all <- AML.to.CH.all %>% mutate(OR_AML = ifelse(Study=='AML GWAS', OR_harmonized, NA))
AML.to.CH.all <- AML.to.CH.all %>% mutate(OR_other = ifelse(!Study=='AML GWAS', OR_harmonized, NA))
AML.to.CH.all <- AML.to.CH.all  %>% mutate(CI_low_AML = ifelse(Study=='AML GWAS', CI_low_harmonized, NA))
AML.to.CH.all <- AML.to.CH.all %>% mutate(CI_low_other = ifelse(!Study=='AML GWAS', CI_low_harmonized, NA))
AML.to.CH.all <- AML.to.CH.all  %>% mutate(CI_hi_AML = ifelse(Study=='AML GWAS', CI_upper_harmonized, NA))
AML.to.CH.all <- AML.to.CH.all %>% mutate(CI_hi_other = ifelse(!Study=='AML GWAS', CI_upper_harmonized, NA))

setnames(AML.to.CH.all, 'p_value', 'P value \n(unadjusted)')
AML.to.CH.all$Trait <- sub('_inclusive',' ',AML.to.CH.all$Trait)


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

names(AML.to.CH.all)
p <- forest(AML.to.CH.all[,c(1,3,21,20,9,2,7,18,19,5,12)],
            est = list(AML.to.CH.all$OR_AML, AML.to.CH.all$OR_other),
            lower = list(AML.to.CH.all$CI_low_AML, AML.to.CH.all$CI_low_other), 
            upper = list(AML.to.CH.all$CI_hi_AML, AML.to.CH.all$CI_hi_other),
            ref_line = 1, sizes = 0.6, xlab = 'Odds Ratio',
            ci_column = 8, xlim = c(0.5,3.5),ticks_at = c(0.5,1,2,3),
            theme = tm)

p_wh <- get_wh(plot = p, unit = "in")
pdf("forest_amlSNPSinCH.pdf", width = p_wh[1], height = p_wh[2]) # 57 for all  # adjust dimensions as needed
grid.draw(p)
dev.off()

}

plot_retro()
