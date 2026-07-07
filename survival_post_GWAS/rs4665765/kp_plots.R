# kp plots for snps of int
rsid <- 'rs4665765'
chr <- 2
Position <- 25362515
snp.list <- data.frame(rsid,chr,Position)

pheno <- read.delim('../../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) %>% filter(!.id=='BirminghamNCL3')
names(pheno)
pheno$Age.at.Diagnose <- as.numeric(pheno$Age.at.Diagnose)
#remove missing
pheno <- pheno[pheno$Age.at.Diagnose >= 0,] 
pheno <- pheno[!is.na(pheno$Age.at.Diagnose),]
range(pheno$Age.at.Diagnose)

table(pheno$.id)

vcf.names <- read.delim(paste0('../vcf/AMLsur_chr',snp.list$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
chr.vcf <- fread(paste0('grep ',snp.list$Position, ' ../vcf/AMLsur_chr',snp.list$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
colnames(chr.vcf) <- names(vcf.names)
ref <- chr.vcf$REF
alt <- chr.vcf$ALT

snp.interest <- as.data.frame(t(chr.vcf))
snp.interest <- snp.interest[-(1:9), ,drop= F]
snp.interest$GT <- substr(snp.interest$V1,1,3)
#remove no call
snp.interest <- snp.interest %>% filter(!GT =='./.')
table(snp.interest$GT)
#name <- paste0('chr14_',snp,'_G_C') 
#snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt))
snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref), 
                                ifelse(snp.interest$GT=='0/1', paste0(ref,'/',alt),paste0(alt,'/',alt)))
#barplot(table(snp.interest$Genotype))

# creating pheno table
#remove APL
apl.idx <- which(pheno$t.15.17.==2)
pheno.noAPL <- pheno[-apl.idx,]
table(pheno$t.15.17., useNA = 'ifany')
table(pheno.noAPL$t.15.17., useNA = 'ifany')

snp.interest$sample <- row.names(snp.interest)
#add cohorts
# only keep snp.interest samples in the pheno file
names(pheno)
snp.interest.noapl <- merge(pheno.noAPL[c(1:4, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
snp.interest.all <- merge(pheno[c(1:4, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
table(snp.interest.all$.id, useNA = 'ifany')

names(snp.interest.all)
#make groups
####open-plots#####

# surv on study - age groups 

p <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.noapl) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time',  legend=c(0.8,0.8),
             legend.title=paste(snp.list$rsid), palette=c("#00BFC4", "#7CAE00", "#F8766D"), 
             risk.table.y.text = T )

cowplot::plot_grid(p$plot, p$table, ncol = 1, rel_heights = c(2,1))
#cowplot::plot_grid(p$plot)
ggsave(file = paste0('AMLsur_KMplots2_',Sys.Date(),'.pdf'), width = 7, height = 5, units = 'in')

#######################
# stratified KP plots
######################
median(snp.interest.noapl$Age.at.Diagnose)
#stratify by age
snp.interest.all <- snp.interest.all %>% mutate(age_group = ifelse(Age.at.Diagnose < median(snp.interest.all$Age.at.Diagnose), "<median", ">median"))
snp.interest.noapl <- snp.interest.noapl %>% mutate(age_group = ifelse(Age.at.Diagnose < median(snp.interest.noapl$Age.at.Diagnose), "<median", ">median"))

# first centre wise
p2 <- snp.interest.noapl %>%  group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=F, data = snp.interest.noapl, facet.by = ".id",
             xlab='Time in years', legend.title='Genotype', pval.coord = c(10, 0.8), 
             palette=c("#66CD00", "#B22222", "#4682B4"),
             title=paste0("Survival stratified by centre for ",snp.list$rsid))

cowplot::plot_grid(p2)
ggsave(file = paste0('AMLsur_KMplots_centre_',Sys.Date(),'.pdf'), width = 8, height = 7, units = 'in')

snp.interest.all %>%  group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=F, data = snp.interest.all, facet.by = ".id",
             xlab='Time in years', legend.title='Genotype', 
             title=paste0("Survival stratified by centre "))


grp <- '>median'
snp.interest.all %>% filter(age_group==grp) %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$age_group==grp), facet.by = ".id",
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title=paste0("Survival Stratified by Age Group ",grp),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)
