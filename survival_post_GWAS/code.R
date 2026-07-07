# plot AMLsur surv plots
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(data.table)
library(purrr)
library(broom)

##############
## chr14 hit #
##############
snp <- c('rs6450183') #89695597
vcf.names <- read.delim('AMLsur_chr14_Allmerged.vcf', skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
chr14.vcf <- fread(paste('grep ',snp, 'AMLsur_chr14_Allmerged.vcf'),sep = '\t', header = F) 
colnames(chr14.vcf) <- names(vcf.names)
ref <- chr14.vcf$REF
alt <- chr14.vcf$ALT

snp.interest <- as.data.frame(t(chr14.vcf))
snp.interest <- snp.interest[-(1:9), ,drop= F]

snp.interest$GT <- substr(snp.interest$V1,1,3)
#remove no call
snp.interest <- snp.interest %>% filter(!GT =='./.')
table(snp.interest$GT)
#name <- paste0('chr14_',snp,'_G_C') 
#snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt))
snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref), 
                                ifelse(snp.interest$GT=='0/1', paste0(ref,'/',alt),paste0(alt,'/',alt)))
barplot(table(snp.interest$Genotype))

# creating pheno table
pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) #%>% filter(!.id =="BirminghamNCL3")
#remove APL
APL.cases <- which(pheno$t.15.17.==-9)
pheno.noAPL <- pheno[APL.cases,]
#read.delim('AMLsur_Birmingham_Dec2024.txt',colClasses = 'character')

table(pheno.noAPL$t.15.17., useNA = 'ifany')

snp.interest$sample <- row.names(snp.interest)
#add cohorts
# only keep snp.interest samples in the pheno file
snp.interest <- merge(pheno.noAPL[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
table(snp.interest$.id, useNA = 'ifany')

table(snp.interest$Genotype)
#dat <- snp.interest.out %>% filter(genotype_G_C %in% c('0/0', '0/1'))
#dat <- dat %>% filter(!.id=='BirminghamNCL3')
table(snp.interest$GT)
table(snp.interest$.id)

#snp.interest$Genotype <- factor(snp.interest$Genotype, levels=c('G/G','G/C or C/C'))
#plot
survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time in months',
             legend.title="Genotype", palette=c("red", "blue", "green"), 
             title=paste0("Kaplan-Meier curve for ",'chr14:',chr14.vcf$POS),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE, xscale = "y_m",xlim = c(0, 6),break.x.by=1)



survfit(Surv(OSdays/365.25, OSstatus) ~Genotype , data = snp.interest) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )

coxph(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest) %>% 
  tbl_regression(exp = TRUE) 

coxph(Surv(OSdays/365.25, OSstatus) ~.id, data = snp.interest) %>% 
  tbl_regression(exp = TRUE)


##############
## chr10 hit #
##############
snp <- c('132739739')
vcf.names <- read.delim('AMLsur_chr10_Allmerged.vcf', skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
chr.vcf <- fread(paste('grep ',snp, 'AMLsur_chr10_Allmerged.vcf'),sep = '\t', header = F) 
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
snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt))
barplot(table(snp.interest$Genotype))

# creating pheno table
pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) #%>% filter(!.id =="BirminghamNCL3")
#remove APL
APL.cases <- which(pheno$t.15.17.==-9)
pheno.noAPL <- pheno[APL.cases,]
#read.delim('AMLsur_Birmingham_Dec2024.txt',colClasses = 'character')

table(pheno.noAPL$t.15.17., useNA = 'ifany')

snp.interest$sample <- row.names(snp.interest)
#add cohorts
# only keep snp.interest samples in the pheno file
snp.interest <- merge(pheno.noAPL[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
table(snp.interest$.id, useNA = 'ifany')

table(snp.interest$Genotype)
#dat <- snp.interest.out %>% filter(genotype_G_C %in% c('0/0', '0/1'))
#dat <- dat %>% filter(!.id=='BirminghamNCL3')
table(snp.interest$GT)
table(snp.interest$.id)

snp.interest$Genotype <- factor(snp.interest$Genotype, levels=c(paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt)))
#plot
survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time in months',
             legend.title="Genotype", palette=c("red", "blue"), 
             title=paste0("Kaplan-Meier curve for ",'chr10:',chr.vcf$POS),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE, xscale = "y_m",xlim = c(0, 13),break.x.by=2)


survfit(Surv(OSdays/365.25, OSstatus) ~Genotype , data = snp.interest) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )

coxph(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest) %>% 
  tbl_regression(exp = TRUE) 

coxph(Surv(OSdays/365.25, OSstatus) ~.id, data = snp.interest) %>% 
  tbl_regression(exp = TRUE)


###################
## Check any hit ## use this for km plots
###################

snp.list <- read.table('snp_list_holly.txt', header = T)
#custom
rsid <- 'rs6450183'
chr <- 5
Position <- 49849445 
snp.list <- data.frame(rsid,chr,Position)


pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T)# %>% filter(!.id=='BirminghamNCL5 ')
#n <- 1
pdf(file = paste0('AMLsur_KMplots_',Sys.Date(),'.pdf'),width = 12, height = 7,onefile = T )
#pdf(file = 'AMLsur_plots.pdf',width = 12, height = 7, onefile = T )

for (n in 1:nrow(snp.list)) {
snp <- snp.list[n,]

vcf.names <- read.delim(paste0('vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
chr.vcf <- fread(paste0('grep ',snp$Position, ' vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
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
APL.cases <- which(pheno$t.15.17.==-9)
pheno.noAPL <- pheno[APL.cases,]
table(pheno.noAPL$t.15.17., useNA = 'ifany')
pheno.eln2 <- subset(pheno.noAPL, pheno.noAPL$ELN.22==2)
table(pheno.eln2$ELN.22, useNA = 'ifany')

snp.interest$sample <- row.names(snp.interest)
#add cohorts
# only keep snp.interest samples in the pheno file
snp.interest.all <- merge(pheno[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
snp.interest.eln2 <- merge(pheno.eln2[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
table(snp.interest.all$.id, useNA = 'ifany')
table(snp.interest.eln2$.id, useNA = 'ifany')

table(snp.interest.all$Genotype)
#dat <- snp.interest.out %>% filter(genotype_G_C %in% c('0/0', '0/1'))
#dat <- dat %>% filter(!.id=='BirminghamNCL3')
table(snp.interest.all$GT)

#snp.interest$Genotype <- factor(snp.interest$Genotype, levels=c(paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt)))
#plot

sur.all <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.all) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time in years',
             legend.title=paste(snp$rsid), palette=c("red", "blue", "forestgreen"), 
             title=paste0("All AML -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

sur.eln2 <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.eln2) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
             xlab='Time in years',
             legend.title=paste(snp$rsid), palette=c("red", "blue", "forestgreen"), 
             title=paste0("ELN2 only -",' chr:' ,snp$chr,snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

plot <- cowplot::plot_grid(sur.all$plot,sur.eln2$plot, sur.all$table, sur.eln2$table, rel_heights = c(6,2))

print(plot)
}
dev.off()


###############################
## Check any hit plot 2 only ##
###############################
snp.list <- read.table('snp_list_holly.txt', header = T)
#custom
rsid <- 'rs4665765'
chr <- 2
Position <- 25362515
snp.list <- data.frame(rsid,chr,Position)

pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) %>% filter(!.id=='BirminghamNCL3')
names(pheno)
pheno$Age.at.Diagnose <- as.numeric(pheno$Age.at.Diagnose)
#remove missing
pheno <- pheno[pheno$Age.at.Diagnose >= 0,] 
pheno <- pheno[!is.na(pheno$Age.at.Diagnose),]
range(pheno$Age.at.Diagnose)

#pheno$.id[pheno$.id %in% c('BirminghamNCL5', 'BirminghamNCL3')] <- 'UK3'

### close ##########
table(pheno$.id)
#n <- 1
#pdf(file = paste0('AMLsur_KMplots_',Sys.Date(),'.pdf'),width = 12, height = 7,onefile = T )
#pdf(file = 'AMLsur_plots.pdf',width = 12, height = 7, onefile = T )
n <- 1
snp <- snp.list[n,]

vcf.names <- read.delim(paste0('vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
chr.vcf <- fread(paste0('grep ',snp$Position, ' vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
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
APL.cases <- which(pheno$t.15.17.==-9)
pheno.noAPL <- pheno[APL.cases,]
table(pheno.noAPL$t.15.17., useNA = 'ifany')
pheno.eln2 <- subset(pheno.noAPL, pheno.noAPL$ELN.22==2)
table(pheno.eln2$ELN.22, useNA = 'ifany')
  
snp.interest$sample <- row.names(snp.interest)
  #add cohorts
  # only keep snp.interest samples in the pheno file
names(pheno)
snp.interest.all <- merge(pheno[c(1:4, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
snp.interest.eln2 <- merge(pheno.eln2[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
table(snp.interest.all$.id, useNA = 'ifany')
table(snp.interest.eln2$.id, useNA = 'ifany')
  
names(snp.interest.all)
#make groups
  
#snp.interest$Genotype <- factor(snp.interest$Genotype, levels=c(paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt)))
####open-plots#####

# surv on study - age groups 
snp.interest.all %>%  filter(age_group=='>median') %>% survfit(Surv(OSdays/365.25, OSstatus) ~.id, data = .) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.int.homo,
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)
  

snp.int.homo %>% group_by(.id) %>% do(tidy(coxph(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)))
snp.int.homo %>% group_by(.id) %>% do(tidy(survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)))

########Age plts####
# KP for sex & age groups#
####################
# check best age cutoff
age.cut <- surv_cutpoint(snp.interest.all, time = "OSdays", event = "OSstatus", variables = "Age.at.Diagnose")
df_cat <- surv_categorize(age.cut)
fit <- survfit(Surv(OSdays, OSstatus) ~ Age.at.Diagnose, data = df_cat)

ggsurvplot(fit, data = df_cat, pval = TRUE, risk.table = TRUE)

# age
median(snp.interest.all$Age.at.Diagnose)
ggplot(snp.interest.all, aes(Age.at.Diagnose)) + geom_density(fill='steelblue', alpha=0.5)+
  geom_vline(aes(xintercept = median(Age.at.Diagnose)), color = "red", linetype = "dashed")

#density by event status
ggplot(snp.interest.all, aes(x = Age.at.Diagnose, fill = factor(OSstatus))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("gray70", "red"), labels = c("Censored", "Event")) +
  geom_vline(aes(xintercept = median(Age.at.Diagnose)), color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Age Distribution by Event Status",
       x = "Age", fill = "Status")

#den for study
ggplot(snp.interest.all, aes(x = Age.at.Diagnose, fill = .id)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Age Distribution by Group", x = "Age", fill = "Group")+
  geom_vline(aes(xintercept = median(Age.at.Diagnose)), color = "red", linetype = "dashed") 

#stratify by age
pheno <- snp.interest.all %>% mutate(age_group = ifelse(Age.at.Diagnose < median(snp.interest.all$Age.at.Diagnose), "<median", ">median"))

# KP by age all
pheno %>%  survfit(Surv(OSdays/365.25, OSstatus) ~age_group, data = .) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, pheno, 
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title="Survival Stratified by Age Group",
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# KP by age groups
snp.interest.all %>%  survfit(Surv(OSdays/365.25, OSstatus) ~age_group, data = .) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, snp.interest.all, facet.by = '.id',
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title="Survival Stratified by Age Group",
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

ggsave('rs4665765/surv_age.pdf', height = 8, width = 10, units = 'in')

#sex
snp.interest.all %>%  survfit(Surv(OSdays/365.25, OSstatus) ~Gender, data = .) %>% 
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, snp.interest.all, facet.by = '.id',
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title="Survival Stratified by Gender",
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

ggsave('rs4665765/surv_gender.pdf', height = 8, width = 10, units = 'in')

table(snp.interest.all$Gender)
snp.interest.all  %>% survfit(Surv(OSdays/365.25, OSstatus) ~Gender, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.interest.all,
             xlab='Time in years', 
             title=paste0("Survival Stratified by Gender "),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

#snp by gender - full
sex <- 'Male'
male <- snp.interest.all %>% filter(Gender==sex) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex),
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("Survival for ", sex,'s'),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

sex <- 'Female'
female <- snp.interest.all %>% filter(Gender==sex) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex),
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("Survival for ", sex,'s'),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

plot_grid(male$plot, female$plot, male$table, female$table, ncol= 2, rel_heights = c(3,1))
ggsave('rs4665765/gender_surv.pdf', height = 7, width = 12, units = 'in')

#snp by gender - cohorts
sex <- 'Female'
snp.interest.all %>% filter(Gender==sex) %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex), facet.by = ".id",
             xlab='Time in years', legend.title=paste(snp$rsid), title=paste0("Survival for ", sex,'s'), 
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)
ggsave('rs4665765/cohort_female.pdf', height = 8.5, width = 11, units = 'in')

# snp by gender & age - full plt
sex <- 'Male'
age <- '>median'
male.old <- snp.interest.all %>% filter(Gender==sex & age_group==age) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex & 
                                                                     snp.interest.all$age_group==age), 
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("Survival for ", sex,'s of age ',age, '(',median(snp.interest.all$Age.at.Diagnose),'y)' ),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

male <- plot_grid(male.young$plot, male.old$plot, male.young$table, male.old$table, ncol= 2, rel_heights = c(3,1))
female <- plot_grid(female.young$plot, female.old$plot, female.young$table, female.old$table, ncol= 2, rel_heights = c(3,1))
plot_grid(male, female, ncol = 1)

ggsave('rs4665765/gender_age.pdf', height = 12, width = 12, units = 'in')

# facet plots
sex <- "Female"
snp.interest.all %>% filter(Gender==sex) %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex), 
             facet.by = 'age_group',
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("Survival of ",sex,"s - age group ", '(median ',median(snp.interest.all$Age.at.Diagnose),')' ),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

ggsave('rs4665765/age_surv.pdf', height = 6, width = 10, units = 'in')

# snp by gender & age - cohorts
sex <- 'Male'
age <- '>median'
snp.interest.all %>% filter(Gender==sex & age_group==age) %>% group_by('.id') %>% survfit(Surv(OSdays/365.25, OSstatus) 
                                                                                          ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$Gender==sex & 
                                                                     snp.interest.all$age_group==age), 
             facet.by = ".id",
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste("Survival for", sex, age,'age group'),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

ggsave(paste0('rs4665765/surv',sex,'&',age,'.pdf'), height = 8, width = 10, units = 'in')

snp.interest.all %>% group_by(age_group) %>% survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)

# all groups plt
snp.interest.all %>% group_by(age_group) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.interest.all,facet.by = "age_group",
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title="Survival Stratified by Age Group",
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# groups sep
grp <- '<median'
snp.interest.all %>% filter(age_group==grp) %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$age_group==grp), facet.by = ".id",
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title=paste0("Survival Stratified by Age Group ",grp),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# full plt all
grp <- '>median'
one <- snp.interest.all %>% filter(age_group==grp) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$age_group==grp),
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title=paste0("Survival Stratified by Age Group ",grp),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

grp <- '<median'
two <- snp.interest.all %>% filter(age_group==grp) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$age_group==grp),
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title=paste0("Survival Stratified by Age Group ",grp),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

plot_grid(two$plot, one$plot, two$table, one$table, ncol= 2, rel_heights = c(3,1))

#full plt group
grp <- '<median'
stdy <- 'Germany'
snp.interest.all %>% filter(age_group==grp & .id==stdy) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.interest.all, snp.interest.all$age_group==grp 
                                                                   & snp.interest.all$.id==stdy),
             xlab='Time in years', legend.title=paste("Median Age=",median(pheno$Age.at.Diagnose)), 
             title=paste0("Survival Stratified by Age Group ",grp, " -",stdy),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)


############
# KP for hom
############
# no study
survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.all)

snp.interest.all %>% group_by(age_group) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.interest.all, facet.by = "age_group",
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

hom.all <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.homo) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE,
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)
hom.all

# separate studies
snp.int.homo %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.int.homo, facet.by = '.id',
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# separate studies
snp.interest.all %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.interest.all, facet.by = '.id',
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

snp.interest.all %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.interest.all,
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Homozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# check one study
std <- 'US1'
snp.int.homo %>% filter(.id==std) %>% survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)
snp.int.homo %>% filter(.id=='US1') %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.int.homo, snp.int.homo$.id==std),
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Hom-",std,' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

#############
# KP for het
############
survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.het)
het.all <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.het) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE,
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Heterozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)
het.all

# separate studies
snp.int.het %>% group_by(.id) %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = snp.int.het, facet.by = '.id',
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Heterozygous -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)

# check one study
std <- 'Finland'
snp.int.het %>% filter(.id==std) %>% survdiff(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)
snp.int.het %>% filter(.id=='US1') %>% survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, data = subset(snp.int.het, snp.int.het$.id==std),
             xlab='Time in years', legend.title=paste(snp$rsid), 
             title=paste0("NoAPL Het-",std,' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
             risk.table.height=.2, risk.table.y.text.col = T,
             risk.table.y.text = FALSE)



# print hom/het all
cowplot::plot_grid(hom.all$plot, het.all$plot, hom.all$table, het.all$table,rel_heights = c(6,2) )

# group and nest data
hom_grouped <- snp.int.homo %>% group_by(.id) %>% nest()
hom_km <- hom_grouped %>% mutate(km_fit=map(data, ~ survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.homo)))

hom_km %>% mutate(plot = map2(km_fit, group, ~ ggsurvplot(.x, data = snp.int.homo[snp.int.homo$.id == .y, ],
                                            title = paste("Group", .y),
                                            risk.table = TRUE)))


  
survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.het) %>%
  ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE,
               xlab='Time in years', facet.by = '.id',
               legend.title=paste(snp$rsid), palette=c("red", "blue"), 
               title=paste0("NoAPL AML -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
               risk.table.height=.2, risk.table.y.text.col = T,
               risk.table.y.text = FALSE)
  
  
  survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.int.homo) %>%
    ggsurvplot(data = snp.int.homo, conf.int=F, pval=TRUE, risk.table=TRUE,
               xlab='Time in years',  facet.by = '.id',
               legend.title=paste(snp$rsid), palette=c("red", "blue"), 
               title=paste0("NoAPL AML -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
               risk.table.height=.2, risk.table.y.text.col = T,
               risk.table.y.text = FALSE)

  

  
  plot <- cowplot::plot_grid(sur.all$plot,sur.eln2$plot, sur.all$table, sur.eln2$table, rel_heights = c(6,2))
  
print(plot)
  
#cox model
snp.int.het %>% group_by(.id) %>% do(tidy(coxph(Surv(OSdays/365.25, OSstatus) ~Genotype, data = .)))

########################
## GT calls for holly ##
########################
library(parallel)

#snp.list <- read.table('snp_list_holly.txt', header = T)
#add new SNPs
rsid <- c('rs6450183', 'rs16829165')
chr <- c(5, 3)
Position <- c('49849445', '118663904')
snp.list <- data.frame(rsid,chr,Position)

pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) 

infer_GT <- function(n) {
  snp <- snp.list[n,]
  vcf.names <- read.delim(paste0('vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
  chr.vcf <- fread(paste0('grep ',snp$Position, ' vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
  colnames(chr.vcf) <- names(vcf.names)
  ref <- chr.vcf$REF
  alt <- chr.vcf$ALT
  
  snp.interest <- as.data.frame(t(chr.vcf))
  snp.interest <- snp.interest[-(1:9), ,drop= F]
  #snp.interest$sample_id <- row.names(snp.interest)
  
  snp.interest$GT <- substr(snp.interest$V1,1,3)
  #remove no call
  #snp.interest <- snp.interest %>% filter(!GT =='./.')
  table(snp.interest$GT)
  #snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt))
  #snp.interest$Genotype <- 'No_call'
  #snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref), 
   #                               ifelse(snp.interest$GT=='0/1', paste0(ref,'/',alt),paste0(alt,'/',alt)))
  
  snp.interest$Genotype <- ifelse(snp.interest$GT=='0/0',paste0(ref,'/',ref), 
                                  ifelse(snp.interest$GT=='0/1', paste0(ref,'/',alt),
                                         ifelse(snp.interest$GT=='1/1',paste0(alt,'/',alt), paste0('No_call') )))
  #snp.interest <- snp.interest[c(1,3)]
  names(snp.interest)[3] <- paste0(snp$rsid)
  names(snp.interest)[1] <- paste0(snp$rsid,'_dosage')
  return(snp.interest[c(1,3)])
}

# run inferGT
test <- parallel::mclapply(seq_len(nrow(snp.list)), infer_GT, mc.cores = 2)
out <- as.data.frame(do.call(cbind,test))
out$sample_id <- row.names(out)
dt <- Sys.Date()
write.table(out, paste0('AMLsur_GTforHolly_',dt,'.txt'), quote = F, row.names = F, col.names = T, sep = '\t')

# plot KM
plot_km <- function(snp.list) {
pheno <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T)# %>% filter(!.id=='BirminghamNCL5 ')
#n <- 1
pdf(file = paste0('AMLsur_KMplots_',Sys.Date(),'.pdf'),width = 12, height = 7,onefile = T )
#pdf(file = 'AMLsur_plots.pdf',width = 12, height = 7, onefile = T )

for (n in 1:nrow(snp.list)) {
  snp <- snp.list[n,]
  
  vcf.names <- read.delim(paste0('vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
  chr.vcf <- fread(paste0('grep ',snp$Position, ' vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
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
  APL.cases <- which(pheno$t.15.17.==-9)
  pheno.noAPL <- pheno[APL.cases,]
  table(pheno.noAPL$t.15.17., useNA = 'ifany')
  pheno.eln2 <- subset(pheno.noAPL, pheno.noAPL$ELN.22==2)
  table(pheno.eln2$ELN.22, useNA = 'ifany')
  
  snp.interest$sample <- row.names(snp.interest)
  #add cohorts
  # only keep snp.interest samples in the pheno file
  snp.interest.all <- merge(pheno[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
  snp.interest.eln2 <- merge(pheno.eln2[c(1:2, 6:7)], snp.interest, by.x = 'IID', by.y = 'sample')
  table(snp.interest.all$.id, useNA = 'ifany')
  table(snp.interest.eln2$.id, useNA = 'ifany')
  
  table(snp.interest.all$Genotype)
  #dat <- snp.interest.out %>% filter(genotype_G_C %in% c('0/0', '0/1'))
  #dat <- dat %>% filter(!.id=='BirminghamNCL3')
  table(snp.interest.all$GT)
  
  #snp.interest$Genotype <- factor(snp.interest$Genotype, levels=c(paste0(ref,'/',ref),paste0(ref,'/',alt,' or ',alt,'/',alt)))
  #plot
  
  sur.all <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.all) %>% 
    ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
               xlab='Time in years',
               legend.title=paste(snp$rsid), palette=c("red", "blue", "forestgreen"), 
               title=paste0("All AML -",' chr' ,snp$chr,':',snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
               risk.table.height=.2, risk.table.y.text.col = T,
               risk.table.y.text = FALSE)
  
  sur.eln2 <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.eln2) %>% 
    ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
               xlab='Time in years',
               legend.title=paste(snp$rsid), palette=c("red", "blue", "forestgreen"), 
               title=paste0("ELN2 only -",' chr:' ,snp$chr,snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
               risk.table.height=.2, risk.table.y.text.col = T,
               risk.table.y.text = FALSE)
  
  plot <- cowplot::plot_grid(sur.all$plot,sur.eln2$plot, sur.all$table, sur.eln2$table, rel_heights = c(6,2))
  
  print(plot)
}
dev.off()
}
plot_km(snp.list)

### other ####
survfit(Surv(OSdays/365.25, OSstatus) ~Genotype , data = snp.interest) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )

coxph(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest) %>% 
  tbl_regression(exp = TRUE) 

coxph(Surv(OSdays/365.25, OSstatus) ~.id, data = snp.interest) %>% 
  tbl_regression(exp = TRUE)
