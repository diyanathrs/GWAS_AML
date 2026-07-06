# cond analysis
#library(locuszoomr)
#library(EnsDb.Hsapiens.v75)
library(dplyr)
library(data.table)

# SNP ids should be same as in merged plink files (1:13380,1:13380:C:G )
# do chr-wise

## Run cojo for each sumstat file
# use unfil meta
sumstats <- list.files('../../AMLmeta_results/unfiltered', pattern = 'meta')
sumstats <- sumstats[c(4,5,8)]
sum.idx <- 3
pheno <- sub('_NCL_PCspeAMLHRC.meta','',sumstats[sum.idx])
pheno
sum.dat <- read.table(paste0('../../AMLmeta_results/unfiltered/',sumstats[sum.idx]),header = T) # change here for pheno
# for complex only use N=Nmax
max(sum.dat$N)
sum.dat <- subset(sum.dat, sum.dat$N==max(sum.dat$N))
head(sum.dat)
# check snp in sum.dat
# load ukb snps
#hrcfilein <- paste0("cut -f1-5 ~/AML/publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
#hrc.rsid <- fread(hrcfilein, header=TRUE)
#hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
#head(hrc.rsid)
#map snp_ids
#sum.dat <- merge(sum.dat, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
#sum.dat[which(sum.dat$ID=='rs115890122'),]


#### input format for GCTA CoJo ####
#SNP A1 A2 freq b se p N 
#rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
#rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
#rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
###################
# Should I do Cojo for assoc.res at GWAS level or meta file. assoc res has beta/SE etc. 
# cojo manual says its ok to use meta sumstats
# how to get freq of SNPs? get from cleaned rds?
# load gwas_attr files
#notes##
# get AF for snps, A1 is freq of effect allele but don't need exact values- use all_maf
# for ref file we can use merged bfile
# calculate beta and se from sumstats and add N manually
# calculate allele freq, take avg, only use status
attr.lst <- list.files('gwas_attr', pattern = paste0('NCL.*_status'), full.names = T,ignore.case = T)
stats.attr <- parallel::mclapply(attr.lst, readRDS, mc.cores = length(attr.lst))
head(stats.attr[[6]])
head(sum.dat)
# check all unfil snps are in attr files
sum.dat$SNP %in% stats.attr[[6]]$rsid

#only keep commmon snps from sumstats
for (i in 1:length(stats.attr)) {
  stats.attr[[i]]$SNP <- paste0(stats.attr[[i]]$chromosome,':',stats.attr[[i]]$position,'_',
                                stats.attr[[i]]$alleleA,'_',stats.attr[[i]]$alleleB)
  print(paste("No of SNPS pre filtering =", nrow(stats.attr[[i]])))
  stats.attr[[i]] <- stats.attr[[i]][stats.attr[[i]]$SNP %in% sum.dat$SNP]
  print(paste("No of SNPS post filtering =", nrow(stats.attr[[i]])))
}

head(stats.attr)
#206048000

# get maf - need to use maf only if EAF is the minor allele
for (i in 1: length(stats.attr)){
  stats.attr[[i]][all_BB>=all_AA, ':=' (all_maf=1-all_maf)]
}

# merge each maf col to cond
names(stats.attr[[1]])
head(sum.dat)
head(stats.attr[[6]][,c(19,10)])

for (i in 1:length(stats.attr)) {
sum.dat <- merge(sum.dat, stats.attr[[i]][,c(19,10)], by='SNP',all.x = T, no.dups = F) }
cond <- copy(sum.dat)
head(cond)
names(cond)

# calculate maf 
cond$maf <- rowSums(cond[c(13:18)])/6
head(cond)
rm(stats.attr)
gc()

# calculate beta se and add N
cond$b <- log(cond$OR)
cond$se <- abs(cond$b) / qnorm(1 - cond$P / 2)
cond$Ntot <- 4710
names(cond)
head(cond)
#change id to plink file style
cond$SNP <- paste0(cond$CHR,':',cond$BP,',',cond$CHR,':',cond$BP,':',cond$A2,':',cond$A1)
names(cond)
nrow(cond)

# check snp of interest - rs115890122 complex
cond[which(cond$ID=='rs4665765'),]

# save data
out.name <- pheno
write.table(cond[,c(1,4,5,19,20,21,7,22)], file = paste0(pheno,'_cond_1maf_in.ma'), quote = F, row.names = F, col.names = T)

stop()

## only for Normal
normal <- read.table('Normal_cond_in.ma', header = T)
head(normal)
normal.out <- normal[-c(grep('32682915', normal$SNP)),]
# remove the snp
write.table(normal.out, file = paste0(out.name,'_cond_in.ma'), quote = F, row.names = F, col.names = T)

# for del5_7
del57 <- read.table('del5_7_cond_in.ma', header = T)
subset(del57, del57$P <5e-6)
del57[grep('', del57$SNP),]

subset(cond, cond$P < 5e-8)
subset(sum.dat, sum.dat$P < 5e-8)
sum.dat[grep('206062950', sum.dat$BP),]

###
normal <- read.table('../../../AMLmeta_results/Normal_NCL_PCspeAMLHRC.meta.gz', header = T)
head(normal)
normal[grep('32671248', normal$BP),]

# complex check
cond[grep('2:206062950',cond$SNP),]
sum.dat[grep('1:159116916', sum.dat$SNP),]
sum.dat[grep('1:159125075', sum.dat$SNP),]

max(sum.dat$N)

sum.dat$SNP

# check snps
