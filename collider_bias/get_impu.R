library(dplyr)
library(data.table)

# get average impu score for R2 filtering of ld prune
# get AF for snps, A1 is freq of effect allele but don't need exact values- use all_maf
# for ref file we can use merged bfile
# calculate allele freq, take avg, only use status
attr.lst <- list.files('~/AML/AML_condAnalysis/cojo/gwas_attr/', pattern = paste0('NCL.*_PanAML'), full.names = T, ignore.case = T)
stats.attr <- parallel::mclapply(attr.lst, readRDS, mc.cores = length(attr.lst))
stats.attr <- stats.attr[-6]
head(stats.attr)

# check all unfil snps are in attr files
max(sapply(stats.attr, nrow))
nrow(stats.attr[[5]])

#check snp of int
stats.attr[[1]][grep('rs4665765', stats.attr[[1]]$rsid)]
test <- stats.attr[[1]][all_BB>=all_AA]
#get the list of snps to reverse beta
head(test)
range(test$all_maf)
####
                
# get maf - need to use maf only if EAF is the minor allele
for (i in 1: length(stats.attr)){
  stats.attr[[i]][all_BB>=all_AA, ':=' (all_maf=1-all_maf)]
  stats.attr[[i]] <- stats.attr[[i]][,c(1,2,3,4,5,6,10)]
}

head(stats.attr[[1]])
# merge each maf col to cond
names(stats.attr[[1]])

# need to do this manually
sum.dat <- merge(stats.attr[[5]], stats.attr[[4]][,c(1,6)], by='rsid', all.x = T, suffixes=c('.gwas5', '.gwas4'))
sum.dat <- merge(sum.dat, stats.attr[[3]][,c(1,6)], by='rsid', all.x = T)
sum.dat <- merge(sum.dat, stats.attr[[2]][,c(1,6)], by='rsid', all.x = T, suffixes=c('.gwas3', '.gwas2'))
sum.dat <- merge(sum.dat, stats.attr[[1]][,c(1,6)], by='rsid', all.x = T)
head(sum.dat)

#get meta final snps
#gwas.meta <- read.table('../../AMLmeta_results/status_NCL_PCspeAMLHRC.meta.gz', header = T)
# merge mean.gwas and other cols to gwas.meta
#gwas.meta.merge <- merge(gwas.meta, sum.dat, by.x='BP', by.y='position')
#gwas.meta.merge[grep('25362515',gwas.meta.merge$BP),]
gwas.meta.merge <- sum.dat
head(gwas.meta.merge)
gwas.meta.merge$mean.info <- rowSums(gwas.meta.merge[,c(6,8:11)])/5
range(gwas.meta.merge$mean.info, na.rm = T)
#save as RDS first
gwas.meta.merge$SNP <- paste0(gwas.meta.merge$chromosome,':',gwas.meta.merge$position, '_',gwas.meta.merge$alleleA,
                              '_' ,gwas.meta.merge$alleleB)

saveRDS(gwas.meta.merge, 'GWAS_info_merge.RDS')
# output only high info snps. use 0.98
gwas.meta.merge <- subset(gwas.meta.merge, gwas.meta.merge$mean.info > 0.98)
range(gwas.meta.merge$mean.info)
tail(gwas.meta.merge)

# take sum.dat.gwas snps for each chr and slurm to plink extract to keep only meta + high r2 snps from merged bfiles
# use each filtered bfiles files to prune LD using plink - 50 ...
write.table(gwas.meta.merge$SNP, 'GWAS_r2filtered.snplist', quote = F, col.names = F, row.names = F)

#delete
test <- readRDS('GWAS_info_merge.RDS')
head(out)