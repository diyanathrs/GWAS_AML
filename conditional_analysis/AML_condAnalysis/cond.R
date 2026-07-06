# cond analysis
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
# load ukb snps
file <- '../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab'
hrcfilein <- paste0("cut -f1-5 ../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
hrc.rsid <- subset(hrc.rsid, hrc.rsid$`#CHROM`==11)
head(hrc.rsid)

# start from here
cond <- read.table('AMLmetaCond_rs11481.meta', header = T)
#cond.rs4930561 <- read.table('AMLmetaCond_rs4930561.meta', header = T)
#cond <- cond.rs4930561
head(cond)
ori <- read.table('../AMLmeta_results/status_NCL_PCspeAMLHRC.meta.gz', header = T)
head(ori)
cond <- subset(cond, cond$N >= 4)
#remove if ref or alt is empty
cond[grep('67931761',cond$SNP),] # present
cond[grep('67820335',cond$SNP),] # this is rs11481 and shouldnt be there
cond.rs4930561[grep('67820335',cond.rs4930561$SNP),] # present
cond.rs4930561[grep('67931761',cond.rs4930561$SNP),]
#in ori data
ori[grep('67931761',ori$SNP),]
ori[grep('67820335',ori$SNP),]
head(cond)

#map snp_ids
cond <- merge(cond, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
head(cond)
#remove non rsid cols 
#cond <- subset(cond, !cond$ID=='.')
#table(cond$A1)
# original status data
#map snp_ids
ori <- merge(ori, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
head(ori)

# get locus zoom data
loc <- locus(data = cond, chrom = 'CHR',pos = 'BP', p = 'P',gene = 'SUV420H1', labs = 'ID', flank = 0.5e6, ens_db = "EnsDb.Hsapiens.v75", LD = 'r2')
loc <- link_LD(loc, token = "fc368a631dc0")
loc <- link_recomb(loc)
#loc2 <- link_eqtl(loc, token = "fc368a631dc0")
#head(loc2$LDexp)

plt.cond <- locus_ggplot(loc, labels = c("index", "rs11481", "rs4930561"), filter_gene_biotype = 'protein_coding', heights = c(2,1))
#rs4930561

ori.loc <- locus(data = ori, chrom = 'CHR',pos = 'BP', p = 'P',gene = 'SUV420H1', labs = 'ID', flank = 0.5e6, ens_db = "EnsDb.Hsapiens.v75", LD = 'r2')
ori.loc <- link_LD(ori.loc, token = "fc368a631dc0")
ori.loc <- link_recomb(ori.loc)
#loc2 <- link_eqtl(loc, token = "fc368a631dc0")
#head(loc2$LDexp)

plt.ori <- locus_ggplot(ori.loc, labels = c("rs4930561", "rs11481"), filter_gene_biotype = 'protein_coding', heights = c(2,1)) 

plot_grid(plt.ori, plt.cond, labels = c('Original data','Cond. on rs11481'))
