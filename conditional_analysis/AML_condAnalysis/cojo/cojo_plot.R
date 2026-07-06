#######################
## post cond analysis##
#######################
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(cowplot)

setwd('cond_out_nofiltr/')

# Plot sup fig 10 - 15
# load ukb snps
hrcfilein <- paste0("cut -f1-5 ../../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, ",",`#CHROM`, ":", POS,':', REF, ":", ALT)]
head(hrc.rsid)


# Status 11q13.2 - sup fig 10 rs4930561
# Status 2p23.3 - sup fig 11 
# sup fig 12 - Normal 6p21.32
# sup fig 13 - del5_7 1q23.2
# sup fig 14 - complex 2p21
# sip fig 15 - complex 2q33.3

pheno <- 'status'
# load cond results
# only use cond_c
cond <- read.table('coj', header = T)
normal.jma <- read.table('cojo_slct_complex_chr2.jma.cojo', header = T)
#cond_c <- read.table('cojo_slct_del57_chr1.cma.cojo', header=T)

# Combine joint + comb
setdiff(names(normal.jma), names(cond))
setnames(normal.jma, c("bJ", "bJ_se", "pJ"), c("bC", "bC_se", "pC"))
# drop LD_r
names(cond)
names(normal.jma)
setdiff(names(cond), names(normal.jma))
normal.jma <- normal.jma[-14]
cond <- rbind(normal.jma, cond)

#map snp_ids
cond <- merge(cond, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
head(cond)
#check snp
cond[which(cond$ID=='rs115890122'),]
# remove snps with pC=NA
cond <- cond[-c(which(is.na(cond$pC))),]
range(cond$pC)
# for normal remove rs79657479
rmv <- which(cond$ID=='rs3916765')
cond <- cond[-rmv,]
cond[(cond$ID=='rs3916765'),]

# Use locus zoom
# prad3b
loc <- locus(data = cond, chrom = 'Chr',pos = 'bp', p = 'pC', gene = 'PARD3B', labs = 'ID', flank = 0.5e6, ens_db = "EnsDb.Hsapiens.v75")
loc <- link_LD(loc, token = "fc368a631dc0")
loc <- link_recomb(loc)
#loc2 <- link_eqtl(loc, token = "fc368a631dc0")
#head(loc2$LDexp)
# see remaining peaks
subset(cond, cond$pC < 5e-5)
ensembldb::listGenebiotypes(EnsDb.Hsapiens.v75)

#trace(locus_plot, edit = T)
locus_plot(loc, labels = c("index"), filter_gene_biotype = c('protein_coding','lincRNA'), heights = c(2,1), label_x = c(2, 2), 
           label_y = c(3,9))

pdf(paste0(pheno,'_cond_supfig.pdf'),width = 7, height = 5.5)
locus_plot(loc, labels = c("index"), filter_gene_biotype = c('protein_coding', 'lincRNA'), heights = c(2,1), label_x = c(-2, 1), 
           label_y = c(4,11))
dev.off()

#########################
#### delete ##
#########################
#check custom snp - normal rs1794275
cond <- read.table('cojo_slct_Normal_chr6.cma.cojo', header = T)
normal.jma <- read.table('cojo_slct_Normal_chr6.jma.cojo', header = T)
cond <- merge(cond, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
normal.jma <- merge(normal.jma, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')

head(cond)
cond[cond$ID=='rs1794275',]

#calculate OR and 95%CI
normal.jma$OR <- exp(normal.jma$b)
normal.jma$OR_lower <- exp(normal.jma$b - 1.96 * normal.jma$se)
normal.jma$OR_upper <- exp(normal.jma$b + 1.96 * normal.jma$se)

normal.jma[normal.jma$ID=='rs1794275',]
-log10(1.91968e-05)

# check complex low maf snp - rs115890122
cplx.in <- read.table('../Complex_cond_in.ma', header = T)
#cplx.in <- merge(cplx.in, hrc.rsid[,c(3,6)], by.x="SNP", by.y='rsid')
head(cplx.in)
cplx.in[grep('206048000', cplx.in$SNP),]
#which(hrc.rsid$ID=='rs115890122')

# del5/7 snps
del57.snps <- c("rs12078864", "rs12078861", "rs12063642" , "rs7519478" ,"rs4406626" ,"rs3806184", "rs3820099", "rs11801302")
#check snp
hrc.rsid[which(hrc.rsid$ID %in% del57.snps),]$rsid

