library(qqman)
library(data.table)
library(dplyr)
library(genom)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(topr)

#####################
## manhatton plots ##
#####################
# load ukb snps
hrcfilein <- paste0("cut -f1-5 HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# load meta
pheno <- 'Normal'
status <- read.table(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC.meta.gz'),header = T)
n.max <- max(status$N)
status <- subset(status, status$N==n.max)
names(status)[3] <- 'rsid'
head(status)
#map snp_ids
merged.out <- merge(status, hrc.rsid, by="rsid", all.x=T)
tail(merged.out)
# check gwas sig hits
merged.out[merged.out$P < 1e-7,]
# more pruning
#merged.out <- merged.out %>%  filter(-log10(P)>1)

if (pheno=='PanAML') {
# add genes and cytoband for status
lead.idx <- c('rs4665765', 'rs11481')
merged.out[grep(lead.idx[1],merged.out$ID),]$rsid <- paste0('rs4665765\n2p23.3',' (EFR3B, DNMT3A)')
merged.out[grep(lead.idx[2],merged.out$ID),]$rsid <- paste0('rs11481\n11q13.2',' (CHKA)')
} else if (pheno=='Normal') {
# add genes and cytoband for Normal - Plot rs3916765 instead of rs3997854
lead.idx <- c('rs3997854', 'rs1794275')
#remove the rest of the rsid
merged.out$rsid <- ''
merged.out$rsid[grep(lead.idx[1],merged.out$ID)] <- paste0('rs3916765\n6p21.32',' (HLA-DQA2)')
merged.out$rsid[grep(lead.idx[2],merged.out$ID)] <- paste0('rs1794275\n6p21.32',' (HLA-DQB1)')

} else if (pheno=='Complex') {
# add genes and cytoband for Complex
lead.idx <- c('rs79918355', 'rs12988876')
merged.out[grep(lead.idx[1],merged.out$ID),]$rsid <- paste0('rs79918355\n2p21',' (EPCAM)') #('2p21','(EPCAM-DT)'
merged.out[grep(lead.idx[2],merged.out$ID),]$rsid <- paste0('rs12988876\n2q33.3',' (PARD3B)') # '(PARD3B)'
#remove 2:47512869_T_C
merged.out[grep('2:47512869_T_C',merged.out$rsid),]$rsid <- ''
} else {
  # for del57
  lead.idx <- c('rs12078864')
  merged.out[grep(lead.idx, merged.out$ID),]$rsid <- paste0('rs12078864\n1q23.2',' (DUSP23)')
}
# del5_7 sig
# what to annotate
#subset(merged.out, merged.out$P<5e-6)

merged.out[which(merged.out$ID=='rs1794275'),]

source('subsample_snps.R')
#test <- subsample_gwas(merged.out, pcol = "P", chrcol = "CHR", poscol = "BP", keep_threshold = 1e-1, max_per_chr = 3000)

## test plot with less points
pdf(paste0(pheno,"_mhHRCres_test.pdf"), width=12, height=7)
qqman::manhattan(subset(merged.out, merged.out$P<5e-3), chr="CHR", bp="BP", snp="rsid", p="P" ,
          suggestiveline = FALSE, col = c('#0b590f','#73a647'), ylim=c(0, 10),
          annotatePval = 5e-6, annotateTop = F) # toggle annotateTop

dev.off()
#ori colors - #161882','#6668d4 
trace(manhattan, edit = T) # change cex size to 1

# real plotting
tiff(paste0(pheno,"_mhHRC.tiff"),  width=12, height=7, units = 'in', res = 500)
qqman::manhattan(merged.out, chr="CHR", bp="BP", snp="rsid", p="P" ,suggestiveline = FALSE, annotatePval = 5e-6, 
          annotateTop = F, col = c('#0b590f','#73a647'),ylim=c(0, 10))
dev.off()



# old script
source('manhattan.plot.r')
png(paste0("status_mhHRCres2.png"), type="cairo", width=14, height=9, units="in", res=300)
manhattan.plot(merged.out,chr = 'CHR', bp = 'BP', snp = 'rsid', p = 'P',  lty=2, lcol="firebrick", cex = 2)
dev.off()

q()
###################
##Plot QQ using gwas attr in cond
########################
# weiyu used tabix data
pheno <- 'normal'
plotQQ <- function(z,color){
  p <- 2*pnorm(-abs(z))
  p <- sort(p)
  expected <- c(1:length(p))
  lobs <- -(log10(p))
  lexp <- -(log10(expected / (length(expected)+1)))
  
  # plots all points with p < 0.05
  p_sig = subset(p,p<0.05)
  points(lexp[1:length(p_sig)], lobs[1:length(p_sig)], pch=21,
         cex=0.8, col=color)
  #, bg=color)
  
  # samples 5,000 points from p > 0.05 (to keep file size down)
  n=5001
  i<- c(length(p)- c(0,round(log(2:(n-1))/log(n)*length(p))),1)
  lobs_bottom=subset(lobs[i],lobs[i] <= 3)
  lexp_bottom=lexp[i[1:length(lobs_bottom)]]
  points(lexp_bottom, lobs_bottom, pch=23, cex=0.8, col=color)#, bg=color)
}

attr.lst <- list.files('../AML_condAnalysis/cojo/gwas_attr', pattern = paste0('del5'), 
                       full.names = T,ignore.case = T)
pheno.attr <- parallel::mclapply(attr.lst, readRDS, mc.cores = length(attr.lst))

for (i in seq_along(attr.lst)) {
gwas <- sub('../AML_condAnalysis/cojo/gwas_attr/', '', attr.lst[i])
gwas <- sub('.rds', '', gwas)
S <- pheno.attr[[i]]
names(S)
setnames(S, 'frequentist_add_pvalue', 'P',skip_absent = T)
#filter on maf and info
pheno.attr[[i]] <- pheno.attr[[i]] %>% dplyr::filter(info > 0.6 & all_maf > 0.02 & !is.na(P))
pheno.attr[[i]]

# calculate lambda
z=qnorm(pheno.attr[[i]]$P/2)
lambda = round(median(z^2,na.rm=T)/qchisq(0.5,df=1),3)
print(lambda)

###############
## plotQQ weiyu
###############
## Plot function ##
xmax <- -log10(S[which.min(P)]$P)
e = -log10(ppoints(length(S$P)))
# plot
png(file = paste0(gwas,'_QQ.png'), width = 6, height = 5, units = 'in', res = 250, type = "cairo")
plot(x=-log10(ppoints(length(S$P))) , y=-log10(S$P),
     type="n",
     xlim = c(0, max(e)+0.2), ylim= c(0, xmax+0.2),
     xlab= expression(Expected ~ ~-log[10](italic(P))),
     ylab= expression(Observed ~ ~-log[10](italic(P))),
     las=1, xaxs="i", yaxs="i", bty="l",
     main=c(substitute(paste("",lambda," = ", lam),list(lam = lambda)),expression()))
abline(0, 1, col= "black", lwd=2)
plotQQ(z, 'red')
## provides legend
legend(x="topleft",legend=c("Expected","Observed"),
       pch=c(NA, 21),
       cex=1, lty=c(1, NA), lwd=c(2, NA),
       col = c("black", "#e41a1c"),
       #pt.bg=c(NA, color),
       bty="n")

dev.off()

}


####################
#### use ggraster ##
#####################
library(ggplot2)
library(ggrastr)

# load ukb snps
hrcfilein <- paste0("cut -f1-5 HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# load meta
pheno <- 'Normal'
status <- read.table(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC.meta.gz'),header = T)
n.max <- max(status$N)
status <- subset(status, status$N==n.max)
names(status)[3] <- 'rsid'
head(status)
#map snp_ids
merged.out <- merge(status, hrc.rsid, by="rsid", all.x=T)
tail(merged.out)
# check gwas sig hits
merged.out[merged.out$P < 1e-7,]
# more pruning
#merged.out <- merged.out %>%  filter(-log10(P)>1)

if (pheno=='PanAML') {
  # add genes and cytoband for status
  lead.idx <- c('rs4665765', 'rs11481')
  merged.out[grep(lead.idx[1],merged.out$ID),]$rsid <- paste0('rs4665765\n2p23.3',' (EFR3B, DNMT3A)')
  merged.out[grep(lead.idx[2],merged.out$ID),]$rsid <- paste0('rs11481\n11q13.2',' (CHKA)')
} else if (pheno=='Normal') {
  # add genes and cytoband for Normal - Plot rs3916765 instead of rs3997854
  lead.idx <- c('rs3997854', 'rs1794275')
  #remove the rest of the rsid
  merged.out$rsid <- ''
  merged.out$rsid[grep(lead.idx[1],merged.out$ID)] <- paste0('rs3916765\n6p21.32',' (HLA-DQA2)')
  merged.out$rsid[grep(lead.idx[2],merged.out$ID)] <- paste0('rs1794275\n6p21.32',' (HLA-DQB1)')
  
} else if (pheno=='Complex') {
  # add genes and cytoband for Complex
  lead.idx <- c('rs79918355', 'rs12988876')
  merged.out[grep(lead.idx[1],merged.out$ID),]$rsid <- paste0('rs79918355\n2p21',' (EPCAM)') #('2p21','(EPCAM-DT)'
  merged.out[grep(lead.idx[2],merged.out$ID),]$rsid <- paste0('rs12988876\n2q33.3',' (PARD3B)') # '(PARD3B)'
  #remove 2:47512869_T_C
  merged.out[grep('2:47512869_T_C',merged.out$rsid),]$rsid <- ''
} else {
  # for del57
  lead.idx <- c('rs12078864')
  merged.out[grep(lead.idx, merged.out$ID),]$rsid <- paste0('rs12078864\n1q23.2',' (DUSP23)')
}
# del5_7 sig
# what to annotate
#subset(merged.out, merged.out$P<5e-6)

merged.out[which(merged.out$ID=='rs1794275'),]

head(merged.out)
merged.out$CHR <- factor(merged.out$CHR)

pdf("manhattan_raster.pdf")

ggplot(subset(merged.out, merged.out$P<5e-4), aes(BP, -log10(P), color = CHR %% 2)) +
  geom_point_rast(size = 0.3, raster.width = 3000, raster.height = 1500) +
  scale_color_manual(values = c("#4C72B0", "#DD8452")) +
  theme_bw()

dev.off()

topr::manhattan(merged.out )
