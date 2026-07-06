## polyfun final plots and sumstats for FUMA

#####################
##start polyfun now #
#####################
library(ggplot2)
library(ggrepel)
library(cowplot)
# load unfiltered sumstat and use it to filter finemap output
sumstat <- read.table('../AMLmeta_results/del5_7_NCL_PCspeAMLHRC.meta.gz', header = T)

list.files(path = 'results_out', pattern = 'susie.out$', full.names = T)
hit <- 'del5_7_hit1_chr1' 
folder <- 'new'

#### names ####
#PanAML 2p23.3 locus

############
#old <- polyfun_plt(hit,'results_out', '(sum.stats)')
#new <- polyfun_plt(hit,'new', '(info=0.3)')
#no.fun <- polyfun_plt.nofun(hit,'results_nofun', '(No Func)')
#add snps of interest
#custom.snps <- c('rs4930561')
#snps.of.int <- dat.susie[dat.susie$SNP==custom.snps,]

min.pip <- 0.05
title <- 'Del(5/7) AML 1q23.2 locus'
#title <- 'Complex karyotype AML 2p21 locus'
#title <- 'Complex karyotype AML 2q33.3 locus'
pdf(file = paste0(hit,'.pdf'),width = 7, height = 5.5)
#png(filename =paste(hit,'.png'), width = 10, height = 6, res = 300, units = 'in' )
polyfun_plt(hit,'new', title, min.pip)
dev.off()

polyfun_plt <- function(hit, folder, title, min.pip){
  dat.susie <- read.table(paste0(folder,'/',hit,'_susie.out'), header = T)
  dat.finemap <- read.table(paste0(folder,'/',hit,'_finemap.out'), header = T)
  
  #check all snps are there in both susie and finemp
  setdiff(dat.finemap$id, dat.susie$id)
  
  # merge N column from sumstats - use bp instead of snp
  dat.susie <- merge(dat.susie, sumstat[,c(2,6)], by='BP')
  dat.finemap <- merge(dat.finemap, sumstat[,c(2,6)], by='BP')
  # remove N!=Nmax
  dat.susie <- subset(dat.susie, dat.susie$N.y==max(dat.susie$N.y))
  dat.finemap <- subset(dat.finemap, dat.finemap$N==max(dat.finemap$N))
  
  dat.susie$CREDIBLE_SET <- as.factor(dat.susie$CREDIBLE_SET)
  # dat.susie$CREDIBLE_SET[dat.susie$CREDIBLE_SET==0] <- ''
  # dat.susie$CREDIBLE_SET <- factor(dat.susie$CREDIBLE_SET, levels = c("Background", "1", "2", "3"))
  # dat.finemap$CREDIBLE_SET <- as.character(dat.finemap$CREDIBLE_SET)
  # dat.finemap$CREDIBLE_SET[dat.finemap$CREDIBLE_SET==0] <- 'Background'
  # dat.finemap$CREDIBLE_SET <- factor(dat.finemap$CREDIBLE_SET, levels = c("Background", "1", "2", "3"))
  
  # merge P and SNPvar from susie
  dat.finemap <- merge(dat.susie[c(1,9)], dat.finemap, by='BP')
  sus.annot <- dplyr::filter(dat.susie, PIP > min.pip & P < 5e-6)
  fin.annot <- dplyr::filter(dat.finemap, PIP > min.pip & P < 5e-6)
  
  ## SUSIE result
  susie.new <- ggplot(dat.susie, aes(BP/1e6)) +
    geom_point(data= subset(dat.susie,dat.susie$CREDIBLE_SET==1), aes(y=-log10(P),size=PIP), alpha=0.8, 
               col='#a1226c') +
    geom_point(data= subset(dat.susie,dat.susie$CREDIBLE_SET==0), aes(y=-log10(P),size=PIP), alpha=0.8, 
               col='gray')+
    geom_label_repel(data = sus.annot, aes( label=SNP, y = -log10(P)), size=3,  max.overlaps = 25, box.padding = 1, direction = 'x')+  theme_bw()+
    ggtitle(paste(title, "- SuSiE"))+ theme(legend.position = 'none', axis.title.x = element_blank())
  
  
  
  susie <- ggplot(dat.susie, aes(BP/1e6)) +
    #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
    geom_point(aes(y=-log10(P), size=PIP, color=CREDIBLE_SET), alpha=0.7, shape=20 ) + 
    geom_label_repel(data = sus.annot, aes( label=SNP, y = -log10(P)), size=3,  max.overlaps = 25, box.padding = 1, direction = 'x')+  theme_bw()+
    ggtitle(paste(title))+ theme(legend.position = 'none', axis.title.x = element_blank())+ 
    scale_color_brewer(palette = 'Set2') 
  
  ## FINEMAP result ##
  finemap <- ggplot(dat.finemap, aes(BP/1e6)) +
    geom_point(aes(y=-log10(P), size=PIP, color=CREDIBLE_SET), alpha=0.7, shape=20) + 
    geom_label_repel(data= fin.annot, aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 25, box.padding = 1, direction = 'x') + theme_bw()+
    ggtitle(paste(title, "- FINEMAP"))+ theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank())+
    scale_color_brewer(palette = 'Set2') 
  
  
  #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
  #plot_grid(finemap, susie)

 fmp <- ggplot(dat.finemap, aes(BP/1e6)) +
    geom_bar(aes(y=PIP, fill=CREDIBLE_SET), stat = 'identity', width=0.02)+ theme_bw()+
   theme(legend.position = 'bottom',axis.title.y = element_blank()) + xlab(paste('Chr',dat.finemap$CHR, '(Mb)')) +
   scale_fill_brewer(palette = 'Set2') + guides(fill=guide_legend(title = 'Credible set'))
 
 sie <- ggplot(dat.susie, aes(BP/1e6)) +
   geom_bar(aes(y=PIP, fill=CREDIBLE_SET), stat = 'identity', width=0.02)+ theme_bw()+
   theme(legend.position = 'none') + xlab(paste('Chr',dat.susie$CHR, '(Mb)')) +
   scale_fill_manual(values = c('gray', '#a1226c')) 
 
  
  plot_grid(susie.new,sie, ncol = 1, rel_heights = c(2,0.8), align = 'v')
  #f <- plot_grid(finemap,fmp, ncol = 1, rel_heights = c(2,0.8), align = 'v')
  
  #plot_grid(s,f)
}

polyfun_plt_custom <- function(hit, folder, title, snps.of.int){
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  
  dat.susie <- read.table(paste0(folder,'/',hit,'_susie.out'), header = T)
  dat.susie$CREDIBLE_SET <- as.character(dat.susie$CREDIBLE_SET)
  snps.of.int <- dat.susie[dat.susie$SNP==custom.snps,]
  
  dat.finemap <- read.table(paste0(folder,'/',hit,'_finemap.out'), header = T)
  dat.finemap$CREDIBLE_SET <- as.character(dat.finemap$CREDIBLE_SET)
  # merge P and SNPvar from susie
  dat.finemap <- merge(dat.susie[c(3,6,9)], dat.finemap, by='BP')
  
  ## SUSIE result
  susie <- ggplot(dat.susie, aes(BP/1e6)) +
    #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
    geom_point(aes(y=-log10(P), size=PIP, fill=CREDIBLE_SET), alpha=0.7, shape=21 ) + 
    geom_label_repel(data= subset(dat.susie, dat.susie$PIP > 0.05), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20
                     ,box.padding = 0.5)+ 
    geom_label_repel(data= snps.of.int, aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20
                     ,box.padding = 0.5)+ theme_bw() +
    ggtitle(paste(hit, "- SUSIE",title))+ theme(legend.position = 'none', axis.title.x = element_blank())
  
  ## FINEMAP result ##
  finemap <- ggplot(dat.finemap, aes(BP/1e6)) +
    geom_point(aes(y=-log10(P), size=PIP, fill=CREDIBLE_SET), alpha=0.7, shape=21) + 
    geom_label_repel(data= subset(dat.finemap, dat.finemap$PIP > 0.05), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20
                     ,box.padding = 0.5,)+ 
    geom_label_repel(data= snps.of.int, aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20
                     ,box.padding = 0.5)+ theme_bw() +
    ggtitle(paste(hit, "- FINEMAP",title))+ theme(legend.position = 'none', axis.title.x = element_blank())
  
  
  #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
  #plot_grid(finemap, susie)
  
  
  fmp <- ggplot(dat.finemap, aes(BP/1e6)) +
    geom_bar(aes(y=PIP, fill=CREDIBLE_SET), stat = 'identity', width=0.02)+ theme_bw()+
    theme(legend.position = 'bottom') + xlab(paste('Chr',dat.finemap$CHR, '(Mb)')) 
  
  sie <- ggplot(dat.susie, aes(BP/1e6)) +
    geom_bar(aes(y=PIP, fill=CREDIBLE_SET), stat = 'identity', width=0.02)+ theme_bw() +
    theme(legend.position = 'bottom') + xlab(paste('Chr',dat.susie$CHR, '(Mb)')) 
  
  
  s <- plot_grid(susie,sie, ncol = 1, rel_heights = c(2,0.8), align = 'v')
  f <- plot_grid(finemap,fmp, ncol = 1, rel_heights = c(2,0.8), align = 'v')
  
  plot_grid(s,f)
  
}

###############
# get new hits
##############
library(dplyr)
files <- list.files(path = 'new', pattern = 'susie.out$', full.names = T)
pheno <- c("CBF1", "Complex 2p21", "Complex 2q33.3", "Del(5/7) 1q23.2", "PanAML 2p23.3", "PanAML 11q13.2", "t57")
all.hits <- lapply(files, read.delim)
names(all.hits) <- pheno
all.hits <- bind_rows(all.hits, .id = 'id')
all.hits$SNPVAR <- formatC(all.hits$SNPVAR, format = "E", digits = 2)
all.hits$P <- formatC(all.hits$P, format = "E", digits = 2)
all.hits$BETA_MEAN <- formatC(all.hits$BETA_MEAN, format = "E", digits = 2)
all.hits$BETA_SD <- formatC(all.hits$BETA_SD, format = "E", digits = 2)
# remove unwanted pheno
hits2tab <- all.hits[all.hits$id %in% c("Complex 2p21", "Complex 2q33.3", "Del(5/7) 1q23.2", "PanAML 2p23.3", "PanAML 11q13.2"),]
hits2tab <- subset(hits2tab, hits2tab$PIP > 0.01)
#write tex
print(xtable(hits2tab[,-c(8,14)], type = "latex"), file = 'finemap_results.tex', include.rownames=FALSE)

#snps.fuma <- subset(all.hits,all.hits$PIP > 0.2)
write.table(snps.fuma, 'Finemap_snps.txt', quote = F, row.names = F, col.names = T)

# check maf for rare variants
ncl.gwas <- readRDS('~/AML/AML_condAnalysis/cojo/gwas_attr/NCL7_Complex.rds')
head(ncl.gwas)
snp <- '2:206048000'
ncl.gwas[grep(snp,ncl.gwas$rsid)]
