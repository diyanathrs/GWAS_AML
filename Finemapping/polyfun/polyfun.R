## plot polyfun results
install.packages('cowplot')
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
 
# correct IIDs for birmingham
bir <- read.delim('../AMLsur_Phenomerged_Dec2024_FINAL.txt')  %>% filter(Origin=="Birmingham") 
bir$IID <- paste0(bir$IID,'_',bir$IID)
setdiff(bir$IID, dat1$IID)

dat3 <- rbind(bir, dat2)

setdiff(dat1$IID, dat3$IID)
dup.dix <- which(duplicated(dat1$IID))
dat1 <- dat1[-dup.dix,]
dat1$IID[duplicated(dat1$IID)]

# merge dat3 anyscp data with dat1
dat4 <- merge(dat3[c(1,16)], dat1, by='IID')
table(dat4$Any.SCT)
dat4 %>% filter(Any.SCT==1) %>% group_by(.id) %>% summarise(n())
noSCT <- subset(dat4, dat4$Any.SCT==0)
#write.table(noSCT$IID, 'AML_subset_noSCT.txt', quote = F, row.names = F, col.names = F)

#####################
##start polyfun now #
#####################
list.files(path = 'results_out/', pattern = '\\.out$')
hit <- 'del5_7_hit1_chr1' 
#folder <- 'new'

old <- polyfun_plt(hit,'results_out', '(sum.stats)')
new <- polyfun_plt(hit,'new', '(info=0.3)')
#no.fun <- polyfun_plt.nofun(hit,'results_nofun', '(No Func)')

plot_grid(old,new)
#plot_grid(old,no.fun)

polyfun_plt <- function(hit, folder, title){
  
dat.susie <- read.table(paste0(folder,'/',hit,'_susie.out'), header = T)
dat.susie$CREDIBLE_SET <- as.character(dat.susie$CREDIBLE_SET)
#nofun
#dat.susie.nofun <- read.table(paste0(folder,'/',hit,'_susie_nofun.out'), header = T)
#dat.susie.nofun$CREDIBLE_SET <- as.character(dat.susie.nofun$CREDIBLE_SET)

dat.finemap <- read.table(paste0(folder,'/',hit,'_finemap.out'), header = T)
dat.finemap$CREDIBLE_SET <- as.character(dat.finemap$CREDIBLE_SET)
# merge P and SNPvar from susie
dat.finemap <- merge(dat.susie[c(3,6,9)], dat.finemap, by='BP')
# nofun
#dat.finemap.nofun <- read.table(paste0('results_nofun/',hit,'_finemap_nofun.out'), header = T)
#dat.finemap.nofun$CREDIBLE_SET <- as.character(dat.finemap.nofun$CREDIBLE_SET)
# merge P and SNPvar from susie
#dat.finemap.nofun <- merge(dat.susie.nofun[c(2,8)], dat.finemap.nofun, by='BP')

## SUSIE result
susie <- ggplot(dat.susie, aes(BP/1e6)) +
  #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
  geom_point(aes(y=-log10(P), fill=PIP, col=CREDIBLE_SET), alpha=0.8, shape=21, size=1.5) + 
  geom_label_repel(data= subset(dat.susie, dat.susie$P < 5e-6), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 30,
                   box.padding = 0.4)+ xlab(paste('Chr',dat.susie$CHR, '(Mb)')) +
  ggtitle(paste(hit, "- SUSIE",title))+ theme(legend.position = 'bottom')

## FINEMAP result ##
finemap <- ggplot(dat.finemap, aes(BP/1e6)) +
  #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
  geom_point(aes(y=-log10(P), fill=PIP, col=CREDIBLE_SET), alpha=0.7, shape=21, size=1.5) + 
  geom_label_repel(data= subset(dat.finemap,dat.finemap$P < 5e-6), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 30, 
                   box.padding = 0.4)+ xlab(paste('Chr',dat.susie$CHR, '(Mb)')) +
  ggtitle(paste(hit, "- FINEMAP", title)) + theme(legend.position = 'bottom')
  
#geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
#plot_grid(finemap, susie)

fmp <- ggplot(dat.finemap, aes(BP/1e6)) +
  geom_point(aes(y=PIP, col=CREDIBLE_SET), alpha=0.8)+  ggtitle(paste(hit, "- FINEMAP",title)) +
  geom_label_repel(data= subset(dat.finemap, dat.finemap$PIP > 0.1), aes(label=SNP, y = PIP), size=3, max.overlaps = 20,
                   box.padding = 0.4)+ theme(legend.position = 'bottom') + xlab(paste('Chr',dat.susie$CHR, '(Mb)')) 

## SUSIE result ##
sie <- ggplot(dat.susie, aes(BP/1e6)) +
  geom_point(aes(y=PIP, col=CREDIBLE_SET), alpha=0.8)+  ggtitle(paste(hit, "- SUSIE",title))+
  geom_label_repel(data= subset(dat.susie, dat.susie$PIP > 0.1), aes(label=SNP, y = PIP), size=3, max.overlaps = 20,
                   box.padding = 0.4)+ theme(legend.position = 'bottom') + xlab(paste('Chr',dat.susie$CHR, '(Mb)')) 

plot_grid(finemap, fmp, susie, sie , rel_widths = c(3,2))

}
polyfun_plt.nofun <- function(hit, folder, title){
  
  dat.susie <- read.table(paste0(folder,'/',hit,'_susie_nofun.out'), header = T)
  dat.susie$CREDIBLE_SET <- as.character(dat.susie$CREDIBLE_SET)
  #nofun
  #dat.susie.nofun <- read.table(paste0(folder,'/',hit,'_susie_nofun.out'), header = T)
  #dat.susie.nofun$CREDIBLE_SET <- as.character(dat.susie.nofun$CREDIBLE_SET)
  
  dat.finemap <- read.table(paste0(folder,'/',hit,'_finemap_nofun.out'), header = T)
  dat.finemap$CREDIBLE_SET <- as.character(dat.finemap$CREDIBLE_SET)
  # merge P and SNPvar from susie
  dat.finemap <- merge(dat.susie[c(2,8)], dat.finemap, by='BP')
  # nofun
  #dat.finemap.nofun <- read.table(paste0('results_nofun/',hit,'_finemap_nofun.out'), header = T)
  #dat.finemap.nofun$CREDIBLE_SET <- as.character(dat.finemap.nofun$CREDIBLE_SET)
  # merge P and SNPvar from susie
  #dat.finemap.nofun <- merge(dat.susie.nofun[c(2,8)], dat.finemap.nofun, by='BP')
  
  ## SUSIE result
  susie <- ggplot(dat.susie, aes(BP/1e6)) +
    #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
    geom_point(aes(y=-log10(P), fill=PIP, col=CREDIBLE_SET), alpha=0.8, shape=21, size=1.5) + 
    geom_label_repel(data= subset(dat.susie, dat.susie$P < 5e-6), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20,
                     box.padding = 0.4)+
    ggtitle(paste(hit, "- SUSIE", title))+ theme(legend.position = 'bottom')
  
  ## FINEMAP result ##
  finemap <- ggplot(dat.finemap, aes(BP/1e6)) +
    #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
    geom_point(aes(y=-log10(P), fill=PIP, col=CREDIBLE_SET), alpha=0.7, shape=21, size=1.5) + 
    geom_label_repel(data= subset(dat.finemap,dat.finemap$P < 5e-6), aes(label=SNP, y = -log10(P)), size=3, max.overlaps = 20, 
                     box.padding = 0.4)+
    ggtitle(paste(hit, "- FINEMAP", title)) + theme(legend.position = 'bottom')
  
  #geom_point(data= subset(status.susie,status.susie$CREDIBLE_SET==1), aes(y=-log10(P)), alpha=0.8, size=2, col='red') +
  #plot_grid(finemap, susie)
  
  fmp <- ggplot(dat.finemap, aes(BP/1e6)) +
    geom_point(aes(y=PIP, col=CREDIBLE_SET), alpha=0.8)+  ggtitle(paste(hit, "- FINEMAP", title)) +
    geom_label_repel(data= subset(dat.finemap, dat.finemap$PIP > 0.3), aes(label=SNP, y = PIP), size=3, max.overlaps = 20,
                     box.padding = 0.4)+ theme(legend.position = 'bottom')
  
  ## SUSIE result ##
  sie <- ggplot(dat.susie, aes(BP/1e6)) +
    geom_point(aes(y=PIP, col=CREDIBLE_SET), alpha=0.8)+  ggtitle(paste(hit, "- SUSIE", title))+
    geom_label_repel(data= subset(dat.susie, dat.susie$PIP > 0.3), aes(label=SNP, y = PIP), size=3, max.overlaps = 20,
                     box.padding = 0.4)+ theme(legend.position = 'bottom')
  
  plot_grid(finemap, fmp, susie, sie, rel_widths = c(3,2))
}


#shape=as.character(CREDIBLE_SET)
## for FUMA
hit.snps <-  read.table(paste0('new/',hit,'_susie.out'), header = T)
plot(hit.snps[-1,]$PIP)
write.table(hit.snps, paste0(hit,'_FUMA.txt'), sep = '\t', quote = F, row.names = F)

# for lead snp
cols <- c('rsID','chr','position')
lead.snp <- data.frame(rsid=hit.snps[1,]$SNP, chr=hit.snps[1,]$CHR, position=hit.snps[1,]$BP)
write.table(lead.snp, 'status_hit1_leadSNP.txt', sep = '\t', quote = F, row.names = F)
