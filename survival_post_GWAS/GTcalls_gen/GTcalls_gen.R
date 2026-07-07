# Generate GT calls for Holly
# used VCFs generated using QCtool with 0.9 threshold
########################
## GT calls for holly ##
########################
library(parallel)
library(survminer)

#snp.list <- read.table('snp_list_holly.txt', header = T)
#add new SNPs
rsid <- c('rs2807011', 'rs9397136', 'rs112121578', 'rs4916944', 'rs16864108', 'rs3006923', 'rs147189138')
chr <- c(1, 6, 11, 7, 2, 1, 11)
Position <- c(31326404, 153556877, 1937116, 75835, 223615152, 243648757, 121149360)
snp.list <- data.frame(rsid,chr,Position)
snp.list <- snp.list %>% arrange(chr, Position)

# check rsid against coords
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS)]
head(hrc.rsid)
# check
hrc.rsid <- hrc.rsid[hrc.rsid$rsid %in% paste0(snp.list$chr, ':',snp.list$Position)]
stopifnot(snp.list$rsid == hrc.rsid$ID)

pheno <- read.delim('../../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) 

infer_GT <- function(n) {
  snp <- snp.list[n,]
  vcf.names <- read.delim(paste0('../vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
  chr.vcf <- fread(paste0('grep -w ',snp$Position, ' ../vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
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
# plot KM
plot_km <- function(snp.list) {
  pheno <- read.delim('../../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T)# %>% filter(!.id=='BirminghamNCL5 ')
  #n <- 1
  pdf(file = paste0('AMLsur_KMplots_',Sys.Date(),'.pdf'),width = 12, height = 7,onefile = T )
  #pdf(file = 'AMLsur_plots.pdf',width = 12, height = 7, onefile = T )
  
  for (n in 1:nrow(snp.list)) {
    snp <- snp.list[n,]
    vcf.names <- read.delim(paste0('../vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'), skip=3, header = T,stringsAsFactors = F, check.names = F, nrows = 1)
    chr.vcf <- fread(paste0('grep -w ',snp$Position, ' ../vcf/AMLsur_chr',snp$chr,'_Allmerged.vcf'),sep = '\t', header = F) 
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
                 risk.table.height=.2, risk.table.y.text.col = T, risk.table.y.text = FALSE)
    
    sur.eln2 <- survfit(Surv(OSdays/365.25, OSstatus) ~Genotype, data = snp.interest.eln2) %>% 
      ggsurvplot(conf.int=F, pval=TRUE, risk.table=TRUE, 
                 xlab='Time in years',
                 legend.title=paste(snp$rsid), palette=c("red", "blue", "forestgreen"), 
                 title=paste0("ELN2 only -",' chr:' ,snp$chr,snp$Position,'_' ,ref,'>',alt,' - ' ,snp$rsid),
                 risk.table.height=.2, risk.table.y.text.col = T, risk.table.y.text = FALSE)
    
    plot <- cowplot::plot_grid(sur.all$plot,sur.eln2$plot, sur.all$table, sur.eln2$table, rel_heights = c(6,2))
    
    print(plot)
  }
  dev.off()
}

# run inferGT
test <- parallel::mclapply(seq_len(nrow(snp.list)), infer_GT, mc.cores = 7)
out <- as.data.frame(do.call(cbind,test))
out$sample_id <- row.names(out)
dt <- Sys.Date()
write.table(out, paste0('AMLsur_GTcalls_',dt,'.txt'), quote = F, row.names = F, col.names = T, sep = '\t')

plot_km(snp.list)
