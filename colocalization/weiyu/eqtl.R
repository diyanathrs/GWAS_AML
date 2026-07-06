# analyse weiyu's eqtl method 
library(GenomicRanges)
library(regioneR)
library(dplyr)
library(xtable)
library(EnsDb.Hsapiens.v75)

# load ukb snps
file <- '../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab'
hrcfilein <- paste0("cut -f1-5 ../../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# snps to check
# 808 - 2p23 status
snps <- c("rs2164808", "rs4930561","rs11481", "rs79918355" , "rs12988876", , "rs12078864"  ,"rs3916765" ,"rs12078864" ,"rs79918355" , "rs12988876"  ,"rs12632224" ,
          "rs11212666" ,"rs2853677" , "rs7705526")
snps <- c("rs2164808"  , "rs11481" , "rs3916765", "rs79918355" , "rs12988876" , "rs12078864")
# keep only snps in above list
hrc.rsid.snps <- hrc.rsid[hrc.rsid$ID %in% snps,]

#load files - weiyu
# status sig1- chr2-  rs2164808
files <- list.files(pattern = '*R4.datainplots.SuppTable2.txt')
eqtl <- lapply(files, read.table, header = T)
for (i in 1:length(eqtl)) {
  eqtl[[i]]$padj <- p.adjust(eqtl[[i]]$eQTLPvalue, method = 'BH')
}

eqtl <- rbindlist(eqtl)
names(eqtl)
#check
stopifnot(snps %in% eqtl$SNPrsID)
eqtlAML <- eqtl[eqtl$SNPrsID %in% snps]
eqtlAML <- eqtlAML[,c(2,4,10,3,5,6,8,1,14,9)]
eqtlAML <- eqtlAML[order(eqtlAML$eQTLPvalue)]
names(eqtlAML) <- c("SNP","SNP Position", "Assessed allele eQTL", "Chr", "Gene", "Gene Symbol", "Gene Position", "P value eQTL", "PBH eQTL", "Z score")
eqtlAML$`P value eQTL` <- formatC(eqtlAML$`P value eQTL`, format = "E", digits = 3)
eqtlAML$`PBH eQTL` <- formatC(eqtlAML$`PBH eQTL`, format = "E", digits = 3)
# required cols - SNP
#SNP, Positiona ,Assessed allele, eQTL, Chr, Gene, Gene Symbol, Gene Positiona, P value, eQTLb PBH eQTLc, Z-score
for (snp in snps){
  eqtlSNP <- subset(eqtlAML, eqtlAML$SNP==snp)
  eqtlSNP[2:nrow(eqtlSNP),1:4] <- ' '
  print(xtable(eqtlSNP, type = "latex"), file = paste0('../latex/',snp,"_AMLeQTL.tex"), include.rownames=FALSE)
}

##############
## FUMA data
##############
fuma <- read.csv('eQTLresult_rs115890122_2025-07-04.csv', header = T)
head(fuma)
fuma$gwasP <- formatC(fuma$gwasP, format = "E", digits = 2)
# sort and remove NA
unique(fuma$IndSigSNP)
snp.ord <- c("rs2164808", "rs11481" , "rs79918355" , "rs12988876", "rs12078864")
fuma <- fuma[fuma$IndSigSNP %in% snp.ord,]
unique(fuma$IndSigSNP)
fuma <- fuma %>% arrange(factor(IndSigSNP, levels = snp.ord))
#fuma <- fuma[order(fuma$IndSigSNP),]
fuma$gwasP[fuma$gwasP == ' NA'] <- ''
print(xtable(fuma, type = "latex"), file = 'fuma_results.tex', include.rownames=FALSE)

################
## eQTL from outside
###################
eqtlAML <- read.csv('eQTLresult_rs12988876_2025-07-04.csv', header = T, tryLogical = F)
head(eqtlAML)
eqtlAML <- eqtlAML[,c(6,8,10,7,3,4,5,1,2,9)]
eqtlAML <- eqtlAML[order(eqtlAML$Nominal_Pval),]

eqtlAML$Nominal_Pval <- formatC(eqtlAML$Nominal_Pval, format = "E", digits = 3)
eqtlAML$FDR <- formatC(eqtlAML$FDR, digits = 2, decimal.mark = 2)
# required cols - SNP
#SNP, Positiona ,Assessed allele, eQTL, Chr, Gene, Gene Symbol, Gene Positiona, P value, eQTLb PBH eQTLc, Z-score
head(eqtlAML)
print(xtable(eqtlAML, type = "latex"), file = paste0('../latex/',snp,"_AMLeQTL.tex"), include.rownames=FALSE)

