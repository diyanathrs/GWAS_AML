# Download eqtl v10
#link <- 'https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar'
# link2 - https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz
# resources - eQTplot, LocusPlot, fastEnloc

library(GenomicRanges)
library(regioneR)
library()
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
# load ukb snps
# load ukb snps
file <- '../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab'
hrcfilein <- paste0("cut -f1-5 ../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# snps to check
snps <- c("rs4930561","rs10789158" , "rs17773014", "rs11481", "rs4665765"  ,"rs3916765" ,"rs12078864" ,"rs79918355" , "rs12988876"  ,"rs12632224" ,
          "rs11212666" ,"rs2853677" , "rs7705526")
# keep only snps in above list
hrc.rsid.snps <- hrc.rsid[hrc.rsid$ID %in% snps,]

#these need to be coloc'ed
gtex_eqtl <- read.delim('GTEx_Analysis_v7_eQTL/Whole_Blood.v7.egenes.txt.gz')
gtex_eqtl2 <- read.delim('GTEx_Analysis_v7_eQTL/Whole_Blood.v7.signif_variant_gene_pairs.txt.gz')
head(gtex_eqtl)
head(gtex_eqtl2)

names(gtex_eqtl2)
names(gtex_eqtl)

which(gtex_eqtl$rs_id_dbSNP147_GRCh37p13==snps[7])
gtex_eqtl[gtex_eqtl$gene_name=='DNAJC27',]
DNAJC27
head(gtex_eqtl)

# using mapped cis-eqtl directly
cis.bon <- read.delim('2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
head(cis.bon)
cis.bon[cis.bon$SNP==snps[n],]
names(cis.bon)
# get eqtl within 500kb region
gr.cis.eqtl <- GRanges(seqnames = cis.bon$SNPChr, ranges = IRanges(cis.bon$SNPPos,width = 1))
mcols(gr.cis.eqtl) <- cis.bon
# get range for snp - from hrc list
n <- 4
gr.snp <- toGRanges(paste0(hrc.rsid.snps[hrc.rsid.snps$ID==snps[n]]$`#CHROM`,':',hrc.rsid.snps[hrc.rsid.snps$ID==snps[n]]$POS,
                           '-',hrc.rsid.snps[hrc.rsid.snps$ID==snps[n]]$POS))
start(gr.snp) <- start(gr.snp)-2.5e5
end(gr.snp) <- end(gr.snp)+2.5e5
gr.snp

library(plyranges)
# get eqtls within 500kb
snp.range.eqtl <- data.frame(find_overlaps(gr.snp, gr.cis.eqtl))
# check specific snp
snp.eqtl <- snp.range.eqtl[snp.range.eqtl$SNP==snps[n],]
nrow(snp.eqtl)
# examine snp.range.eqtl
length(unique(snp.range.eqtl$GeneSymbol))
unique(snp.range.eqtl$GeneSymbol)
table(snp.eqtl$GeneSymbol)

snp <- snps[4]
# just get lowest nominal P val for each gene. 
#Genes that mapped to the snp of interest, use those
get.gene <- function(dat, snp){
  out.first <- snp.range.eqtl[snp.range.eqtl$SNP==snp,]
  # remove out first from snp.range.eqtl
  eqtl.rem <- snp.range.eqtl[-(snp.range.eqtl$SNP==snp),]
  
  gene.lst <- unique(eqtl.rem$Gene)
  out.list <- list()
  for (i in seq_along(gene.lst)) {
  sub <-subset(eqtl.rem, eqtl.rem$Gene==gene.lst[i])
  sub <- sub[order(sub$Pvalue),]
  out.list[[i]] <- sub[1,] #add to list
  }
  out <- do.call(rbind, out.list)
  # add out first too
  out <- rbind(out.first, out)
  return(out)
}

test <- get.gene(snp.range.eqtl, snps[1])

######################
# Locuszoom to plot ##
######################
#load gwas data
pheno <- 'del5_7'
gwas <- read.table(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC.meta.gz'),header = T)
names(gwas) <- c('chrom', 'pos', 'rsid', 'other_allele', 'effect_allele', 'N', 'p','p.R','OR','OR.R','Q','I')
gwas <- subset(gwas, gwas$N==max(gwas$N))
#map snp_ids
gwas <- merge(gwas, hrc.rsid, by="rsid", all.x=T)
gwas$rsid <- gwas$ID
head(gwas)

gene <- c('ALDH3B1')

# use a range instead of a gene
snp.range <- gwas[which(gwas$rsid==snp),]$pos
chrom <- gwas[which(gwas$rsid==snp),]$chrom
end <- snp.range + 1e6
start <- snp.range - 1e6

loc2 <- locus(data = gwas, chrom = 'chrom', pos = 'pos', seqname = chrom, xrange = c(start, end), 
              ens_db = "EnsDb.Hsapiens.v75")
loc2 <- link_eqtl(loc2, token = "fc368a631dc0")
loc2 <- link_LD(loc2, token = "fc368a631dc0")
locus_plot(loc2, labels='index')
eqtl_plot(loc2, tissue = "Whole Blood", eqtl_gene = gene)
# intorregate more on loc2 object
table(loc2$LDexp$Gene_Symbol)
t <- loc2$LDexp[loc2$LDexp$Tissue=='Whole Blood',]
table(t$RS_ID)

# plot eqtl
loc2 <- locus(data = gwas, chrom = 'chrom',pos = 'pos', gene = gene, flank = 1e6,ens_db = "EnsDb.Hsapiens.v75")
loc2 <- link_eqtl(loc2, token = "fc368a631dc0")
loc2 <- link_LD(loc2, token = "fc368a631dc0")
locus_plot(loc2, labels='index')
eqtl_plot(loc2, tissue = "Whole Blood", eqtl_gene = gene)
head(loc2$LDexp)

# overlay eqtl plot
overlay_plot(loc2, eqtl_gene = gene,tissue = "Whole Blood")

# intorregate more on loc2 object
table(loc2$LDexp$Gene_Symbol)
t <- loc2$LDexp[loc2$LDexp$Tissue=='Whole Blood',]
table(t$Query)
unique(t$Tissue)

#################
## use of eQTLplot
##################
#deps <- c("biomaRt", "dplyr", "GenomicRanges", "ggnewscale", "ggplot2", "ggplotify", "ggpubr", "gridExtra", "Gviz", "LDheatmap", "patchwork")
#BiocManager::install('snpStats')
# https://github.com/RitchieLab/eQTpLot?tab=readme-ov-file#installation

#devtools::install_github("SFUStatgen/LDheatmap")
#devtools::install_github("RitchieLab/eQTpLot")

library(eQTpLot)
