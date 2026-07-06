###############################
# create fastenloc input file from susie
############################
susie.lst <- list.files(path = '../polyfun/new',pattern = 'susie.out', full.names = T)
susie.lst
n <- 5
susie <- read.table(susie.lst[n],header = T)
head(susie)
# output needs a tab file Locus_id   SNP_id   beta_gwas   se_gwas
#get SE
susie$SE <- susie$BETA_MEAN/susie$Z
susie$locus_id <- 'panAML_sig1'
names(susie)
write.table(susie[c(16,2,11,15)], file = paste0('fastnloc_out/fastnlocIn_','panAML_sig1'), quote = F, sep = '\t', row.names = F)

# now run system - 
#system(paste0('susie2enloc -dir susie_rst_dir -vcf snp_vcf_file [-tissue tissue_name] | gzip - > fastenloc.susie.annotation.vcf.gz'))

# upper isnt gonna work as it needs susie raw output files additional to above files

####################
##using sum stat as input - hybrid mode
#####################
# map rsids
file <- '../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab'
hrcfilein <- paste0("cut -f1-5 ../publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

#format - Locus_id   SNP_id   beta_gwas   se_gwas
# locus id can be status-hit1? 
# take each sumstat and calculate beta and SE -log(ORs)
sumstat <- read.table('../AMLmeta_results/status_NCL_PCspeAMLHRC.meta.gz', header = T)
sumstat <- subset(sumstat, sumstat$N==(max(sumstat$N)))
head(sumstat)
# map rsid
sumstat <- merge(sumstat, hrc.rsid, by.x="SNP",by.y='rsid', all.x=T)
#calculate beta and SE
sumstat$beta <- log(sumstat$OR)
sumstat$Z <- qnorm(1 - sumstat$P/2)
sumstat$SE <- sumstat$beta/sumstat$Z
head(sumstat)
# separate signals around lead SNPs
lead.snps <- c('rs4665765', 'rs11481')
n <- 2
chrom <- sumstat[which(sumstat$ID==lead.snps[n]),]$CHR
pos <- sumstat[which(sumstat$ID==lead.snps[n]),]$BP
sum.sub <- subset(sumstat, sumstat$CHR==chrom)
sum.sub <- sum.sub[order(sum.sub$BP),]
# subset 0.5mb from both sides
sum.sub <- dplyr::filter(sum.sub, between(BP, pos-1e4, pos+1e4))
# round beta and se
sum.sub$beta <- round(sum.sub$beta,digits = 4)
#sum.sub$SE <- round(sum.sub$SE,digits = 3)

# use only 4 req cols
name <- paste0('panAML_',lead.snps[n])
sumstat$Locus_id <- name
names(sum.sub)
write.table(sum.sub[c(21,15,18,20)], file = paste0('fastnloc_out/fastnlocIn_',name), quote = F, sep = '\t', row.names = F, col.names = F)

#################
# check with other data
stat1.z <- read.table('fastnloc_out/test/Regionstatus1.z', header = T)
head(stat1.z)
stat1.z <- merge(stat1.z, hrc.rsid, by='rsid', all.x=T)
stat1.z <- subset(stat1.z, stat1.z$N==(max(stat1.z$N)))
name <- paste0('panAML1')
stat1.z$Locus_id <- name
names(stat1.z)
write.table(stat1.z[c(20,17,7,8)], file = paste0('fastnloc_out/fastnlocIn_test',name), quote = F, sep = '\t', row.names = F)

########################
#### get prob input vcf manually
########################
# chr1	115746	chr1_115746_C_T_b38	C	T	ENSG00000269981:1@Spleen=2.00812e-01[9.997e-01:4]
# The first five columns represent the chromosome, position, SNP ID, reference allele, and alternative allele, consistent with a standard VCF file. 
# Note that FastENLOC uses only the SNP ID information and does not verify or utilize the position or allele information. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

library(VariantAnnotation)

susie.lst <- list.files(path = '../polyfun/new',pattern = 'susie.out', full.names = T)
susie.lst
n <- 5
susie <- read.table(susie.lst[n],header = T)
head(susie)

# Create a GRanges object
susie.gr <- GRanges(seqnames=susie$CHR, ranges=IRanges(start=susie$BP, width=1))
mcols(susie.gr)$ID <- susie$SNP
head(susie.gr) 

# Build fixed and info columns
#fixed <- DataFrame(REF=DNAStringSet(ref), ALT=DNAStringSetList(strsplit(alt, ",")), QUAL=NA_real_, FILTER="PASS")
#locus_id@tissue=pip[cpip:number_of_snps] - create df for this
top.hit <- order(susie$PIP, decreasing = T)[1]
top.hit <- susie$SNP[top.hit]
susie$CREDIBLE_SET <- susie$CREDIBLE_SET+1
no.cred <- unique(susie$CREDIBLE_SET)

cred <- list()
for (i in seq_along(unique(susie$CREDIBLE_SET))) {
  cred$sum[i] <- sum(subset(susie, susie$CREDIBLE_SET==i)$PIP)
  cred$snps[i] <- nrow(subset(susie, susie$CREDIBLE_SET==i))
}

#top.hit,'@','status','=',
info <- DataFrame(paste0(susie$PIP,'[',cred$sum[susie$CREDIBLE_SET],':',
                     cred$snps[susie$CREDIBLE_SET],']'))
names(info) <- paste0(top.hit,'@','status')
row.names(info) <- susie$SNP


fixed <- DataFrame(REF=DNAStringSet(susie$A2), ALT=DNAStringSetList(strsplit(susie$A1, ",")))

# Create VCF object
vcf <- VCF(rowRanges=susie.gr, info = info, fixed = fixed)


# Define metadata for INFO field (optional but recommended)
info_header <- DataFrame(
  Number = "1",
  Type = "Character",
  Description = "rsid",
  row.names = "ID"
)

header <- VCFHeader(
  reference = "",
  header = DataFrameList(INFO = info_header)
)

# Assign header to VCF object
metadata(vcf)$header <- header
# write
writeVcf(vcf, filename="AMLstatus1.vcf")

#################
### second method
#################
install.packages("vcfR")
library(vcfR)

### create matrix from data
head(susie)
susie.mat <- as.matrix(susie)

vcf_data <- matrix(".", nrow=nrow(susie), ncol=10)  # fill with real genotype data

colnames(vcf_data) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","sample1")

vcf_data[1,] <- c("chr1", "123456", "rs1", "A", "G", ".", "PASS", ".", "GT", "0/1")
vcf_data[2,] <- c("chr1", "123789", "rs2", "G", "A", ".", "PASS", ".", "GT", "1/1")

vcf_obj <- new("vcfR", meta = c("##fileformat=VCFv4.2"), fix = vcf_data[,1:8], gt = vcf_data[,9:10])

# save
write.vcf(vcf_obj, file="output.vcf")
