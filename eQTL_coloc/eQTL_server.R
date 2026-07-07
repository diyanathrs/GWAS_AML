#################
# Get cis-eQTL ##
#################

####. your inputs below  #########################################

snps <- c("rs11481" , "rs2164808") # <<<-- Add your SNP rsid here
distance  <- 250000  # <<<-- Distance from the SNP to look eQTL for

##################################################################
# Run the script from here
# load libraries
require(data.table)
setwd('.')

# load ukb snps
hrcfilein <- paste0("cut -f1-5 ~/AML/publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE, )
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]

geteQTL <- function(snp){
  #get snp pos
  snp.of.int <- hrc.rsid[which(hrc.rsid$ID %in% snp)]
  interval <- paste0(snp.of.int$`#CHROM`,':',snp.of.int$POS - distance,'-',snp.of.int$POS + distance)
  
  eqtl.in <- "cis-eQTLs_full_20180905tabix.gz"
  out <- tempfile()
  tabixcmd <- paste0(" tabix -hf ", file.path(eqtl.in), " ", interval, " > ", out)
  try(system(tabixcmd))
  if (file.info(out)$size==0)stop("No eQTL data found for the region..")
  dsn1 <- fread(out, header = TRUE)
  unlink(out)
  
  # filter for snp of interest
  eQTL.out <- subset(dsn1, dsn1$SNP==snp)
  stopifnot("No eQTL for the SNP.."=nrow(eQTL.out) > 0)
  print(paste("Found",nrow(eQTL.out),"eQTL associations for",snp))
  setnames(eQTL.out, "#Pvalue", "Nominal_Pval")
  eQTL.out <- eQTL.out[,c(1,14,8,9,11,2,3,4,5,6,7,12,13)]
  # sort
  eQTL.out <- eQTL.out[order(eQTL.out$Nominal_Pval,decreasing = F)]
  file.nm <- paste0('eQTL_', snp ,'.csv')
  print(paste("Saving the output as",file.nm , "in eQTL_out folder.."))
  write.csv(eQTL.out, file.path('eQTL_out',file.nm), quote = F, row.names = F)
}

for (snp in snps){
  print(paste("Getting eQTL for", snp))
  geteQTL(snp)
}
## done ##


####. your inputs below  #########################################

snps <- c("rs11481" , "rs2164808") 
distance  <- 250000 

# Run the script from here
# load libraries
require(data.table)
setwd('.')

# load ukb snps
hrcfilein <- paste0("cut -f1-5 ~/AML/publication_plots/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]

geteQTL <- function(snp){
  #get snp pos
  snp.of.int <- hrc.rsid[which(hrc.rsid$ID %in% snp)]
  interval <- paste0(snp.of.int$`#CHROM`,':',snp.of.int$POS - distance,'-',snp.of.int$POS + distance)
  
  eqtl.in <- "cis-eQTLs_full_20180905tabix.gz"
  out <- tempfile()
  tabixcmd <- paste0(" tabix -hf ", file.path(eqtl.in), " ", interval, " > ", out)
  try(system(tabixcmd))
  if (file.info(out)$size==0)stop("No eQTL data found for the region..")
  dsn1 <- fread(out, header = TRUE)
  unlink(out)
  
  # filter for snp of interest
  eQTL.out <- subset(dsn1, dsn1$SNP==snp)
  stopifnot("No eQTL for the SNP.."=nrow(eQTL.out) > 0)
  print(paste("Found",nrow(eQTL.out),"eQTL associations for",snp))
  setnames(eQTL.out, "#Pvalue", "Nominal_Pval")
  eQTL.out <- eQTL.out[,c(1,14,8,9,11,2,3,4,5,6,7,12,13)]
  # sort
  eQTL.out <- eQTL.out[order(eQTL.out$Nominal_Pval,decreasing = F)]
  file.nm <- paste0('eQTL_', snp ,'.csv')
  print(paste("Saving the output as",file.nm , "in eQTL_out folder.."))
  write.csv(eQTL.out, file.path('eQTL_out',file.nm), quote = F, row.names = F)
}

for (snp in snps){
  print(paste("Getting eQTL for", snp))
  geteQTL(snp)
}
