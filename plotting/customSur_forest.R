#plot forest plts
require(remotes)

#install_version("metafor", version = "3.4-0")
x <- c("data.table","metafor")
lapply(x, require, character.only=T)


# snps to plot
snps <- c('rs28645857', 'rs2880743', 'rs115692085', 'rs3791334', 'rs146310501', 'rs184188725', 'rs16829165' )
snps <- 'rs3791334'

####################
### for Sur meta ###
####################
# don't need to redo OS res
#load study level data and save to temp rds
gwas.study.ori <- readRDS('os_assocRes/temp_AML_OSsnps.Rds')
#check if snps are there
snps.2.chk <- setdiff(snps, gwas.study.ori$ID)
if (!length(snps.2.chk) < 1) {
  # load ukb snps
  hrcfilein <- paste0("cut -f1-5 HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
  hrc.rsid <- fread(hrcfilein, header=TRUE)
  hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
  head(hrc.rsid)
  # only keep snps on int
  hrc.rsid <- hrc.rsid[hrc.rsid$ID %in% c(snps),]
  head(hrc.rsid)
  snps.2.chk <- hrc.rsid[hrc.rsid$ID %in% snps.2.chk]$rsid
  
  study.lst <- list.files('os_assocRes/', pattern = '*gz')
  gwas.study <- parallel::mclapply(study.lst, function(x) {
    fread(paste0('zcat os_assocRes/', x), header = T) }, mc.cores = 8)
  #change study names
  study.lst <- sub('assoc_HRC.gz','',study.lst)
  study.lst <- sub('BirminghamNCL5|BirNCL5','UK3',study.lst)
  study.lst <- study.lst %>% sub('_[^_]*_','_', .) %>% sub('_[^_]*_','_',.)
  names(gwas.study) <- sub('assoc_HRC.gz','',study.lst)
  gwas.study <- rbindlist(gwas.study, idcol = "source")
  gwas.study$snp <- paste0(gwas.study$chromosome,':',gwas.study$position,'_',gwas.study$alleleA,'_',gwas.study$alleleB)
  #keep only the snps of interest
  gwas.study <- gwas.study[gwas.study$snp %in% snps.2.chk]
  gc()
  #merge rsids again as us1 us2 are missing rsids
  gwas.study <- merge(gwas.study, hrc.rsid[,c(3,6)], by.x='snp', by.y='rsid')
  table(gwas.study$source)
  #split col
  gwas.study <- separate(gwas.study, source, into = c('cohort', 'pheno') ,sep='_')
  test <- rbind(gwas.study.ori, gwas.study)
  #save.table
  saveRDS(test, paste0('os_assocRes/','temp_AML_OSsnps.Rds')) 
}

gwas.study <- setDT(readRDS('os_assocRes/temp_AML_OSsnps.Rds'))
#convert EAF to 2 digits
gwas.study$studyMAF <- as.numeric(format(round(gwas.study$studyMAF, digits=2)))
rm(gwas.study.ori)

names(gwas.study)
setnames(gwas.study, c('beta', 'SE'), c('frequentist_add_beta_1','frequentist_add_se_1'))
stdy.odr <- c("UK1","UK2","UK3", "Germany", "Hungary", "Finland", "US1", "US2")
gwas.study$cohort <- factor(gwas.study$cohort, levels = stdy.odr)
gwas.study <- gwas.study[order(gwas.study$cohort)]
#check
snps %in% gwas.study$ID
dir.create('custom_forestPlts2')

i <- 7
for (i in seq_along(snps)) {
  print(snps[i])
  plot.snp <- gwas.study[gwas.study$ID==snps[i]]
  for (subtype in unique(plot.snp$pheno)) {
    print(subtype)
    plot.sub <- plot.snp[plot.snp$pheno==subtype]
    # change direction according to HR
    mean.hr <- mean(as.numeric(plot.sub$crudeHR), na.rm=T)
    if (mean.hr < 1) {
  # change directon <<<< remove this
  plot.sub$frequentist_add_beta_1 <- -(plot.sub$frequentist_add_beta_1)
  plot.sub$studyMAF <- 1-plot.sub$studyMAF
  #flip alllele
  setnames(x = plot.sub, c('alleleA','alleleB'), c('alleleB','alleleA'))
    }
  # round info score
  if (max(plot.sub$info_score) >= 0.995){
    plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
  } else {
    plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
  }
  range(plot.sub$info_score.2)
  plot.sub$cohort
  #plot.sub <- plot.sub[c(1,2,8,3:7)
  head(plot.sub)
  
  fixed <- rma.uni(yi = plot.sub$frequentist_add_beta_1, sei=plot.sub$frequentist_add_se_1, method="FE", slab=plot.sub$cohort)
  random <- rma.uni(yi = plot.sub$frequentist_add_beta_1, sei=plot.sub$frequentist_add_se_1, method="DL", slab=plot.sub$cohort)
  fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
  fixedtext3 <- paste(" (p", sep="")
  #fixedtext3 <- bquote(paste(" (p",.[het])) #test
  randomI <- paste("=",format(random$I2,digits=2),"%)",sep="")
  snpname <- paste0(unique(plot.sub$ID), collapse=",")
  i.annot <- plot.sub[, .(paste(No.event, totN ,sep="/"),
                          ifelse(Genotyped=="Genotyped", "Genotyped", format(plot.sub$info_score.2,digits=2)), 
                          paste(alleleB, alleleA ,sep="/"), format(studyMAF, digits=1))]
  i.annot <- plot.sub[, .(paste(No.event, totN ,sep="/"),
                          ifelse(Genotyped=="Genotyped", "Genotyped", format(plot.sub$info_score.2,digits=2)), 
                          paste(alleleB, alleleA ,sep="/"), studyMAF)]
  
  #i.annot <- plot.sub[, paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/")]
  pdf(paste0('custom_forestPlts2/',snps[i],'_', subtype,"_forestAMLHRC.pdf"), onefile=F, width=8, height=4.8)
  par(mar=c(4, 4, 1, 2)) # use 4422 for smaller plots

  x.range <- c(min(as.numeric(plot.sub$crudeHR_L95))-18, max(as.numeric(plot.sub$crudeHR_H95))+9)
  x.range <- c(-16, 11) #for small use -17, 12
  forest(fixed, slab=plot.sub[, cohort], refline=1, transf=exp,
         xlim=x.range,at=c(0, 1, 3),
         ilab=i.annot,
         ilab.xpos=c(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.18,
                     par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.32,
                     par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.44,
                     par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.54),
         mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
         xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
         cex = 1, colout = 'black', pch = 19, lwd=1.7, col = 'black', border = 'black', showweights = T)
  
  addpoly(random, row=-0.4, transf=exp, cex = 1, col = 'black', border = 'black',
          mlab=paste0("Random-effect (p=", format(random$pval,digits=3, scientific = TRUE,nsmall=3), ")"))
  pos <- c(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.18,
           par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.32,
           par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.44,
           par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.54)
  op <- par(cex=1, font=2)
  text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.85 + par("usr")[3],"GWAS", pos=4)
  text(pos, 9.3, c('No of \n events',"Info \n score","EA/RA","EAF" ), adj = c(0.5,0))
  #text(2.5, 9.3, "Weight", adj = c(0,0))
  text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.85 + par("usr")[3] ,"Weight        HR [95% CI]", pos=2)
  dev.off()
  } 
  }

  
  