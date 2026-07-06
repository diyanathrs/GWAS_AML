#plot forest plts
require(remotes)

#install_version("metafor", version = "3.4-0")
x <- c("data.table","metafor")
lapply(x, require, character.only=T)

# load ukb snps
hrcfilein <- paste0("cut -f1-5 HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# snps to plot
#ch.snps <- c("rs12632224", "rs11212666", "rs2853677", "rs7705526", "rs10131341")
snps <- c("rs11481", "rs4665765"  ,"rs3916765" ,"rs12078864" ,"rs79918355" , "rs12988876"  ,"rs12632224" ,
          "rs11212666" ,"rs2853677" , "rs7705526")

snps <- 'rs6450183'
# only keep snps on int
hrc.rsid <- hrc.rsid[hrc.rsid$ID %in% c(snps),]
head(hrc.rsid)

####################
### for Sur meta ###
####################
# don't need to redo OS res
redo_res <- T
if (isTRUE(redo_res)) {
#load study level data
study.lst <- list.files('OSassocRes_normal/', pattern = '*gz')
gwas.names <- c('UK3', 'Finland', 'Germany', 'Hungary', 'UK1', 'UK2', 'US1', 'US2')

gwas.study <- parallel::mclapply(study.lst, function(x) {
  fread(paste0('zcat OSassocRes_normal/', x), header = T)
}, mc.cores = 8)

names(gwas.study) <- gwas.names
gwas.study <- rbindlist(gwas.study, idcol = "source")
gwas.study$snp <- paste0(gwas.study$chromosome,':',gwas.study$position,'_',gwas.study$alleleA,'_',gwas.study$alleleB)
head(gwas.study)
head(hrc.rsid)
#keep only the snps of interest
gwas.study <- gwas.study[gwas.study$snp %in% hrc.rsid$rsid]
gc()
#merge rsids again as us1 us2 are missing rsids
gwas.study <- merge(gwas.study, hrc.rsid[,c(3,6)], by.x='snp', by.y='rsid')
table(gwas.study$source)
#save.table
saveRDS(gwas.study,'AML_OSsnps.Rds') }
#################
## start here ###
#################
gwas.study <- readRDS('AML_OSsnps.Rds')
names(gwas.study)
setnames(gwas.study, c('beta', 'SE'), c('frequentist_add_beta_1','frequentist_add_se_1'))
stdy.odr <- c("UK1","UK2","UK3", "Germany", "Hungary", "Finland", "US1", "US2")
gwas.study$source <- factor(gwas.study$source, levels = stdy.odr)
gwas.study <- gwas.study[order(gwas.study$source)]
#check
snps %in% gwas.study$ID

## plot forest
#study.res[, ':=' (effAllele=alleleB, refAllele=alleleA)]
#study.res[alleleA == minorAllele, frequentist_add_beta_1 := -1*frequentist_add_beta_1]
#import collider correction
collider <- readRDS('../AMLsur/collider_bias/AML_collider_corrected.Rds')
head(collider)

i <- 'rs3916765'
#snps <- ch.snps
for (i in ch.snps) {
print(i)
plot.sub <- gwas.study[gwas.study$ID==i]
# round info score
if (max(plot.sub$info_score) >= 0.995){
  plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
} else {
  plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
}
range(plot.sub$info_score.2)
plot.sub$source
#plot.sub <- plot.sub[c(1,2,8,3:7)
head(plot.sub)
snp.collider <- subset(collider, collider$SNP==unique(plot.sub$snp))

fixed <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, source])
random <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, source])
collider.corr <- rma.uni(snp.collider[, "beta.adj"], sei=snp.collider[, "se.adj"], method="FE")
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
#fixedtext3 <- bquote(paste(" (p",.[het])) #test
randomI <- paste("=",format(random$I2,digits=2),"%)",sep="")
snpname <- paste0(unique(plot.sub[,ID]), collapse=",")
#i.annot <- plot.sub[, .(format(info_score,digits=2), alleleB ,paste(format(studyMAF, digits=2),sep="/"))]
#i.annot <- plot.sub[, .(paste(totN, No.event, sep="/"), paste(alleleB), paste0(format(studyMAF, digits=2)))]
i.annot <- plot.sub[, .(paste(totN, No.event ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", format(info_score.2,digits=2)),alleleB,
                        paste(format(studyMAF, digits=2)))]

#i.annot <- plot.sub[, paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/")]
pdf(paste0('os_forestNew/',snpname, "_forest_addmafHRC.pdf"), onefile=F, width=9, height=5.2)
par(mar=c(4,4,2,2))
pos <- c(-5.5, -3.5, -1.8, -0.5)

forest(fixed, slab=plot.sub[, source], refline=1, transf=exp, showweight=T, xlim = c(-8, 6), alim = c(0,3), at = c(0,1,2),
       ilab=i.annot, ilab.xpos = pos, header = c("GWAS", "Weight       HR [95% CI]"),
       mlab=paste0("Fixed-effect (p=", format(fixed$pval, digits=3, scientific = TRUE,nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")), cex = 0.9)
       #xlab=bquote(paste(.(fixedtext3),.(fixedtext2)," I"^{2})))

addpoly(random, row=-0.5, transf=exp, mlab=paste0("Random-effect (p=", format(random$pval, digits=3,scientific = TRUE,nsmall=3), ")"))
addpoly(collider.corr, row=-1.5, transf=exp, mlab=paste0("Adjusted for index bias (p=", format(snp.collider$p.adj, digits=3,scientific = TRUE,nsmall=3), ")"))
text(pos,10, c('No case/\n events','Info \n score' ,'Effect \nallele','EAF'), font=2, cex=0.9)
#text(par("usr")[1]+ 4.6*(par("usr")[2]-par("usr")[1])/5.8,  plot.sub[, .N] + 1.5 ,  "Weight",pos=2)
dev.off()
}
