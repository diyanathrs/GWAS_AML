#plot forest plts
require(remotes)
install_version("metafor", version = "3.4-0")
x <- c("data.table","metafor")
lapply(x, require, character.only=T)

# load ukb snps
hrcfilein <- paste0("cut -f1-5 HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
hrc.rsid <- fread(hrcfilein, header=TRUE)
hrc.rsid[, rsid:=paste0(`#CHROM`, ":", POS, "_", REF, "_", ALT)]
head(hrc.rsid)

# snps to plot
#rs11481 (11q13.2) b rs4665765 (2p23.3) c rs3916765 (6p21.32) d rs12078864 (1q23.2) 
#rs79918355 (2p21) f rs12988876 (2q33.3) g rs12632224 (3q25.33) h rs11212666 (11q23) 
#rs2853677 (5p15.3) J rs7705526 (5p15.3). 

## new snps from CH - c("rs12632224" , "rs11212666" , "rs2853677" , "rs7705526" , "rs13130545" , "rs2086132" , "rs8088824", "rs35452836" , "rs79633204" , "rs10131341")
snps <- c("rs11481", "rs4665765"  ,"rs3916765" ,"rs12078864" ,"rs79918355" , "rs12988876"  ,"rs12632224" ,
          "rs11212666" ,"rs2853677" , "rs7705526")

# only keep snps on int
hrc.rsid <- hrc.rsid[hrc.rsid$ID %in% snps,]

# manual input of case/con numbers
#pheno <- c('status', 'normal', 'complex', 'del57')
study <- paste('GWAS', seq_along(1:6))
all.cases <- c(1119, 931, 991, 977, 351 ,341)
normal <- c(387, 177, 286, 465, 128, 137)
complx <- c(75, 61, 89, 82 ,11, 51)
del5_7 <- c(72, 58, 68, 57, 17, 64)
con <- c(2671, 2477, 1612, 3728, 1055, 1395)

case.con <- data.frame(study, all.cases,con, normal, complx, del5_7)
#######################
## for etiological meta
########################
# load from meta results
# change names
pheno <- 'status'
study.res <- fread(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC_hitsByStudy.txt'), header=TRUE)[
  grepl("^NCL1_2", study), study:="GWAS 1"][grepl("^NCL3", study), study:="GWAS 2"][
    grepl("^NCL4", study), study:="GWAS 3"][grepl("^NCL5", study), study:="GWAS 4"][
      grepl("^NCL6", study), study:="GWAS 5"][grepl("^NCL7", study), study:="GWAS 6"]

#map rsids again
study.res <- merge(study.res, hrc.rsid[,c(3,6)], by.x="rsid1", by.y='rsid')
study.res[all_BB>=all_AA, ':=' (cases_maf=1-cases_maf, controls_maf=1-controls_maf)]
# format info score
#study.res$info_score.2 <- format(floor(study.res$info_score*100)/100, nsmall=2)

#################
## hit 1 of status
################
n <- 1

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 
# for normal flip the allele and OR
#-(plot.sub$frequentist_add_beta_1)

fixed <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

#i.annot <- plot.sub[, .(format(info_score,digits=2), alleleA ,paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/"))]
i.annot <- plot.sub[, .(paste(all.cases, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleB,
                        paste(format(cases_maf, digits=2), format(controls_maf, digits=2),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6)
par(mar=c(4, 4, 1, 2)) # use 4422 for smaller plots
pos <- c(-3.3, -2.1, -1, 0)
#pos <- c(-4, -3, -1.5, -0.2) # for status and normal
# pos <- c(-4, -2.6, -1.5, -0.2) pos for del5_7

forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-5, 3.5),
       at=c(0.5, 1, 2),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos, #header = c("GWAS", "OR [95% CI]"), 
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3],"GWAS", pos=4)
text(pos, 7.8, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ),adj = c(0.5,0.5))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()

#################
## hit 2 of status
################
n <- 2

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 
# for normal flip the allele and OR
#-(plot.sub$frequentist_add_beta_1)

fixed <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

#i.annot <- plot.sub[, .(format(info_score,digits=2), alleleA ,paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/"))]
i.annot <- plot.sub[, .(paste(all.cases, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleB,
                        paste(format(cases_maf, digits=2), format(controls_maf, digits=2),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6)
par(mar=c(4, 4, 1, 2)) # use 4422 for smaller plots
pos <- c(-3.3, -2.1, -1, 0)
#pos <- c(-4, -3, -1.5, -0.2) # for status and normal
# pos <- c(-4, -2.6, -1.5, -0.2) pos for del5_7

forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-5, 3.5),
       at=c(0.5, 1, 2),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos,
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3],"GWAS", pos=4)
text(pos, 7.8, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ),adj = c(0.5,0.5))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()

###########
#normal only - hit3
############
# change names
pheno <- 'Normal'
study.res <- fread(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC_hitsByStudy.txt'), header=TRUE)[
  grepl("^NCL1_2", study), study:="GWAS 1"][grepl("^NCL3", study), study:="GWAS 2"][
    grepl("^NCL4", study), study:="GWAS 3"][grepl("^NCL5", study), study:="GWAS 4"][
      grepl("^NCL6", study), study:="GWAS 5"][grepl("^NCL7", study), study:="GWAS 6"]

#map rsids again
study.res <- merge(study.res, hrc.rsid[,c(3,6)], by.x="rsid1", by.y='rsid')
# format info score
study.res$info_score.2 <- format(floor(study.res$info_score*100)/100, nsmall=2)
# for del57
#study.res$info_score.2 <- format(round(study.res$info_score, digits = 2), nsmall=2)

n <- 3

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 

fixed <- rma.uni(plot.sub[, -(frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, -(frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

#i.annot <- plot.sub[, .(format(info_score,digits=2), alleleA ,paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/"))]
i.annot <- plot.sub[, .(paste(normal, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleA,
                        paste(format(1-cases_maf, digits=2), format(1-controls_maf, digits=2),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6, onefile = F)
par(mar=c(4,4,1,2))
pos <- c(-4.9, -3.1, -1.6, -0.1)
forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-7.5, 5.5),
       at=c(0.5,1,2, 3),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos,
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3],"GWAS", pos=4)
text(pos,7.8, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.82 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()

#################
## hit 4 of del5_7
################
pheno <- 'del5_7'
study.res <- fread(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC_hitsByStudy.txt'), header=TRUE)[
  grepl("^NCL1_2", study), study:="GWAS 1"][grepl("^NCL3", study), study:="GWAS 2"][
    grepl("^NCL4", study), study:="GWAS 3"][grepl("^NCL5", study), study:="GWAS 4"][
      grepl("^NCL6", study), study:="GWAS 5"][grepl("^NCL7", study), study:="GWAS 6"]

#map rsids again
study.res <- merge(study.res, hrc.rsid[,c(3,6)], by.x="rsid1", by.y='rsid')
# format info score
#study.res$info_score.2 <- format(floor(study.res$info_score*100)/100, nsmall=2)
#plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
# for del57
#study.res$info_score.2 <- format(round(study.res$info_score, digits = 2), nsmall=2)
n <- 4

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 
# for normal flip the allele and OR
#-(plot.sub$frequentist_add_beta_1)

fixed <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

i.annot <- plot.sub[, .(paste(del5_7, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleB,
                        paste(format(cases_maf, digits=2), format(controls_maf, digits=2),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6)
par(mar=c(4, 4, 2, 2)) # use 4422 for smaller plots
pos <- c(-4, -2.6, -1.5, -0.2) 
# pos <- c(-4, -2.6, -1.5, -0.2) pos for del5_7

forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-6, 5),
       at=c(0.5, 1, 2.5),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos,
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3],"GWAS", pos=4)
text(pos, 6.75, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ),adj = c(0.5,0.5))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()



#################
## hit 5 of complex
################
pheno <- 'Complex'
study.res <- fread(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC_hitsByStudy.txt'), header=TRUE)[
  grepl("^NCL1_2", study), study:="GWAS 1"][grepl("^NCL3", study), study:="GWAS 2"][
    grepl("^NCL4", study), study:="GWAS 3"][grepl("^NCL5", study), study:="GWAS 4"][
      grepl("^NCL6", study), study:="GWAS 5"][grepl("^NCL7", study), study:="GWAS 6"]

#map rsids again
study.res <- merge(study.res, hrc.rsid[,c(3,6)], by.x="rsid1", by.y='rsid')
n <- 5

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
#plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 
# for normal flip the allele and OR
#-(plot.sub$frequentist_add_beta_1)

fixed <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

i.annot <- plot.sub[, .(paste(complx, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleB,
                        paste(format(cases_maf, digits=1), format(controls_maf, digits=1),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6)
par(mar=c(4, 4, 2, 2)) # use 4422 for smaller plots
pos <- c(-4.9, -3.2, -1.9, -0.2) 

forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-7.5, 7.5),
       at=c(0.8, 3, 4.5),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos,
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3],"GWAS", pos=4)
text(pos, 6.75, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ),adj = c(0.5,0.5))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()

#################
## hit 6 of complex
################
pheno <- 'Complex'
study.res <- fread(paste0('../AMLmeta_results/',pheno,'_NCL_PCspeAMLHRC_hitsByStudy.txt'), header=TRUE)[
  grepl("^NCL1_2", study), study:="GWAS 1"][grepl("^NCL3", study), study:="GWAS 2"][
    grepl("^NCL4", study), study:="GWAS 3"][grepl("^NCL5", study), study:="GWAS 4"][
      grepl("^NCL6", study), study:="GWAS 5"][grepl("^NCL7", study), study:="GWAS 6"]

#map rsids again
study.res <- merge(study.res, hrc.rsid[,c(3,6)], by.x="rsid1", by.y='rsid')
n <- 6

plot.sub <- study.res[!is.na(frequentist_add_beta_1) & ID==snps[n], ]
#plot.sub$info_score.2 <- format(floor(plot.sub$info_score*100)/100, nsmall=2)
plot.sub$info_score.2 <- format(round(plot.sub$info_score, digits = 2), nsmall=2)
plot.sub$alleleA
# add case con numbers
plot.sub <- merge(plot.sub, case.con, by='study') 
# for normal flip the allele and OR
#-(plot.sub$frequentist_add_beta_1)

fixed <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
random <- rma.uni(plot.sub[, (frequentist_add_beta_1)], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
#fixedtext1 <- paste("Pooled OR (p",sep="")
fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
fixedtext3 <- paste(" (p", sep="")
randomI <- paste("=",format(random$I2, digits=2),"%)",sep="")
#snpname <- paste(unique(plot.sub[,ID]), '6:32685550_A_G', sep =";")
snpname <- paste(unique(plot.sub[,ID]), unique(plot.sub[,rsid1]), sep =";")

i.annot <- plot.sub[, .(paste(complx, con ,sep="/"),
                        ifelse(Genotyped=="Genotyped", "Genotyped", info_score.2),alleleB,
                        paste(format(cases_maf, digits=1), format(controls_maf, digits=1),sep="/"))]

pdf(file=paste0(snps[n], pheno,"_Figure2_2_main.pdf"), width=8, height=4.6)
par(mar=c(4, 4, 2, 2)) # use 4422 for smaller plots
pos <- c(-4.9, -3.2, -1.9, -0.2) 

forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, 
       xlim=c(-7.5, 7.4),
       at=c(0.8, 3, 4.5),
       #clim=c(0.5, 3.5),
       ilab=i.annot,
       ilab.xpos=pos,
       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE, nsmall=3), ")"),
       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")),
       cex = 1.1, lwd=1.5, pch=19)

addpoly(random, row=-0.5, transf=exp, cex = 1.1, lwd=1.5,
        mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")"))
op <- par(cex=1.1, font=2)
#par(op)
text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3],"GWAS", pos=4)
text(pos, 6.75, c('No \n case/con',"Info \n score","Effect \n allele","EAF \n case/con" ),adj = c(0.5,0.5))
text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.81 + par("usr")[3] ,"OR [95% CI]", pos=2)
dev.off()

