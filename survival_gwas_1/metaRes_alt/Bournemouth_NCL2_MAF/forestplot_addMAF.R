args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

#phenoprefix <- "TTFT_DX_RX_Status"

#dataprefix <- "CLL_finalQced4imp"

setwd("./Hitsmeta")

x <- c("data.table","metafor")
lapply(x, require, character.only=T)
meta.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, "_hits.lst"), header=TRUE)
hits.snps <- meta.res[, SNP]
study.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, "_hitsByStudy.txt"), header=TRUE)
#[
#    grepl("^NCL1_2", study), study:="NCL1_2"][grepl("^NCL3", study), study:="NCL3"][grepl("^NCL4", study), study:="NCL4"]
## check minor allele
## study.res[, minor:=alleleB][all_BB>=all_AA, minor:=alleleA]

## if minor allele == alleleB do nothing
## if minor allele == alleleA give minus
## study.res[minor==alleleA, frequentist_add_beta_1:=frequentist_add_beta_1*-1]
## stopifnot(study.res[is.na(study),.N]==0, unique(study.res, by=c("rsid","minor"))[, .N, by="rsid"][N>1, .N]==0)

setnames(study.res, c("beta", "SE"), c("frequentist_add_beta_1", "frequentist_add_se_1"))
study.res[alleleB==minorAllele, ':=' (effAllele=alleleB, refAllele=alleleA)][alleleA==minorAllele, ':=' (effAllele=alleleA, refAllele=alleleB)]

pdf(paste0(phenoprefix, "_NCL_", dataprefix, "_forest_addmafHRC.pdf"), onefile=TRUE, width=12, height=8)
#par(mar=c(4,4,1,2), mfrow=c(2,2))
for (n in seq_along(hits.snps)){
    plot.sub <- study.res[rsid1==hits.snps[n] & !is.na(totN), ]
    if (any(exp(plot.sub[,frequentist_add_beta_1])<=0.1)){next}
    fixed <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
    random <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
                                        #fixedtext1 <- paste("Pooled OR (p",sep="")
    fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
    fixedtext3 <- paste(" (p", sep="")
    randomI <- paste("=",format(random$I2,digits=2),"%)",sep="")
    snpname <- paste0(c(paste0(unique(plot.sub[,rsid]), collapse=","), paste0(unique(plot.sub[,snp]), collapse=",")), collapse=",")
    par(mar=c(4,4,1,2))
    ## i.annot <- plot.sub[, paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/")]
    i.annot <- plot.sub[, .(paste(totN, No.event, sep="/"), paste(effAllele, refAllele, sep="/"), paste0(format(studyMAF, digits=2)))] 
                        #paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/")]
    forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, showweight=T,
           ilab=i.annot,
           ilab.xpos=c(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.18,
                       par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25,
                       par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.32),
           mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE,nsmall=3), ")"),
           xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")))

     ## forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, showweight=F,
     ##       ilab=i.annot, ilab.xpos=c(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25),
     ##       mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE,nsmall=3), ")"),
     ##       xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")))
    addpoly(random, row=-0.5, transf=exp,
            mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")")
            )
    op <- par(cex=1.2, font=4)
    text(par("usr")[1], plot.sub[, .N] + 1.5,"Study", pos=4)
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.18,  plot.sub[, .N] + 1.5, "No/\nevents")
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25,  plot.sub[, .N] + 1.5, "Eff/Ref")
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.32,  plot.sub[, .N] + 1.5, "EAF")
    
    text(par("usr")[1]+ 5.2*(par("usr")[2]-par("usr")[1])/5.8,  plot.sub[, .N] + 1.5 ,  "Weight(%)",pos=2)
    text(par("usr")[2],  plot.sub[, .N] + 1.5 ,"HR [95% CI]", pos=2)
#text((par("usr")[1]+par("usr")[2])/2, (par("usr")[4]-par("usr")[3])/2 + 1.25 , paste0(phenoprefix), pos=4)
    
   ##  text(par("usr")[1], (par("usr")[4]-par("usr")[3])/2 + 0.25 ,"Study",pos=4)

##     text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.2, (par("usr")[4]-par("usr")[3])/2 + 0.3, "No.Subatrisk/No.events")
##     text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25, (par("usr")[4]-par("usr")[3])/2 + 0.3, "EffAllele\RefAllele")
##     text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.23, (par("usr")[4]-par("usr")[3])/2 + 0.3, "EAF")
    
##     text(par("usr")[1]+4.7*(par("usr")[2]-par("usr")[1])/5.8,  (par("usr")[4]-par("usr")[3])/2 + 0.25 ,  "Weight(%)",pos=2)
##     text(par("usr")[2], (par("usr")[4]-par("usr")[3])/2 + 0.25 ,"OR [95% CI]", pos=2)
## #text((par("usr")[1]+par("usr")[2])/2, (par("usr")[4]-par("usr")[3])/2 + 1.25 , paste0(phenoprefix), pos=4)
    par(op)
    
    op <- par(cex=1.5, font=2)
                                        #text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.4, (par("usr")[4]-par("usr")[3])/2 + 1.25 , paste0(phenoprefix, "_", dataprefix ), pos=4)
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.4, plot.sub[, .N] + 2.5, paste0(phenoprefix, "_", dataprefix ), pos=4)
    par(op)
}

dev.off()
