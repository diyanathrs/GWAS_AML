args=(commandArgs(TRUE))


#phenoprefix <- "status"

dataprefix <- "PCspeAMLHRC"
phenoprefix <- args[1]
phenoprefix
## pngnames <- c('status', 'Normal', 'CBF', 'Trans', 'Complex', 't15_17', 'del5_7', 'Trisomies', 'Any.Monosomy','monosomal.karyotype')
## cat(pngnames, file="phenolist.lst", sep="\n")

## dataprefix <- "AML1_5tidyup"
#dataprefix <- args[1]
dataprefix

lzoutputfolder <- "ForResults"
if (!dir.exists(lzoutputfolder)){
    dir.create(file.path(getwd(), lzoutputfolder))}


x <- c("data.table","metafor")
lapply(x, require, character.only=T)
meta.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, "_hits.lst"), header=TRUE)
hits.snps <- subset(meta.res, N >= 1)$SNP
study.res <- fread(paste0(phenoprefix, "_NCL_", dataprefix, "_hitsByStudy.txt"), header=TRUE)[
    grepl("^NCL1_2", study), study:="NCL1_2"][grepl("^NCL3", study), study:="NCL3"][grepl("^NCL4", study), study:="NCL4"][grepl("^NCL5", study), study:="NCL5"]

study.res[info_score <= 0.6,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

study.res[cases_maf <= 0.02,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

study.res[controls_maf <= 0.02,
                    ':=' (frequentist_add_beta_1=NA, frequentist_add_se_1=NA, frequentist_add_pvalue=NA)
                    ]

study.res[all_BB>=all_AA, ':=' (cases_maf=1-cases_maf, controls_maf=1-controls_maf)]

## check minor allele
## study.res[, minor:=alleleB][all_BB>=all_AA, minor:=alleleA]

## if minor allele == alleleB do nothing
## if minor allele == alleleA give minus
#study.res[minor==alleleA, frequentist_add_beta_1:=frequentist_add_beta_1*-1]
#stopifnot(study.res[is.na(study),.N]==0, unique(study.res, by=c("rsid","minor"))[, .N, by="rsid"][N>1, .N]==0)

setwd("./ForResults")

#par(mar=c(4,4,1,2), mfrow=c(2,2))
for (n in seq_along(hits.snps)){
    plot.sub <- study.res[!is.na(frequentist_add_beta_1) & rsid1==hits.snps[n], ]
    ###if (any(exp(plot.sub[,frequentist_add_beta_1])<=0.1) | any(is.infinite(exp(plot.sub[,frequentist_add_beta_1]))) ){next}
    if (any(is.infinite(exp(plot.sub$frequentist_add_beta_1)))){next}
    pdf(paste0(phenoprefix, hits.snps[n], dataprefix, ".pdf"), onefile=TRUE, width=9, height=6)
    ## plot.sub <- study.res[rsid1==hits.snps[n], ]
    ## if (any(exp(plot.sub[,frequentist_add_beta_1])<=0.1)){next}
    ## if (any(exp(plot.sub[, frequentist_add_beta_1])<=0.1) | any(exp(plot.sub[,frequentist_add_beta_1]) > 50)){next}
    
    fixed <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="FE", slab=plot.sub[, study])
    random <- rma.uni(plot.sub[, frequentist_add_beta_1], sei=plot.sub[, frequentist_add_se_1], method="DL", slab=plot.sub[, study])
                                        #fixedtext1 <- paste("Pooled OR (p",sep="")
    fixedtext2 <- paste("=",format(fixed$QEp,digits=3),";",sep="")
    fixedtext3 <- paste(" (p", sep="")
    randomI <- paste("=",format(random$I2,digits=2),"%)",sep="")
    snpname <- paste0(unique(plot.sub[,rsid]), collapse=",")
    par(mar=c(4,4,1,2))
    i.annot <- plot.sub[, .(format(info_score,digits=2), paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/"))]
    #i.annot <- plot.sub[, paste(format(cases_maf, digits=2),format(controls_maf, digits=2),sep="/")]
    forest(fixed, slab=plot.sub[, study], refline=1, transf=exp, showweight=T,
           ilab=i.annot, ilab.xpos=c(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.14,
                                     par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25),
           mlab=paste0("Fixed-effect (p=", format(fixed$pval,digits=3,scientific = TRUE,nsmall=3), ")"),
           xlab=bquote(paste(.(snpname), .(fixedtext3)[het],.(fixedtext2)," I"^{2},.(randomI),sep="")))
    addpoly(random, row=-0.5, transf=exp,
            mlab=paste0("Random-effect (p=", format(random$pval,digits=3,scientific = TRUE,nsmall=3), ")")
            )
    op <- par(cex=1.2, font=4)
    
    text(par("usr")[1], (par("usr")[4]-par("usr")[3])*0.8 + par("usr")[3],"Study",pos=4,cex=0.8)
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.14, (par("usr")[4]-par("usr")[3])*0.8 + par("usr")[3], "Info_socre",cex=0.8)
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25, (par("usr")[4]-par("usr")[3])*0.8 + par("usr")[3], "EAF \n(Case/Con)",cex=0.8)
    text(par("usr")[1]+4.7*(par("usr")[2]-par("usr")[1])/5.8,  (par("usr")[4]-par("usr")[3])*0.8 + par("usr")[3],  "Weight(%)",pos=2,cex=0.8)
    text(par("usr")[2], (par("usr")[4]-par("usr")[3])*0.8 + par("usr")[3] ,"OR [95% CI]", pos=2,cex=0.8)
    
    ## text(par("usr")[1], (par("usr")[4]-par("usr")[3])/2 + 0.25 ,"Study",pos=4)
    ## text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.25, (par("usr")[4]-par("usr")[3])/2 + 0.3, "MAF\n Case/Con")
    ## text(par("usr")[1]+4.7*(par("usr")[2]-par("usr")[1])/5.8,  (par("usr")[4]-par("usr")[3])/2 + 0.25 ,  "Weight(%)",pos=2)
    ## text(par("usr")[2], (par("usr")[4]-par("usr")[3])/2 + 0.25 ,"OR [95% CI]", pos=2)

    par(op)
    
    op <- par(cex=1.5, font=2)
    text(par("usr")[1]+ (par("usr")[2]-par("usr")[1])*0.4, (par("usr")[4]-par("usr")[3])/2 + 1.25 , paste0(phenoprefix, "_", dataprefix, "\n", hits.snps[n]), pos=4)
    par(op)

dev.off()
}
