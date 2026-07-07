##source("/home/weiyu/AML/AML_HRCpcspecRocket151118/shinyDevel/AMLQ7_8/jimfun.R")
source("/home/weiyu/AML/AML_HRCpcspecRocketNatComRevTrim/ManuscriptFigures/jimfunR1.R")
retired <- data.table::fread("zcat /home/weiyu/UCSC_SNPs/RsMergeArch.bcp.gz")


pheno <- "status"
dataver <-  "HRC"
snpid <- "rs4665765"

pheno <- "status"
dataver <-  "HRC"
snpid <- "rs11481"

# Hollys snps 
pheno <- "status"
dataver <-  "HRC"
snpid <- "rs16829165"

#pheno <- "del5_7l"
#dataver <-  "HRC"
#snpid <- "rs12078864"



snp_pos <- snpformatcheck(x=snpid)
snp_chr <- sapply(strsplit(snp_pos, split=":"), "[[", 1)
snp_pos <- as.integer(unlist(strsplit(sapply(strsplit(snp_pos, split=":"), "[[", 2), "-"))[1])

check_interval <- paste0(snp_chr, ":", max(1, snp_pos - 250000), "-", snp_pos + 250000)



geteQTL <- function(x="blood", interval){
    require(data.table)
    eqtl.dir <- ifelse(x=="eQTLgen", "/home/weiyu/AML/AML_HRCpcspecRocket151118/eQTL",
                       "/home/weiyu/AML/AML_RocketHPC/eQTL"
                       )
    eqtl.res <-c("WholeBloodGTExV7_eQTL_V7.header.gz",
                 "LymphocytesGTExV7_eQTL_V7.header.gz",
                 "cis-eQTLs_full_20180905tabix.gz")
    ## eqtlfile <- ifelse(x=="blood", "WholeBloodGTExV7_eQTL_V7.header.gz",
    ##                    "LymphocytesGTExV7_eQTL_V7.header.gz")

     eqtlfile <- switch(x,
                       "blood"="WholeBloodGTExV7_eQTL_V7.header.gz",
                       "Lymphocyte"="LymphocytesGTExV7_eQTL_V7.header.gz",
                       "eQTLgen"="cis-eQTLs_full_20180905tabix.gz"
                       )

    
    out <- tempfile()
    tabixcmd <- paste0(" tabix -hf ", file.path(eqtl.dir, eqtlfile), " ", interval, " > ", out)
    try(system(tabixcmd))
    if (file.info(out)$size==0)stop("no data found")
    dsn1 <- fread(out, header = TRUE)
    unlink(out)
    return(dsn1)
}

geteGene <- function(x="Wholeblood", chr, gene){
    require(data.table)
    eqtl.dir <- ifelse(x=="eQTLgen", "/home/weiyu/AML/AML_HRCpcspecRocket151118/eQTL",
                       "/home/weiyu/AML/AML_RocketHPC/eQTL")
    eqtl.res <-c("WholeBloodGTExV7_eQTL_V7.header.gz",
                 "LymphocytesGTExV7_eQTL_V7.header.gz",
                 "cis-eQTLs_full_20180905tabix.gz")
    ## eqtlfile <- ifelse(x=="Wholeblood", "WholeBloodGTExV7_eQTL_V7.header.gz",
    ##                    "LymphocytesGTExV7_eQTL_V7.header.gz")
    
    eqtlfile <- switch(x,
                       "Wholeblood"="WholeBloodGTExV7_eQTL_V7.header.gz",
                       "Lymphocyte"="LymphocytesGTExV7_eQTL_V7.header.gz",
                       "eQTLgen"="cis-eQTLs_full_20180905tabix.gz"
                       )


    out <- tempfile()
    tabixcmd <- paste0(" tabix -hf ", file.path(eqtl.dir, eqtlfile), " ", chr,
                       " | parallel --no-notice --pipe 'grep -w ", paste0("\"", gene, "\""), "' > ", out)
    print(tabixcmd)
    try(system(tabixcmd))
    ##if (file.info(out)$size==0)stop("no data found")
    if (file.info(out)$size==0){
        return(NULL)
        unlink(out)
    } else {
        dsn1 <- fread(out, header = FALSE)
        unlink(out)
        return(dsn1)
    }
}


geneIn <- geteQTL(interval=check_interval)
geneIds <- unique(geneIn$gene_id)
geneIds <- stringr::str_extract(geneIds, "ENSG[[:digit:]]+")

outgetx <- paste0("eQTLGen_", snpid, "_500kb.rds")

if (!file.exists(outgetx)){
    geneIdsRes <- lapply(geneIds, function(y) geteGene(x="eQTLgen", chr= snp_chr, gene=y))
    saveRDS(geneIdsRes, file=outgetx)
} else {
    geneIdsRes <- readRDS(outgetx)
}

ddexist <- sapply(geneIdsRes, function(x) !is.null(x))

genetopeQTL <- rbindlist(lapply(geneIdsRes[ddexist], function(x) x[which.min(V1)]),
                         use.names = TRUE, fill = TRUE)

genetopeQTL[V4 - snp_pos >0 & abs(V4-snp_pos) > 250000, ':=' (start=max(1, snp_pos - 5000), end=V4 +250000)]
genetopeQTL[V4 - snp_pos >0 & abs(V4-snp_pos) <=250000, ':=' (start=max(1, V4-250000), end=V4 +250000)]

genetopeQTL[V4 - snp_pos <0 & abs(V4-snp_pos) > 250000, ':=' (start=max(1, V4-250000), end=snp_pos + 5000)]
genetopeQTL[V4 - snp_pos <0 & abs(V4-snp_pos) <=250000, ':=' (start=max(1, V4-250000), end= V4 +250000)]



resdir <- paste("eQTLgen_R4N4HRC", pheno, snpid, sep="_")
if (!dir.exists(resdir)){dir.create(resdir)}
setwd(resdir)

## genekeep <- genetopeQTL[abs(V2- snp_pos) <= 250000 | V11=="AKR1B10", ]

## genekeep <- genetopeQTL[abs(V2- snp_pos) <= 250000, ]
genekeep <- genetopeQTL

genessss <- genekeep$V9


## geneIntervals <- paste0(genekeep$V1, ":", pmax(1, genekeep$V2 - 250000), "-", genekeep$V2 + 250000)
geneIntervals <- paste0(genekeep$V3, ":", genekeep$start, "-", genekeep$end)

resstore <- vector("list", length=dim(genekeep)[1])
names(resstore) <- genessss


## pdf(paste0(dataver, "_", pheno, "_", snpid, "_eQTLGen_R2.pdf"), width=10, height=5, onefile = TRUE)
## ##par(mfrow=c(1,1))
## par(mar=c(5, 4, 4, 1), mfrow=c(1, 2))  

for (n in seq_along(genessss)){

    ## our results
    ourRes <- tabixfun(dataver=dataver, datalevel="meta", phenoprefix=pheno, chrposstring=geneIntervals[n])[N==4, ]
    ## GTex results
    gtex <- geteQTL(x="eQTLgen", interval=geneIntervals[n])[GeneSymbol== genekeep[n, V9], ]
    ## 
    setnames(gtex, c("SNPChr", "SNPPos"), c("#chr", "pos_b37"))
    
    if (prod(dim(gtex))==0){next}
    
    merged.res <- merge(ourRes, gtex,
                        by.x=c("#CHR", "BP"), by.y=c("#chr", "pos_b37"),
                        suffixes=c("our", "getx"), all = FALSE)
    if (!(snp_pos %in% merged.res$BP)){next}

    subset.data <- tempfile()
    bcftools.cmd <- paste0("/home/weiyu/Downloads/bcftools-1.9build/bin/bcftools view --samples-file /home/weiyu/AML/AML_HRCpcspecRocket151118/shinyDevel/Eur.lst --regions ",
                           geneIntervals[n],
                           " -O z -o ",
                           subset.data,
                           " --threads 8 ",
                           " /home/weiyu/1KGdata/ALL.chr",
                           snp_chr, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
                           )
    try(system(bcftools.cmd))

    plinkld <- paste0("plink ",
                      " --vcf ", subset.data,
                      " --keep-fam /home/weiyu/AML/AML_HRCpcspecRocket151118/shinyDevel/Eur.lst",
                      " --r2 ",
                      " --ld-window-kb 1000 ",
                      " --ld-window 99999 ",
                      " --ld-window-r2 0 ",
                      " --ld-snp ", snpid,
                      " --out ", paste0(subset.data, snpid)
                      )
    try(system(plinkld))
    ## ld with eqtl top
    #esqnpid <- system(paste0("zgrep -w '", merged.res[which.min(pval_nominal), SNPrsid], "' ",
    #                         subset.data, " | cut -f3 "), intern = TRUE)
    esqnpid.cmd <- paste0("zgrep -w '", merged.res[order(`#Pvalue`), ][1:3, rsid], "' ",
                             subset.data, " | cut -f3 ")
    esqnpid <- as.character(unlist(sapply(esqnpid.cmd, system, intern = TRUE))[1])
    plinkld1 <- paste0("plink ",
                      " --vcf ", subset.data,
                      " --keep-fam /home/weiyu/AML/AML_HRCpcspecRocket151118/shinyDevel/Eur.lst",
                      " --r2 ",
                      " --ld-window-kb 1000 ",
                      " --ld-window 99999 ",
                      " --ld-window-r2 0 ",
                      " --ld-snp ", esqnpid,
                      " --out ", paste0(subset.data, esqnpid, "eqtl")
                      )
    try(system(plinkld1))
    
    unlink(subset.data)

    myld <- read.table(paste0(subset.data, snpid, ".ld"), header = TRUE)
    setDT(myld)

    merged.res[myld, RSQR:= i.R2, on=c(BP="BP_B")]
    merged.res$ldcol <- as.character(cut(merged.res$RSQR,
                                         breaks=c(0,0.2,0.4,0.6,0.8,1),
                                         labels=c('navy','lightskyblue','green','orange','red'),
                                         include.lowest=TRUE)
                                     )

    merged.res[rsid==snpid, ':=' (ldcol='purple3', shape=23)]
    
    merged.res[which.min(`#Pvalue`), ':=' (shape=24)]
    merged.res[is.na(shape), shape:=21]
    merged.res[is.na(RSQR), ldcol:='grey50']

    ourcols <-  c('grey50', 'navy','lightskyblue','green','orange','red', 'purple3')

    merged.res[, ldcol1 := factor(ldcol, levels=ourcols)]

    eqtlld <- read.table(paste0(subset.data, esqnpid, "eqtl.ld"), header = TRUE)
    setDT(eqtlld)

    merged.res[eqtlld, eRSQR:= i.R2, on=c(BP="BP_B")]
    merged.res$eldcol <- as.character(cut(merged.res$eRSQR,
                                         breaks=c(0,0.2,0.4,0.6,0.8,1),
                                         labels=c('navy','lightskyblue','green','orange','red'),
                                         include.lowest=TRUE)
                                     )

    merged.res[rsid==esqnpid, ':=' (eldcol='purple3', eshape=23)]
    merged.res[rsid==snpid, ':=' (eshape=24)]
    merged.res[is.na(eshape), eshape:=21]
    merged.res[is.na(eRSQR), eldcol:='grey50']
 
    merged.res[, eldcol1 := factor(eldcol, levels=ourcols)]

    cexsize <- c(rep(0.7, 4), rep(0.9, 2), rep(1.5))

    resstore[[n]] <-  merged.res

    setnames(merged.res, c("#Pvalue"), "pval_nominal")

    pdf(paste0("N4_", dataver, "_", pheno, "_", genessss[n],
                "_", snpid, "_eQTLgen_R4.pdf"), width=5, height=5, onefile = TRUE)
    par(mgp=c(2, 0.8, 0))
    plot(x= -log10(merged.res$pval_nominal),
         y=-log10(merged.res$P),
         xlim=c(0, max(2, -log10(merged.res$pval_nominal))+ 0.2 ),
         col=scales::alpha(merged.res$ldcol, 0.8), type="n",
         xlab="", ylab="",
         cex=0.9, cex.axis=0.8, las=1
         )
    

    merged.dsn1 <- split(merged.res,
                         factor(merged.res$ldcol, levels= ourcols)
                         )
    for(i in 1:length(merged.dsn1)){
        points(x=merged.dsn1[[i]][, -log10(pval_nominal)],
               y=merged.dsn1[[i]][, -log10(P)],
               bg=merged.dsn1[[i]][, scales::alpha(ldcol, 0.8)],
               pch=merged.dsn1[[i]][, shape],
               cex=ifelse(merged.dsn1[[i]]$shape==24, 1.5, cexsize[i])
               )
    }

    text(x=merged.dsn1[["purple3"]][, -log10(pval_nominal)],
         y=merged.dsn1[["purple3"]][,  -log10(P)],
         labels=merged.dsn1[["purple3"]]$rsid, pos=4, cex=0.8
         )
    text(x=merged.res[shape==24, -log10(pval_nominal)],
         y=merged.res[shape==24, -log10(P)],
         labels=merged.res[shape==24, rsid],
         pos=2, cex=0.8, adj=c(1, 1),
         srt=-45
         
         )


     ## plot(x= -log10(merged.res$pval_nominal),
    ##  y=-log10(merged.res$P),
    ##  col=scales::alpha(merged.res$eldcol, 0.8), type="n",
    ##  xlab="", ylab="",
    ##  cex=0.9
    ##  )

    ## merged.dsn1 <- split(merged.res,
    ##                      factor(merged.res$eldcol, levels= ourcols)
    ##                      )
    ## for(i in 1:length(merged.dsn1)){
    ##     points(x=merged.dsn1[[i]][, -log10(pval_nominal)],
    ##            y=merged.dsn1[[i]][, -log10(P)],
    ##            bg=merged.dsn1[[i]][, scales::alpha(eldcol, 0.8)],
    ##            pch=merged.dsn1[[i]][, eshape],
    ##            cex=ifelse(merged.dsn1[[i]]$eshape==24, 1.5, cexsize[i])
    ##            )
    ## }

    ## text(x=merged.dsn1[["purple3"]][, -log10(pval_nominal)],
    ##      y=merged.dsn1[["purple3"]][,  -log10(P)],
    ##      labels=merged.dsn1[["purple3"]]$rsid, pos=2, cex=0.8,
    ##      adj=c(1, 1),
    ##      srt=-45
    ##      )
    ## text(x=merged.res[eshape==24, -log10(pval_nominal)],
    ##      y=merged.res[eshape==24, -log10(P)],
    ##      labels=merged.res[eshape==24, rsid],
    ##      pos=4, cex=0.8
         
    ##      )

    ## mtext(text=paste0(genessss[n],"\nourResN=", ourRes[, .N], ";eQTLgenN=", gtex[, .N], ";N=", merged.res[, .N]
    ##                   ),
    ##       side =3, cex=1.2, outer = TRUE, line=-3
    ##       )
    ## mtext(text="eQTLgen -log10(p values)", side=1, outer = TRUE, line=-2)
    ## mtext(text=paste0(dataver, "_", pheno, " -log10(p values)"), side=2, outer = TRUE, line=-1.8)

    mtext(text=paste0(genessss[n]), side =3, cex=1.2, outer = TRUE, line=-3)
    mtext(text=expression("eQTLgen -log"[10]~"("~italic("P")~"value)" ), side=1, outer = TRUE, line=-3)
    mtext(text=expression(" -log"[10]~"("~italic("P")~"value)"), side=2, outer = TRUE, line=-2.5)

    kkpos <- "topright"
    if (-log10(merged.res[which.min(P)][, pval_nominal]) > (max(-log10(merged.res$pval_nominal)))/2){kkpos <- "topleft"}
    legend(kkpos, pch=21, legend=rev(c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")),
               pt.bg= rev(c('navy','lightskyblue','green','orange','red')),
               border=NULL, title="r2", bty="n", cex=0.6
                                        , x.intersp=0.5, y.intersp=0.8#, horiz = TRUE
           
           )

    dev.off()

    


}

paste0("N4_", dataver, "_", pheno, "_", snpid, "eQTLgen_R4.datainplots.txt")
names(resstore[[1]])

resstore.all <- rbindlist(resstore, use.names = TRUE, fill = TRUE)
names(resstore.all)
## resstore.all[, SNP1:=SNP][, SNP:=SNPrsid]
resstore.all[, ':=' (eQTLPvalue=pval_nominal,
                     SNPrsID=SNPgetx, Chr=`#CHR`, SNPPoshg19=BP,
                     GeneID=Gene, EffectAllele= AssessedAllele,
                     ORforAML=format(OR, digits=3), pvalueforassociationAML=P,
                     r2withsentinelSNP=RSQR, SentinelSNP=snpid,
                     EffectAlleleforGWAS=A1
                     )

             ]

resstore.all[OR < 1, ':=' (ORforAML= format(1/OR, digits=3), EffectAlleleforGWAS=A2)]
## resstore.all[, Effect_allele_eQTL:=sapply(strsplit(variant_id, "_"), function(x) x[4])]
## resstore.all[, Tstatistic := slope/slope_se]
## resstore.all[, Effect_allele_AML_GWAS := A1][OR<1, Effect_allele_AML_GWAS := A2]
names(resstore.all)
## setnames(resstore.all,
##          c("gene_name", "pval_nominal", "slope", "P", "RSQR"),
##          c("Gene", "GTEx_eQTL_Pvalue", "Normalised_effect_size", "PvalueforassociationwithAMLinmetaanalysis", "r2withSentinelAMLsnp")
##          )

## resstore.all$Cytoband <- "1p31.3"
## resstore.all$Cytoband <- "6p21.32"
## resstore.all$Cytoband <- "7q33"
## resstore.all$Cytoband <- "17p13.1"
write.table(resstore.all[, c("eQTLPvalue", "SNPrsID", "Chr", "SNPPoshg19", "GeneID",
                             "GeneSymbol", "GeneChr", "GenePos", "Zscore", "EffectAllele",
                             "OtherAllele", "NrCohorts", "NrSamples", "FDR", "ORforAML",
                             "pvalueforassociationAML", "r2withsentinelSNP",
                             "SentinelSNP", "EffectAlleleforGWAS"
                             ), with = FALSE],
            paste0("N4_", dataver, "_", pheno, "_", snpid, "eQTLgen_R4.datainplots.SuppTable2.txt"),
            quote = FALSE, sep="\t",
            row.names = FALSE, col.names = TRUE
            )

write.table(resstore.all[order(Gene, pval_nominal)][, head(.SD, n=10), by=c("Gene")][, c("eQTLPvalue", "SNPrsID", "Chr", "SNPPoshg19", "GeneID",
                             "GeneSymbol", "GeneChr", "GenePos", "Zscore", "EffectAllele",
                             "OtherAllele", "NrCohorts", "NrSamples", "FDR", "ORforAML",
                             "pvalueforassociationAML", "r2withsentinelSNP",
                             "SentinelSNP", "EffectAlleleforGWAS"
                             ), with = FALSE],
            paste0("N4_", dataver, "_", pheno, "_", snpid, "eQTLgen_R4.datainplots.top10bygeneSuppTable2.txt"),
            quote = FALSE, sep="\t",
            row.names = FALSE, col.names = TRUE
            )


q(save="no")




## require(ggplot2)
## ppp_base <- ggplot(data= merged.res, aes(x=-log10(pval_nominal), y=-log10(P)))
## ppp_base + geom_point(aes(colour=ldcol1, fill=ldcol1), alpha=0.8) + scale_fill_manual(values = ourcols, guide = "none")



sapply(1:length(merged.dsn1), function(x) 

       )


##points(x= -log10(merged.res$pval_nominal), x= -log10(merged.res$pval_nominal))





par(mar=c(2,2,2,1))
par(fig=c(0, 0.6, 0.5, 1))
plot(x=merged.res$BP, y=-log10(merged.res$pval_nominal),
     col=scales::alpha(merged.res$ldcol, 0.8), cex=0.6)
par(fig=c(0, 0.6, 0, 0.5))
par(new=T)
plot(x=merged.res$BP, y=-log10(merged.res$P),
     col=scales::alpha(merged.res$ldcol, 0.8), cex=0.6)
par(fig=c(0.6, 1, 0.3, 0.7))
par(new = T)
plot(x= -log10(merged.res$pval_nominal), y=-log10(merged.res$P),
     col=scales::alpha(merged.res$ldcol, 0.8), cex=0.6)
