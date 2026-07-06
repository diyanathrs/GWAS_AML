
gen2dosageSurv <- function(dataChunk, phenofile, cen.col, time.col, ncpu, chrom,...){
    !is.null(dataChunk) || stop("no data provided")
    !is.null(phenofile) || stop("no pheno data provided")
    !is.null(cen.col) || stop("censored variable must be provided")
    !is.null(time.col) || stop("time variable must be provided")
    !is.null(ncpu) || stop("no. cpus musted be specified")
    #argnames <- list(...)
    #if (!("ncpu" %in% argnames)){ncpu=detectCores()}
    #print(c('ncpu =', ncpu))

    ## for rsid==. then assign chr_pos
    dataChunk$V2 <- ifelse(dataChunk$V2==".", paste0(chrom, ":", dataChunk$V4), dataChunk$V2)
    dataChunk$V2 <- paste0("chr", dataChunk$V2)
    if (any(duplicated(dataChunk$V2))){
        ## identify duplicated records, same snp id but different pos
        dup.ids <- dataChunk$V2[duplicated(dataChunk$V2) | duplicated(dataChunk$V2, fromLast=TRUE)]
        dataChunk$V2 <- ifelse(dataChunk$V2 %in% dup.ids,
                               paste0(dataChunk$V2, "_", dataChunk$V5, "_", dataChunk$V6),
                               dataChunk$V2)
    }
    ## removal of "-" in the marker names
    if(any(grepl("-", dataChunk$V2))){
        odds.ids <- dataChunk$V2[grepl("-", dataChunk$V2)]
        dataChunk$V2 <- ifelse(dataChunk$V2 %in% odds.ids, gsub("-", ".", dataChunk$V2), dataChunk$V2)

    }
    ## replace ":" with "_"
     if(any(grepl(":", dataChunk$V2))){
        #odds.ids <- dataChunk$V2[grepl("-", dataChunk$V2)]
        dataChunk$V2 <- gsub(":", "_", dataChunk$V2)

    }
    
    ## dosage
    start.AA <- 1 + 3*(1-1) + 6
    end.AA <- 1 + 3*(samples[, .N] -1) + 6
    AA.columns <- seq(start.AA, end.AA, by=3)

    start.AB <- 2 + 3*(1-1) + 6
    end.AB <- 2 + 3*(samples[, .N] -1) + 6
    AB.columns <- seq(start.AB, end.AB, by=3)

    start.BB <- 3 + 3*(1-1) + 6
    end.BB <- 3 + 3*(samples[, .N] -1) + 6
    BB.columns <- seq(start.BB, end.BB, by=3)

    AA.probs <- as.matrix(dataChunk[, AA.columns])
    AB.probs <- as.matrix(dataChunk[, AB.columns])
    BB.probs <- as.matrix(dataChunk[, BB.columns])

    class(AA.probs) <- "numeric"
    class(AB.probs) <- "numeric"
    class(BB.probs) <- "numeric"

    ## check which allele is minor
    AA.probs.sum <- apply(AA.probs, 1, sum)
    BB.probs.sum <- apply(BB.probs, 1, sum)
    snp.info <- setDT(dataChunk[, 1:6])
    snp.info[AA.probs.sum >= BB.probs.sum, minorAllele:=V6][AA.probs.sum < BB.probs.sum, minorAllele:=V5]

    ## check if snp names is not .
                                        #stopifnot(!any(grepl("\\.", dataChunk$V2)))

    dosageB <- AB.probs + 2*BB.probs
    rownames(dosageB) <- dataChunk[, 2]

    dosageA <- AB.probs + 2*AA.probs
    rownames(dosageA) <- dataChunk[, 2]

    dosageAB <- (dosageA+dosageB)
    dosageAB.sum <- apply(dosageAB, 1, sum)
    mafA <- apply(dosageA, 1, sum)/dosageAB.sum
    mafB <- apply(dosageB, 1, sum)/dosageAB.sum

    snp.info[AA.probs.sum >= BB.probs.sum, MAF :=mafB[AA.probs.sum >= BB.probs.sum]][AA.probs.sum < BB.probs.sum, MAF:=mafA[AA.probs.sum < BB.probs.sum]]
    setnames(snp.info, paste0("V", seq(1,6)), c("chr", "rsid", "snp", "posb37", "alleleA", "alleleB"))

    ## if sum(AA.probs.sum) >= sum(BB.probs.sum) then A is major allele and B is minor allele
    ## if above is TRUE, B dosage is used, above is FALSE, A dosage is used
    ## dosages based on the minor allele
    dosageall.out <- rbind(dosageB[AA.probs.sum >= BB.probs.sum, ], dosageA[AA.probs.sum < BB.probs.sum,])
    ## dosageall.out <- dosageB
    colnames(dosageall.out) <- samples[, ID_1]
    data.out <- as.data.frame(t(dosageall.out))
    data.out$sampleId <- rownames(data.out)

    ## phenotypes
    #phenos <- fread(phenofile, header=TRUE)
    #stopifnot(c("id", cen.col, time.col, "Study") %chin% names(phenos))
    phenos <- fread(phenofile, header=TRUE)[, c(".id", "IID", cen.col, time.col), with=F]
    stopifnot(c(cen.col, time.col) %chin% names(phenos))
    ## keep columns wanted
    #phenos <- phenos[, c("id", cen.col, time.col, "Study"), with=FALSE]
    #surv.sub <- merge(data.out, phenos, by.x=c("sampleId"), by.y="id")
    #markers <- setdiff(names(surv.sub), c("sampleId", cen.col, time.col, "Study"))
    #setDT(surv.sub)
    surv.sub <- merge(data.out, phenos, by.x=c("sampleId"), by.y="IID")
    markers <- setdiff(names(surv.sub), c("sampleId", "sex", "IID", ".id" , cen.col, time.col))
    setDT(surv.sub)

    ## check minor allele counts, assume event is recoded as 1
    ma.counts0 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==1"))), markers, with=FALSE], 2, sum)
    abc <- data.table(macountEvent=ma.counts0, snp=names(ma.counts0))
    ma.counts1 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==0"))), markers, with=FALSE], 2, sum)
    abc1 <- data.table(macountCensored=ma.counts1, snp=names(ma.counts1))
    ma.summary <- merge(abc, abc1, by="snp")
    snp.info[ma.summary, macountEvent:=i.macountEvent, on=c(rsid="snp")][
        ma.summary, macountCensored:=i.macountCensored, on=c(rsid="snp")]

    snps.count2 <- names(ma.counts0[ma.counts0 >=2])
    if (length(snps.count2)>0) {
        registerDoParallel(cores=ncpu)
        res.out <- foreach(n=seq_along(snps.count2), .combine=rbind, .export="snps.count2") %dopar% { 
            surv.sub.1 <- surv.sub[, c("sampleId", ".id", cen.col, time.col, snps.count2[n]), with=FALSE]
            model <- paste0("Surv(", time.col,", ", cen.col, "==1) ~ strata(.id) + ", snps.count2[n])
            cox.res <- coxph(as.formula(model), data=surv.sub.1, method="breslow")
            ahr.p <-  format(summary(cox.res)$coefficients[snps.count2[n],c("Pr(>|z|)")], digit=3, nsmall=3)
            ahrci <- format(summary(cox.res)$conf.int[snps.count2[n], c("exp(coef)","lower .95","upper .95")], digit=3)
            ahr.n <- summary(cox.res)$n
            ahr.nevent <- summary(cox.res)$nevent

            model <- paste0("Surv(", time.col,", ", cen.col, "==1) ~  ", snps.count2[n])
            cox.res <- coxph(as.formula(model), data=surv.sub.1, method="breslow")
            cru.p <-  format(summary(cox.res)$coefficients[snps.count2[n],c("Pr(>|z|)")], digit=3, nsmall=3)
            cruci <- format(summary(cox.res)$conf.int[snps.count2[n], c("exp(coef)","lower .95","upper .95")], digit=3)
            cru.n <- summary(cox.res)$n
            cru.nevent <- summary(cox.res)$nevent
        
            adjusted.hr <- data.table(SNP=snps.count2[n], totN=ahr.n, No.event=ahr.nevent, totNc=cru.n,
                                      No.eventc=cru.nevent, crudeHR=cruci[1], crudeHR_L95=cruci[2],
                                      crudeHR_H95=cruci[3], crude.pvalue=cru.p,
                                  aHR=ahrci[1], aHR_L95=ahrci[2], aHR_H95=ahrci[3], a.pvalue=ahr.p)
            
        }
        stopImplicitCluster()
        setDT(res.out)
        summary.out <- merge(snp.info, res.out, by.x="rsid", by.y="SNP", all.x=TRUE)
        setorder(summary.out, posb37)
    }# all hrs 
}
