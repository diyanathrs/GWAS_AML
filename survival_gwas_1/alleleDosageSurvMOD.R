
gen2dosageSurv <- function(dataChunk, phenofile, cen.col, time.col, ncpu, chrom,...){
    !is.null(dataChunk) || stop("no data provided")
    !is.null(phenofile) || stop("no pheno data provided")
    !is.null(cen.col) || stop("censored variable must be provided")
    !is.null(time.col) || stop("time variable must be provided")
    !is.null(ncpu) || stop("no. cpus musted be specified")
    #argnames <- list(...)
    #if (!("ncpu" %in% argnames)){ncpu=detectCores()}
    #print(c('ncpu =', ncpu))

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
    ## odds.ids <- dataChunk$V2[grepl("-", dataChunk$V2)]
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

    phenos <- fread(phenofile, header=TRUE)[, c(".id", "IID", cen.col, time.col), with=F]
    stopifnot(c(cen.col, time.col) %chin% names(phenos))

    phenos1 <- phenos[match(samples[, ID_1], phenos$IID), ]
    stopifnot(all(phenos1$IID == samples[, ID_1]))
    Sampevent1 <- phenos1$IID[which(eval(parse(text=paste0("phenos1$", cen.col, "==1"))))]
    Sampevent0 <- phenos1$IID[which(eval(parse(text=paste0("phenos1$", cen.col, "==0"))))]
    macount1 <- apply(dosageall.out[, Sampevent1], 1, sum)
    macount0 <- apply(dosageall.out[, Sampevent0], 1, sum)
    macount1.order <- macount1[match(snp.info$rsid, names(macount1))]
    macount0.order <- macount0[match(snp.info$rsid, names(macount0))]
    stopifnot(all(names(macount1.order)==snp.info$rsid), all(names(macount0.order)==snp.info$rsid))
    snp.info[, ':=' (macountEvent=macount1.order,  macountCensored=macount0.order)]


    snps.count2 <- names(macount1.order)[which(macount1.order >=2)]

    if (length(snps.count2)>0) {
        registerDoParallel(cores=ncpu)
        res.out <- foreach(n=seq_along(snps.count2), .combine=rbind, .export="snps.count2") %dopar% {
            tryCatch(
            {
                surv.sub.1 <- data.table(.id=phenos1[[".id"]], status=phenos1[[cen.col]], time=phenos1[[time.col]], snp=dosageall.out[snps.count2[n], ])
                setnames(surv.sub.1, c("status", "time", "snp"), c(cen.col, time.col, snps.count2[n]))
                ## check if allele counts==0 in any one of study
                ## if (all(table(dosageall.out[snps.count2[n], ] !=0, phenos1$.id)["TRUE", ] !=0) ){
                model <- paste0("Surv(", time.col,", ", cen.col, "==1) ~ strata(.id) + ", snps.count2[n])
                cox.res <- coxph(as.formula(model), data=surv.sub.1, method="breslow")
                ahr.p <-  format(summary(cox.res)$coefficients[snps.count2[n],c("Pr(>|z|)")], digit=3, nsmall=3)
                ahrci <- format(summary(cox.res)$conf.int[snps.count2[n], c("exp(coef)","lower .95","upper .95")], digit=3)
                ahr.n <- summary(cox.res)$n
                ahr.nevent <- summary(cox.res)$nevent
            ## } else {
            ##    ahr.p <-  NA
            ##    ahrci <- NA
            ##    ahr.n <- NA
            ##    ahr.nevent <- NA   
            ## }

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
            },
            error= function(e){
                cat(paste0(snps.count2[n], ", Error ", conditionMessage(e), sep="\n"))
            },
            warning = function(cond){
                ## cat(paste0(snps.wanted, ", Warning ", conditionMessage(cond), sep="\n"))}
                NULL
            }
            )                               #trycatch
            
        }                                   # res.out
        stopImplicitCluster()
    }                                      #snps.count2 loop


    if (is.null(res.out) | length(snps.count2)==0){
        
        res.out <- data.table(SNP="faked",  totN=NA, No.event=NA,
                          crudeHR=NA, crudeHR_L95=NA,
                          crudeHR_H95=NA, crude.pvalue=NA,
                          beta=NA,
                          SE=NA)
    }
    setDT(res.out)
    summary.out <- merge(snp.info, res.out, by.x="rsid", by.y="SNP", all.x=TRUE)
    setorder(summary.out, posb37)
    return(summary.out)


}
