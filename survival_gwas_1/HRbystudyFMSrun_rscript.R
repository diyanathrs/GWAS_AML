args=(commandArgs(TRUE))

## if(length(args)==0){
##     print("No arguments supplied.")
##     ##supply default values
## } else {
##     for(i in 1:length(args)){
##         eval(parse(text=args[[i]]))
##     }
## }

#chrom <- args[1]
## ncpu <- as.numeric(Sys.getenv("NSLOTS"))
## ncpu
## chrom

chrom <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
chrom

ncpu <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
ncpu

prefix <- args[6]
prefix

genfileName <- paste0(prefix, "_chr", chrom, "subset.gen")
gensample <- paste0(prefix, "_chr", chrom, "subsetMOD.sample")
## genfileName <- paste0("OxfCLL4imp_chr", chrom, "subset.gen")
## gensample <- paste0("OxfCLL4imp_chr", chrom, "subsetMOD.sample")
#phenofile <- "HullB1_OS_Dx_Death_LFU.txt"
#cen.col <- "OS_Dx_Death_LFU_Status"
#time.col <- "OS_Dx_Death_LFU"
phenofile <- args[1]
phenofile
cen.col <- args[2]
cen.col

time.col <- args[3]
time.col

data.dir <- args[4]
data.dir

study <- args[5]
study




cat("Study is ", study, "\n")
cat("phenofile is ", phenofile, "\n")
cat("chrom", chrom, "data\n")
cat("gen file is", genfileName, "\n")
cat("sample file is", gensample, "\n" )
cat("censored variable used is", cen.col, "\n")
cat("time variable used is", time.col, "\n")
cat("No. cpus are ", ncpu, "\n")

##quit(save="no")

#subtype
kkk <- c("survival", "data.table", "doParallel")
lapply(kkk, require, character.only=TRUE)


index <- 0
chunkSize <- 20000

## checking no. columns in gen are the same as 3*no.samples
## samples <- fread(paste0(args[5], "/", gensample), header=TRUE)[-1, ]
## genfile <- file(paste0(args[5], "/", genfileName), open="r")
samples <- fread(paste0(data.dir, "/", gensample), header=TRUE)[-1, ]
genfile <- file(paste0(data.dir, "/", genfileName), open="r")

sapply(strsplit(readLines(genfile, n=1)," "), length)
stopifnot(sapply(strsplit(readLines(genfile, n=1)," "), length)==(samples[, .N]*3+6))
close(genfile)

## the type of columns
cols <- rep("numeric", samples[, .N]*3+6 )
cols[c(1,2,3,5,6)] <- "character"

## read the data for a chunk
con <- file(paste0(data.dir, "/", genfileName), open="r")


dataChunk <- read.table(con,  nrows=chunkSize, skip=0, header=FALSE, colClasses=cols)

#source("../HRbystudy.R")

#summary.out <- HRbystudy(dataChunk, phenofile, cen.col, time.col, ncpu, chrom, study)

## !is.null(dataChunk) || stop("no data provided")
## !is.null(phenofile) || stop("no pheno data provided")
## !is.null(cen.col) || stop("censored variable must be provided")
## !is.null(time.col) || stop("time variable must be provided")
## !is.null(study) || stop("study must be specified")

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
## dosageall.out <- rbind(dosageB[AA.probs.sum >= BB.probs.sum, ], dosageA[AA.probs.sum < BB.probs.sum,])

## alt allele as the effect allele
dosageall.out <- dosageB

## dosageall.out <- dosageB
colnames(dosageall.out) <- samples[, ID_1]
data.out <- as.data.frame(t(dosageall.out))
data.out$sampleId <- rownames(data.out)

## phenotypes
## phenos <- fread(phenofile, header=TRUE)
## stopifnot(c("id", cen.col, time.col, "Study") %chin% names(phenos))
phenos <- fread(phenofile, header=TRUE)[, c(".id", "IID", cen.col, time.col), with=F]
stopifnot(c(cen.col, time.col) %chin% names(phenos))
    ## keep columns wanted
    #phenos <- phenos[, c("id", cen.col, time.col, "Study"), with=FALSE]
    #surv.sub <- merge(data.out, phenos, by.x=c("sampleId"), by.y="id")
    #markers <- setdiff(names(surv.sub), c("sampleId", cen.col, time.col, "Study"))
    #setDT(surv.sub)
dsn1 <- merge(data.out, phenos, by.x=c("sampleId"), by.y="IID")
markers <- setdiff(names(dsn1), c("sampleId", "sex", "IID", ".id" , cen.col, time.col))

    #surv.sub <- merge(data.out, phenos, by.x=c("sampleId"), by.y="IID")
    #markers <- setdiff(names(surv.sub), c("sampleId", "sex", "IID", ".id" , cen.col, time.col))
setDT(dsn1)
dsn1.clean <- dsn1[complete.cases(dsn1[, c(".id", cen.col, time.col), with=F]), ]

#################
## study level ##
#################

surv.sub <- dsn1.clean[.id==study, ]

## check minor allele counts, assume event is recoded as 1
ma.counts0 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==1"))), markers, with=FALSE], 2, sum)
abc <- data.table(macountEvent=ma.counts0, snp=names(ma.counts0))
ma.counts1 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==0"))), markers, with=FALSE], 2, sum)
abc1 <- data.table(macountCensored=ma.counts1, snp=names(ma.counts1))
ma.summary <- merge(abc, abc1, by="snp")
snp.info[ma.summary, macountEvent:=i.macountEvent, on=c(rsid="snp")][
    ma.summary, macountCensored:=i.macountCensored, on=c(rsid="snp")]
## study.maf
study.maf <- apply(surv.sub[, markers, with=FALSE], 2, sum)/(2*surv.sub[, .N])
study.maf <- study.maf[match(snp.info$rsid, names(study.maf))]
## event maf
study.event.maf <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==1"))), markers, with=FALSE], 2, sum)/(2*surv.sub[eval(parse(text=paste0(cen.col, "==1"))), .N])
study.event.maf <- study.event.maf[match(snp.info$rsid, names(study.event.maf))]
## censor maf
study.censor.maf <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==0"))), markers, with=FALSE], 2, sum)/(2*surv.sub[eval(parse(text=paste0(cen.col, "==0"))), .N])
study.censor.maf <- study.censor.maf[match(snp.info$rsid, names(study.censor.maf))]
    

snp.info[, ':=' (studyMAF = study.maf, eventMAF = study.event.maf, censoredMAF = study.censor.maf)]

    #stopifnot(all(names(study.event.maf)==names(study.maf)), all(names(study.censor.maf)==names(study.maf)))
    #snps.count2 <- names(ma.counts0[ma.counts0 >=2])
    #snps.count2 <- snp.info[eventMAF >=0.01 & censoredMAF>=0.01 , rsid]
#snps.count2 <- snp.info[macountEvent >=1 & macountCensored >=1, rsid]
snps.count2 <- snp.info[macountEvent >=1, rsid]
if (length(snps.count2)>0) {
    registerDoParallel(cores=ncpu)
    res.out <- foreach(n=seq_along(snps.count2), .combine=rbind, .export="snps.count2") %dopar% { 
        tryCatch(
        {snps.wanted <- snps.count2[n]
            model <- paste0("Surv(", time.col,", ", cen.col, "==1) ~  ", snps.wanted)
            cox.res <- coxph(as.formula(model), data=surv.sub, method="breslow")
            cru.p <-  format(summary(cox.res)$coefficients[snps.wanted ,c("Pr(>|z|)")], digit=3, nsmall=3)
            cruci <- format(summary(cox.res)$conf.int[snps.wanted, c("exp(coef)","lower .95","upper .95")], digit=3)
            cru.n <- summary(cox.res)$n
            cru.nevent <- summary(cox.res)$nevent
            crude.hr <- data.table(SNP=snps.wanted,  totN=cru.n, No.event=cru.nevent,
                                   crudeHR=cruci[1], crudeHR_L95=cruci[2],
                                   crudeHR_H95=cruci[3], crude.pvalue=cru.p,
                                   beta=summary(cox.res)$coefficients[snps.wanted , c("coef")],
                                   SE=summary(cox.res)$coefficients[snps.wanted , c("se(coef)")])
        },
        error= function(e){
            cat(paste0(snps.wanted, ", Error ", conditionMessage(e), sep="\n"))},
        warning = function(cond){
            NULL
            ##cat(paste0(snps.wanted, ", Warning ", conditionMessage(cond), sep="\n"))
        }
        )
    }
    stopImplicitCluster()
    gc()
    
}

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

write.table(summary.out, paste0(study, "_", cen.col, "_chr", chrom, "_HRCSurRes.out"), append=FALSE, sep="\t",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

repeat {
        index <- index + 1
        print(paste('Processing rows:', index * chunkSize))
 
        if (nrow(dataChunk) != chunkSize){
                print('Processed all rows!')
                break}

        dataChunk <- read.table(con,  nrows=chunkSize, skip=0, header=FALSE, colClasses=cols)
        #summary.out <- HRbystudy(dataChunk, phenofile, cen.col, time.col, ncpu, chrom, study)

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
        ## dosageall.out <- rbind(dosageB[AA.probs.sum >= BB.probs.sum, ], dosageA[AA.probs.sum < BB.probs.sum,])
        dosageall.out <- dosageB
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
        dsn1 <- merge(data.out, phenos, by.x=c("sampleId"), by.y="IID")
        markers <- setdiff(names(dsn1), c("sampleId", "sex", "IID", ".id" , cen.col, time.col))

    #surv.sub <- merge(data.out, phenos, by.x=c("sampleId"), by.y="IID")
    #markers <- setdiff(names(surv.sub), c("sampleId", "sex", "IID", ".id" , cen.col, time.col))
        setDT(dsn1)
        dsn1.clean <- dsn1[complete.cases(dsn1[, c(".id", cen.col, time.col), with=F]), ]

#################
## study level ##
#################

        surv.sub <- dsn1.clean[.id==study, ]

        ## check minor allele counts, assume event is recoded as 1
        ma.counts0 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==1"))), markers, with=FALSE], 2, sum)
        abc <- data.table(macountEvent=ma.counts0, snp=names(ma.counts0))
        ma.counts1 <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==0"))), markers, with=FALSE], 2, sum)
        abc1 <- data.table(macountCensored=ma.counts1, snp=names(ma.counts1))
        ma.summary <- merge(abc, abc1, by="snp")
        snp.info[ma.summary, macountEvent:=i.macountEvent, on=c(rsid="snp")][
            ma.summary, macountCensored:=i.macountCensored, on=c(rsid="snp")]
        ## study.maf
        study.maf <- apply(surv.sub[, markers, with=FALSE], 2, sum)/(2*surv.sub[, .N])
        study.maf <- study.maf[match(snp.info$rsid, names(study.maf))]
        ## event maf
        study.event.maf <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==1"))), markers, with=FALSE], 2, sum)/(2*surv.sub[eval(parse(text=paste0(cen.col, "==1"))), .N])
        study.event.maf <- study.event.maf[match(snp.info$rsid, names(study.event.maf))]
        ## censor maf
        study.censor.maf <- apply(surv.sub[eval(parse(text=paste0(cen.col, "==0"))), markers, with=FALSE], 2, sum)/(2*surv.sub[eval(parse(text=paste0(cen.col, "==0"))), .N])
        study.censor.maf <- study.censor.maf[match(snp.info$rsid, names(study.censor.maf))]
    

        snp.info[, ':=' (studyMAF = study.maf, eventMAF = study.event.maf, censoredMAF = study.censor.maf)]

    #stopifnot(all(names(study.event.maf)==names(study.maf)), all(names(study.censor.maf)==names(study.maf)))
    #snps.count2 <- names(ma.counts0[ma.counts0 >=2])
    #snps.count2 <- snp.info[eventMAF >=0.01 & censoredMAF>=0.01 , rsid]
    # snps.count2 <- snp.info[macountEvent >=1 & macountCensored >=1, rsid]
        snps.count2 <- snp.info[macountEvent >=1, rsid]
        if (length(snps.count2)>0) {
            registerDoParallel(cores=ncpu)
            res.out <- foreach(n=seq_along(snps.count2), .combine=rbind, .export="snps.count2") %dopar% { 
                tryCatch(
                {snps.wanted <- snps.count2[n]
                    model <- paste0("Surv(", time.col,", ", cen.col, "==1) ~  ", snps.wanted)
                    cox.res <- coxph(as.formula(model), data=surv.sub, method="breslow")
                    cru.p <-  format(summary(cox.res)$coefficients[snps.wanted ,c("Pr(>|z|)")], digit=3, nsmall=3)
                    cruci <- format(summary(cox.res)$conf.int[snps.wanted, c("exp(coef)","lower .95","upper .95")], digit=3)
                    cru.n <- summary(cox.res)$n
                    cru.nevent <- summary(cox.res)$nevent
                    crude.hr <- data.table(SNP=snps.wanted,  totN=cru.n, No.event=cru.nevent,
                                           crudeHR=cruci[1], crudeHR_L95=cruci[2],
                                           crudeHR_H95=cruci[3], crude.pvalue=cru.p,
                                           beta=summary(cox.res)$coefficients[snps.wanted , c("coef")],
                                           SE=summary(cox.res)$coefficients[snps.wanted , c("se(coef)")])
                },
                error= function(e){
                    cat(paste0(snps.wanted, ", Error ", conditionMessage(e), sep="\n"))},
                warning= function(cond){
                    NULL
                    ## cat(paste0(snps.wanted, ", Warning ", conditionMessage(cond), sep="\n"))
                }
                )
            }
            stopImplicitCluster()
            gc()
        }

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

        write.table(summary.out, file=paste0(study, "_", cen.col, "_chr", chrom, "_HRCSurRes.out"),
                    append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
close(con)
       

