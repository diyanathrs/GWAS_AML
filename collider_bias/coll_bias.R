library(devtools)
#install_github("qingyuanzhao/mr.raps")
#install_github("SiyangCai/ColliderBias")

library(ColliderBias)
library(dplyr)

# Load the test dataset
data(testData)

# Adjust for collider bias using instrumental effect regression,
# and weak instrument bias using CWLS.
test <- methodCB(testData$dbeta, testData$dse, testData$ybeta, testData$yse, method = "cwls")
test
test$b.raw
#naive
naive <- lm(testData$ybeta ~ testData$xbeta, weights = 1/testData$yse^2, data = testData)
summary(naive)
coef(naive)["testData$xbeta"]
coef(naive)["old$B.risk"]

# import beta and se of prog.gwas - snp of int - dnmt3A  rs4665765 (2:25362515_T_C)
# use noAPL prog
# gwas sumstats or meta sumstats to use? only meta sumstats can be used
# dudbridge used all autosomal snps with R2 => 0.98 imp and --indep-pairwise 250 25 0.1
# before Ld pruning in rocket, make sure to filter combined HRC panel for meta only snps
# Get chr-wise snp ids
ld.prune.lst <- list.files('ld_pruned/', pattern = paste0('*.in'), full.names = T)
ld.pruned <- do.call("rbind", lapply(ld.prune.lst, function(fn) read.table(fn, header = F)))
head(ld.pruned)
grep('2:25362515_T_C',ld.pruned$V1, value = T)

prog.gwas <- read.table('../assoc_res/OSstatus_NCL_AMLsur_noAPL.meta', header = T) 
#prog.gwas.filp <- read.table('OSstatus_NCL_AMLsur_noAPL.meta', header = T) 

risk.gwas <- read.table('../../AMLmeta_results/status_NCL_PCspeAMLHRC.meta.gz', header = T) 
grep('2:25362515_T_C',risk.gwas$SNP, value = T)
setdiff(names(prog.gwas), names(risk.gwas))

# check common snps between 2 gwas
length(intersect(risk.gwas$SNP, prog.gwas$SNP))

AML.snps <- merge(risk.gwas, prog.gwas, by='SNP', all.x=T) # use this if not want to isolate snp

#calculate Beta and SE
AML.snps$B.risk <- log(AML.snps$OR.x)
AML.snps$B.prog <- log(AML.snps$OR.y)
z <- qnorm(1 - AML.snps$P.x/2)
# Compute SE
AML.snps$SE.risk <- abs(AML.snps$B.risk) / z
z <- qnorm(1 - AML.snps$P.y/2)
# Compute SE
AML.snps$SE.prog <- abs(AML.snps$B.prog) / z
head(AML.snps)

flip_beta = F
# do this if only want to flip beta based on maf
## add maf data to switch beta of maf > 0.5
if (isTRUE(flip_beta)) {
gwas.maf <- readRDS('GWAS_info_merge.RDS')
head(gwas.maf)
head(AML.snps)

AML.snps <- merge(AML.snps, gwas.maf[c(3,13,14,17)], by='SNP')
#check maf > 0.5 snps
test <- which(AML.snps$all_maf > 0.5)
#test[grep('2:25362515_T_C', test$SNP),]
#change beta.prog if maf higher
AML.snps$B.prog[test] <- -(AML.snps$B.prog[test])
}

#only flip our snp of int
#AML.snps$B.risk[grep('2:25362515_T_C', AML.snps$SNP)] <- -(AML.snps$B.risk[grep('2:25362515_T_C', AML.snps$SNP)])

# remove NA
AML.snps <- AML.snps[!is.na(AML.snps$SE.prog),]
AML.snps <- AML.snps[!is.na(AML.snps$SE.risk),]
#remove 0
AML.snps <- AML.snps %>% filter(SE.prog > 0)
AML.snps <- AML.snps %>% filter(SE.risk > 0)
range(AML.snps$SE.prog)
range(AML.snps$SE.risk)

#load ld pruned snps list
prune.idx <- which(AML.snps$SNP %in% ld.pruned$V1)
head(prune.idx)
#check
table(ld.pruned$V1 %in% risk.gwas$SNP)
table(ld.pruned$V1 %in% prog.gwas$SNP)

# perform cwls
cwls_out <- methodCB(xbeta = AML.snps$B.risk, xse = AML.snps$SE.risk, ybeta = AML.snps$B.prog, yse = AML.snps$SE.prog, method = "CWLS",
                    prune = prune.idx)
cwls_out$b
cwls_out$b.raw

length(cwls_out$ybeta.adj)
#trace(methodCB, edit = T)

AML.snps$beta.adj <- cwls_out$ybeta.adj
AML.snps$se.adj <- cwls_out$yse.adj
AML.snps$p.adj <- cwls_out$yp.adj
AML.snps$adj.HR <- exp(AML.snps$beta.adj)
AML.snps$ci_low <- exp(AML.snps$beta.adj - 1.96 * AML.snps$se.adj)
AML.snps$ci_high <- exp(AML.snps$beta.adj + 1.96 * AML.snps$se.adj)

# check OR and P trends
range(AML.snps$P.x)
names(AML.snps)
aml.plt <- AML.snps %>% filter(if_any(c(P.x, P.y, p.adj), ~ . < 1e-2))
range(aml.plt$p.adj)

AML.snps[c(1,30:33)][grep('2:25362515_T_C', AML.snps$SNP),]
AML.snps[grep('2:25362515_T_C', AML.snps$SNP),]

#naive
naive <- lm(AML.snps$B.prog ~ AML.snps$B.risk, weights = 1/AML.snps$SE.prog^2, data = AML.snps)
summary(naive)
coef(naive)["AML.snps$B.risk"]

library(ggplot2)
# p <- ggplot(AML.snps, aes(x = B.risk, y = B.prog)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", aes(weight = 1/(SE.risk^2)), se = FALSE) +
#   labs(title = "Naive regression of ybeta on dbeta (no correction)",
#        x = "Observed dbeta (noisy)",
#        y = "Observed ybeta") +
#   theme_minimal()
# 
# ggsave('test.png')

########################
# check old col bias rds
##########################
old <- readRDS('AML_collider_corrected.Rds')
# get 52K IV result again
cwls_old <- methodCB(xbeta = old$B.risk, xse = old$SE.risk, ybeta = old$B.prog, yse = old$SE.prog, method = "CWLS")
cwls_old$b
cwls_old$b.raw

cwls_old[grep('2:25362515_T_C', cwls_old$SNP),]

#naive
naive <- lm(old$B.prog ~ old$B.risk, weights = 1/old$SE.prog^2, data = old)
summary(naive)
coef(naive)["old$B.risk"]
naive.adj <- lm(old$beta.adj ~ old$B.risk, weights = 1/old$se.adj^2, data = old)
coef(naive.adj)["old$B.risk"]

p.adjust(3.89e-06, n=6, method = 'BH')
