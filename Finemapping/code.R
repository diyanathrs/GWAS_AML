library(data.table)
library(susieR)
set.seed(1)


t15.sum.st <- read.table('t15_17_hit1_chr13.z', sep = ' ', header = T)
head(t15.sum.st)

maf  <- t15.sum.st$maf
bhat <- t15.sum.st$beta
shat <- t15.sum.st$se
z    <- bhat/shat

t15.ld <- as.matrix(fread('HRCout_t15_17_hit1_chr13.ld'))

fit1 <- susie_rss(z,t15.ld,n = 800,min_abs_corr = 0.1,refine = FALSE,
                  verbose = TRUE)

## Finemap output
par(mar = c(4,4,1,1))
finemap <- read.table("status_hit1_chr2.cred7",header = TRUE)
pip   <- rep(0,1001)
cs1   <- finemap$cred1
rows1 <- which(!is.na(cs1))
cs1   <- cs1[rows1]
pip[cs1] <- pip[cs1] + finemap$prob1[rows1]
plot(1:1001,pip,pch = 20,cex = 0.8,ylim = c(0,0.1),
     xlab = "SNP",ylab = "finemap PIP")
points(cs1,pip[cs1],pch = 1,cex = 1,col = "cyan")
points(vars,pip[vars],pch = 2,cex = 0.8,col = "tomato")