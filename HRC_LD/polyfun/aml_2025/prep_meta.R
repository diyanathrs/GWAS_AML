args <- commandArgs(trailingOnly = TRUE)

print(args[1])

dat <- read.table(args[1], head=T)
head(dat)
n.max <- max(dat$N)

#dat2 <- dat[dat$N==n.max, ]
dat <- dat[-6] 
# remove Odds ratio 0s
dat <- subset(dat, dat$OR > 0)

write.table(dat, 'sumstat.meta', sep='\t', row.names=F, quote=F)
