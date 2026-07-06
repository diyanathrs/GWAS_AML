# code to get case and control iids
# cases are coded as 2, controls are 1. 0 is trimmed

require(data.table)
dat <- read.table('AML_casecon.txt',header=T)

head(dat)

dat[]
