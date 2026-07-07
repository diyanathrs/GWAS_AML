## subset samples for AMLsub analysis
library(dplyr)

dat1 <- read.delim('../AMLsur_Phenomerged_Dec2024_subgroups.txt', header = T) #%>% filter(!.id %in% c("BirminghamNCL3", "BirminghamNCL5"))
dat2 <- read.delim('AMLsur_final/AMLsur_Phenomerged_Dec2024_FINAL.txt') # %>% filter(!Origin=="Birmingham")
dat2$.id <- sub('Birmingham.NCL3','BirminghamNCL3',dat2$.id)
dat2$.id <- sub('Birmingham.NCL5','BirminghamNCL5',dat2$.id)

setdiff(unique(dat2$.id) , unique(dat1$.id))
study <- 'BirminghamNCL3'
setdiff(subset(dat1, dat1$.id==study)$IID ,subset(dat2, dat2$.id==study)$IID)
setdiff(subset(dat2, dat2$.id==study)$IID ,subset(dat1, dat1$.id==study)$IID)

#fix bir3 and 5 names in dat2
bir35 <- dat2 %>% filter(.id %in% c('BirminghamNCL3', 'BirminghamNCL5'))
dat2 <- dat2 %>% filter(!.id %in% c('BirminghamNCL3', 'BirminghamNCL5'))
bir35$IID <- paste0(bir35$IID,"_",bir35$IID)
dat3 <- rbind(dat2, bir35)

setdiff(dat3$IID, dat1$IID)
setdiff(dat1$IID, dat3$IID)

# add any sct to dat1
setdiff(unique(dat3$.id) , unique(dat1$.id))
setdiff(names(dat3), names(dat1))
dat4 <- merge(dat1, dat3[c(1,16)], by = 'IID', all.x = T)
dat4[duplicated(dat4$IID), ]

names(dat4)
st=format(Sys.time(), "%Y%b")
#write.table(dat4, paste0('AMLsur_final/','AMLsurvALL_',st,'.txt'), quote = F, row.names = F, col.names = T, sep = '\t')

aml_all25 <- read.delim('AMLsur_final/AMLsurvALL_2025Feb.txt', header = T)
setdiff(dat1$IID, aml_all25$IID)
aml_all25[duplicated(aml_all25$IID), ]
dat1[duplicated(dat1$IID), ]

############
### old ####
############

#investigating MIAML-239
dat3[duplicated(dat3$IID), ]
dat2[duplicated(dat2$IID), ]
dat1[duplicated(dat1$IID), ]

grep('GSM1502350_Malek_MIAML-239-N-BU_GenomeWideSNP_6_.Malek_GSM1502350_Malek_MIAML-239-N-BU_GenomeWideSNP_6_.Malek', dat1$IID)

setdiff(subset(dat2, dat2$.id=='US2')$IID, subset(dat1, dat1$.id=='US2')$IID)


table(dat1$.id)
table(dat2$.id)

setdiff(dat1$IID, dat2$IID)


table(dat1$t.15.17.,useNA = 'ifany')
dat1$APL <- 'noAPL'
dat1$APL[dat1$t.15.17.==2] <- 'APL'
table(dat1$APL,useNA = 'ifany')
noapl <- subset(dat1, dat1$APL=='noAPL')
noapl <- noapl[!duplicated(noapl$IID) ,]
#write.table(noapl$IID, 'AML_subset_noAPL.txt', quote = F, row.names = F, col.names = F)
old.noapl <- subset(dat1, dat1$t.15.17.==-9)
table(old.noapl$.id)
table(noapl$.id)

# ELN2 cases
eln2 <- subset(dat1, dat1$ELN.22==2)
eln2 <- eln2[!duplicated(eln2$IID) ,]
#write.table(eln2$IID, 'AML_subset_ELN2.txt', quote = F, row.names = F, col.names = F)

table(noapl$.id)
table(eln2$.id)

# noAPL in ELN2 groups
ELN2_noAPL <- subset(eln2, eln2$APL=='noAPL')
#table(eln2$t.15.17.)
table(ELN2_noAPL$.id)
table(eln2$.id)

#write.table(ELN2_noAPL$IID, 'AML_subset_ELN2noAPL.txt', quote = F, row.names = F, col.names = F)

#dat2$IID <- paste0(dat2$IID,"_",dat2$IID)

### Normal subgroup ##
table(dat1$Normal,useNA = 'ifany')
dat1$APL <- 'noAPL'
dat1$APL[dat1$t.15.17.==2] <- 'APL'
table(dat1$APL,useNA = 'ifany')
noapl <- subset(dat1, dat1$APL=='noAPL')
noapl <- noapl[!duplicated(noapl$IID) ,]
# subset normal now
table(noapl$Normal, useNA = 'ifany')
normal <- subset(noapl, noapl$Normal=='2')

write.table(normal$IID, 'AMLsubset_Normal.txt', quote = F, row.names = F, col.names = F)

table(normal$.id)


## data for anySCT (transplants)
# noSCT+noAPL
# noSCT+ELN2 (noAPL too)
# noSCT+Normal (noAPL too)

table(aml_all25$Any.SCT, useNA = 'ifany')
aml_all25$APL <- 'noAPL'
aml_all25$APL[aml_all25$t.15.17.==2] <- 'APL'
table(aml_all25$APL,useNA = 'ifany')
noapl <- subset(aml_all25, aml_all25$APL=='noAPL')
noapl <- noapl[!duplicated(noapl$IID) ,]
# subset noSCP now
table(noapl$Any.SCT, useNA = 'ifany')
noapl$transplant <- 0
noapl$transplant[noapl$Any.SCT >= 1] = 1
noSCT.noAPL <- subset(noapl, noapl$transplant == '0')
table(noSCT.noAPL$.id)
ELN2noSCT <- subset(noSCT.noAPL, noSCT.noAPL$ELN.22=='2')
table(ELN2noSCT$.id)
NormalnoSCT <- subset(noSCT.noAPL, noSCT.noAPL$Normal=='2')
table(NormalnoSCT$.id)

write.table(noSCT.noAPL$IID, 'AMLsubset_noSCT_noAPL.txt', quote = F, row.names = F, col.names = F)
write.table(ELN2noSCT$IID, 'AMLsubset_ELN2noSCT.txt', quote = F, row.names = F, col.names = F)
write.table(NormalnoSCT$IID, 'AMLsubset_NormalnoSCT.txt', quote = F, row.names = F, col.names = F)
