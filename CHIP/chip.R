library(dplyr)
library(parallel)
library(readxl)


###########################
## using filtered meta ####
###########################

# dataset1
chip.dat1 <- read_xlsx('CHiP GWAS results.xlsx',sheet = 2)

# combine meta subtypes
meta.lst <- list.files(pattern = '*meta.gz')
meta.all <- mclapply(meta.lst, read.table, header=T, mc.cores = 12)

meta.lst.names <- sub("_NCL_PCspeAMLHRC.meta.gz","",meta.lst)

for (i in seq_len(length(meta.all))) {
  meta.all[[i]]$subtype <- paste0(meta.lst.names[i])
  print(head(meta.all[[i]]))
}

meta.all <- do.call(rbind,meta.all)
#save meta res
#saveRDS(meta.all,'metaRes_all.rds')

# check for chip snps
head(chip.dat1)
chip.dat1$rsid <- paste0(chip.dat1$Chr,':',chip.dat1$Pos)
meta.all$rsid <- paste0(meta.all$CHR,':',meta.all$BP)
head(meta.all)
head(chip.dat1)

# get matching snps for chip.dat1
test <- merge(chip.dat1[46], meta.all, by='rsid')

# dataset2
chip.dat2 <- read_xlsx('CHiP GWAS results.xlsx',sheet = 2)
chip.dat2$chr <- sub('chr','',chip.dat2$CHR)
chip.dat2$rsid <- paste0(chip.dat2$chr,':',chip.dat2$POS)
head(chip.dat2)
merged.dat2 <- merge(chip.dat2, meta.all, by='rsid')


############################
## use unfilted meta hits ##
############################
# dataset1
chip.dat1 <- read.delim('chip_input_April25.txt', header = T)
head(chip.dat1)

# combine meta subtypes
meta.lst <- list.files(path = '../AMLmeta_results/unfiltered',pattern = '*.meta', full.names = T)
meta.lst <- meta.lst[c(7,8)]
meta.all <- mclapply(meta.lst, read.table, header=T, mc.cores = 2)

meta.lst.names <- sub("_NCL_PCspeAMLHRC.meta","",meta.lst)
meta.lst.names <- sub(".*unfiltered/","",meta.lst.names)

for (i in seq_len(length(meta.all))) {
  meta.all[[i]]$subtype <- paste0(meta.lst.names[i])
  print(head(meta.all[[i]]))
}

meta.all <- do.call(rbind,meta.all)

head(chip.dat1)
chip.dat1$rsid <- paste0(chip.dat1$CHR,':',chip.dat1$POS)
meta.all$rsid <- paste0(meta.all$CHR,':',meta.all$BP)
head(meta.all)
head(chip.dat1)

# get matching snps for chip.dat1
merged.dat1 <- merge(chip.dat1, meta.all, by='rsid')


# add SE and CI
calculate_CI_from_p_OR <- function(OR, p_value, conf.level = 0.95) {
  # Convert OR to log scale
  log_OR <- log(OR)
  
  # Compute Z-score from p-value (two-tailed test)
  Z <- abs(qnorm(p_value / 2, lower.tail = FALSE))
  
  # Calculate standard error
  SE_log_OR <- abs(log_OR / Z)
  
  # Compute confidence interval on log scale
  log_CI_lower <- log_OR - qnorm(1 - (1 - conf.level) / 2) * SE_log_OR
  log_CI_upper <- log_OR + qnorm(1 - (1 - conf.level) / 2) * SE_log_OR
  
  # Convert back to OR scale
  CI_lower <- exp(log_CI_lower)
  CI_upper <- exp(log_CI_upper)
  
  # Return a named vector
  return(c(CI_lower, CI_upper))
}

merged.dat1 <- merged.dat1 %>%
  rowwise() %>%
  mutate(
    CI_lower = calculate_CI_from_p_OR(OR.y, P.y)[1],
    CI_upper = calculate_CI_from_p_OR(OR.y, P.y)[2]
  ) %>%
  ungroup()


#write.table(merged.dat1, 'CHIP_AML_Apr25.txt', quote = F, row.names = F, col.names = T, sep = '\t')

##############
## find snps##
##############
status.meta <- read.table('unfiltered/status_NCL_PCspeAMLHRC.meta', header = T)

# combine meta subtypes
meta.lst <- list.files(path = 'unfiltered',pattern = '*.meta', full.names = T)
meta.all <- mclapply(meta.lst, read.table, header=T, mc.cores = 11)

meta.lst.names <- sub("_NCL_PCspeAMLHRC.meta","",meta.lst)
meta.lst.names <- sub("unfiltered/","",meta.lst.names)
for (i in seq_len(length(meta.all))) {
  meta.all[[i]]$subtype <- paste0(meta.lst.names[i])
  print(head(meta.all[[i]]))
}

meta.all <- do.call(rbind,meta.all)

rsid <- 'rs117525685'
head(status.meta)
idx <- grep('10:50675568_', meta.all$SNP)
out <- meta.all[idx,]
out <- out %>%
  rowwise() %>%
  mutate(
    CI_lower = calculate_CI_from_p_OR(OR, P)[1],
    CI_upper = calculate_CI_from_p_OR(OR, P)[2]
  ) %>%
  ungroup()
out$rsid <- rsid


###############
#forest plots##
###############
rm(meta.all)
library(ggplot2)
library(gt)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

dat <- read.delim('CHIP_AML_Apr25_new.txt', header = T)[1:16]
head(dat)
names(dat)
# change gwas names
dat$GWAS <- str_replace_all(dat$GWAS, c("CH"="CH-All", "DNMT3A"= "CH-DNMT3A", "LARGE" = "CH-Large", "SMALL" = "CH-Small",
                            "TET2" ="CH-TET2", "status"="Pan-AML", "Normal"="Normal-AML"))
# #slice table into two
# ch <- dat[c(2,3,5,6,7,10:13,15,16)]
# names(ch)
# names(ch) <- c('subtype','rsid','pos','EA','OA','OR','CI_lower','CI_upper','P','locus','gene')
# aml <- dat[c(18:20,22,26:28)]
# names(aml)
# names(aml) <- c('EA','OA','P','OR','subtype','CI_lower','CI_upper')
# aml$pos <- ch$pos
# aml$rsid <- ch$rsid
# aml$gene <- ch$gene
# aml$locus <- ch$locus

#separate lin vs kar
table(dat$study)
table(subset(dat, dat$study=='Lin_2025')$GWAS)
table(subset(dat, dat$study=='Kar_2022')$GWAS)

subset(dat, dat$study=='Lin_2025')$EA
subset(dat, dat$study=='Kar_2022')$EA

#remove complex and del5
dat <- subset(dat, !dat$GWAS %in% c('Complex', 'del5_7'))
#remove empty 
dat <- dat[-which(dat$SNP==''),]

# take one rsid at a time
snps <- unique(dat$rsid)
for (n in seq_along(snps)) {
dat.2 <- dat[dat$rsid==snps[n],]
table(dat.2$GWAS)
#dat.2 <- dat.2 %>% dplyr::arrange(desc(rsid))
dat.2$axis <- paste(dat.2$rsid, dat.2$GWAS)
dat.2$P.chr <- dat.2$P
# add colors
table(dat.2$axis)
#rbPal <- colorRampPalette(c('blue','red'))
#col.pal <- brewer.pal(name = 'Set3',n = length(unique(dat.2$rsid)))

#dat <- dat %>% mutate(Col=ifelse(SNP=='rs11212666',col.pal[1],ifelse(SNP=='rs2853677',col.pal[2],col.pal[3])))

## plot forest plot bars
len <- length(unique(dat.2$rsid))
arr.col <- rep(c("gray","black"),20)

#dat.2$P  <- dat.2[order(dat.2$P), ]

p <- dat.2 |>  ggplot(aes(y = reorder(axis, -P), col=rsid)) + 
  theme_classic() + scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.5)) +coord_cartesian(xlim = c(0.5, 1.5))

p <- p + geom_point(aes(x=OR),shape=15, size=2.5) +
  geom_linerange(aes(xmin=L95CL, xmax=U95CL)) +  geom_vline(xintercept = 1, linetype="dashed", linewidth = 0.3) +
  labs(x="Hazard Ratio") #+ facet_wrap(~rsid, strip.position='left', scales='free_y', ncol=1)

p

p_mid <- p + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank()) + guides(col="none")

p_mid 

# wrangle results into pre-plotting table form
out <- dat.2 |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(OR, L95CL, U95CL),
    ~ str_pad(
      round(.x, 3),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),# add an "-" between HR estimate confidence intervals
estimate_lab = paste0(OR, " (", L95CL, "-", U95CL, ")")) |>
  # round p-values to two decimal places, except in cases where p < .001
  mutate(P.chr = case_when(
    P.chr < .001 ~ as.character(P.chr),
    round(P.chr, 2) == .05 ~ as.character(round(P.chr,3)),
    P.chr < .01 ~ str_pad( # if less than .01, go one more decimal place
      as.character(round(P.chr, 3)),
      width = 4,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
      as.character(round(P.chr, 2)),
      width = 4,
      pad = "0",
      side = "right"
    )
  ))  |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      rsid = "Rsid",
      estimate_lab = "HR (95% CI)",
      U95CL = "conf.high",
      L95CL = "conf.low",
      P.chr = "P-value",
      GWAS = "Subtype",
      gene="Gene",
      axis= "SNP",
      study="Study",
      EA="E/A"
    )
  ) #|>
#  mutate(SNP = fct_rev(fct_relevel(SNP, "rsid")))

glimpse(out)

## plot txt
head(out)

p_left <- out |>
  ggplot(aes(y=reorder(axis, -P))) +
  geom_text(aes(x = 1, label = rsid),  hjust = 1, fontface = ifelse(out$rsid == "Rsid", "bold", "plain"))+
  geom_text(aes(x = 0, label = GWAS), hjust = 0, fontface = ifelse(out$GWAS == "Subtype", "bold", "plain"))+
  geom_text(aes(x = 0.4, label = gene), hjust = 0, fontface = ifelse(out$gene == "Gene", "bold", "plain"))+
  theme_void() + coord_cartesian(xlim = c(0, 1))

p_left

# p right
p_right <- out |> ggplot(aes(y=reorder(axis, -P))) +  geom_text(aes(x = 0.1, label = estimate_lab), hjust = 0,
                                                   fontface = ifelse(out$estimate_lab == "HR (95% CI)", "bold", "plain"))+
  geom_text(aes(x = 1, label = study), hjust=1, fontface = ifelse(out$study == "Study", "bold", "plain")) + theme_void() 

p_right <- p_right +
  geom_text(aes(x = 0, label = reorder(EA, -P)), hjust = 0.5, fontface = ifelse(out$EA == "E/A", "bold", "plain"))+
  geom_text(aes(x = 0.55, y = reorder(axis, -P), label = P.chr), hjust = 0, fontface = ifelse(out$P.chr == "P-value", "bold", "plain")) +
  coord_cartesian(xlim = c(0, 1))

p_right

# plot all
layout <- c(
  area(t = 0, l = 0, b = 100, r = 13), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 12, l = 14, b = 100, r = 22), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l =22, b = 100, r = 38) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
p_left + p_mid + p_right + plot_layout(design = layout)
#geom_segment(x= -20,y=20.5, color = "black", size=1, xend= 0,yend=20.5)

ggsave(paste0(snps[n],"_CHIP_forest.pdf"), width = 10, height = 5,units = 'in')
#new + geom_segment(x= -20,y=20.5, color = "black", size=1, xend= 0,yend=20.5)

}

#############################
### attempt to use facet_wrap
#############################
dat <- read.delim('CHIP_AML_Apr25_new.txt', header = T)[1:17]
head(dat)
names(dat)
# change gwas names
dat$GWAS <- str_replace_all(dat$GWAS, c("CH"="All-CH", "DNMT3A"= "DNMT3A-CH", "LARGE" = "Large-CH", "SMALL" = "Small-CH",
                                        "TET2" ="TET2-CH", "status"="Pan-AML", "Normal"="CN-AML"))

#separate lin vs kar
table(dat$study)
table(subset(dat, dat$study=='Lin_2025')$GWAS)
table(subset(dat, dat$study=='Kar_2022')$GWAS)

subset(dat, dat$study=='Lin_2025')$EA
subset(dat, dat$study=='Kar_2022')$EA

#change to kar.etal
dat$study <- gsub('Kar_2022','Kar et al',dat$study)
dat$study <- gsub('Lin_2025','AML GWAS',dat$study)

#remove complex and del5 , also large, small
dat <- subset(dat, !dat$GWAS %in% c('Complex', 'del5_7' , 'Large-CH', 'Small-CH'))
#remove two cond. snps "rs2356817", "rs13356700" 
dat.2 <- subset(dat, dat$rsid %in% c("rs12632224" , "rs11212666" , "rs2853677" , "rs7705526" , "rs13130545" , 
                                    "rs2086132" , "rs8088824", "rs35452836" , "rs79633204" , "rs10131341"))

#remove empty 
#dat.2 <- dat[-which(dat$SNP==''),]
table(dat.2$GWAS)
#dat.2 <- dat.2 %>% dplyr::arrange(desc(rsid))
dat.2$axis <- paste(dat.2$rsid, dat.2$GWAS)
#dat.2$P <- as.character(dat.2$P)
# add colors
table(dat.2$axis)
#rbPal <- colorRampPalette(c('blue','red'))

# use factors to plot rs2356817 and rs13356700 last
# order in the main body - 
# Rs12632224, rs11212666, Rs2853677, Rs7705526, Rs13130545, Rs2086132, Rs8088824, Rs35452836, Rs79633204, Rs10131341

dat.2$rsid <- factor(dat.2$rsid, levels = c("rs12632224" , "rs11212666" , "rs2853677" , "rs7705526" , "rs13130545" , 
                                            "rs2086132" , "rs8088824", "rs35452836" , "rs79633204" , "rs10131341"))
# second factor for gwas
unique(dat.2$GWAS)
dat.2$GWAS <- factor(dat.2$GWAS, levels = c("CN-AML","Pan-AML", "TET2-CH", "DNMT3A-CH", "All-CH"))
#dat.2$P <- factor(dat.2$P, ordered = T)
#dat.2$P <- arrange(dat.2$P)
glimpse(dat.2)

## plot forest plot bars
len <- length(unique(dat.2$rsid))
arr.col <- rep(c("gray","black"),20)

#facet_wrap(~rsid, strip.position='left',  ncol=1, scales='free_y', drop=T)
p <- dat.2 |>  ggplot(aes(y = GWAS, col=rsid)) +  theme_classic()  + 
  facet_wrap(rsid ~.,   strip.position='left', scales='free_y', shrink=T, ncol=1)

p <- p + geom_point(aes(x=OR),shape=15, size=2.5) +
  geom_linerange(aes(xmin=L95CL, xmax=U95CL)) +  geom_vline(xintercept = 1, linetype="dashed", linewidth = 0.3) + labs(x="Odds Ratio")

p

p_mid <- p + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank()) + guides(col="none")

p_mid 

# wrangle results into pre-plotting table form
out <- dat.2 |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(OR, L95CL, U95CL),
    ~ str_pad(
      round(.x, 3),
      width = 4,
      pad = "0",
      side = "right")),
  # add an "-" between HR estimate confidence intervals
  estimate_lab = paste0(OR, " (", L95CL, "-", U95CL, ")")) 


glimpse(out)

## plot txt
head(out)

p_left <- out |>
  ggplot(aes(y=GWAS)) +
  geom_text(aes(x = 0, label = GWAS), hjust = 0)+
  geom_text(aes(x = 1, label = gene), hjust = 1, fontface = "italic")+
  theme_void() + coord_cartesian(xlim = c(0, 1))

# facet_grid(rsid ~.,  space='free', scales='free') 
p_left <- p_left + facet_wrap(rsid ~.,   strip.position='left', scales='free_y', shrink=T, ncol=1) +
  theme(strip.background = element_blank(), strip.text= element_blank())

p_left

# p right
p_right <- out |> ggplot(aes(y=GWAS)) +  geom_text(aes(x = 0.35, label = estimate_lab), hjust = 0)+ theme_void() 

p_right <- p_right +
  geom_text(aes(x = 0, label = EA), hjust = 1)+
  geom_text(aes(x = 0.1, y = GWAS, label = case_con), hjust = 0)+
  geom_text(aes(x = 0.67, y = GWAS, label = P), hjust = 0) +
  geom_text(aes(x = 0.88, y = GWAS, label = study), hjust = 0) + coord_cartesian(xlim = c(0, 1))

p_right <- p_right +facet_grid(rsid ~.,  space='free', scales='free') + 
  theme(strip.background = element_blank(), strip.text= element_blank())

p_right

# plot all
layout <- c(
  area(t = 6, l = 0, b = 250, r = 8), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 5, l = 9, b = 250, r = 17), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 6, l = 18, b = 250, r = 36) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
p_left + p_mid + p_right + plot_layout(design = layout) 
#geom_segment(x= -20,y=20.5, color = "black", size=1, xend= 0,yend=20.5)

# add headers
#head.right <- ggplot(aes(y=GWAS))geom_text(aes(x = 0.2, label = 'test'), hjust = 0)+ theme_void() 
ggsave('AML_CHIPvariants_Fig4.pdf',width = 12, height = 11,units = 'in')

###########
## deleted
###########
#flip OR of AML if OA!=OA 
ch$EA==aml$EA
ch$OA==aml$OA
flip.or <- !ch$EA==aml$EA

aml$OR[flip.or] <- 1/aml$OR[flip.or]
aml$CI_lower[flip.or] <- 1/aml$CI_lower[flip.or]
aml$CI_upper[flip.or] <- 1/aml$CI_upper[flip.or]

#change Effect allele for aml
aml$EA[flip.or] <- aml$OA[flip.or]
aml$OA[flip.or] <- ch$OA[flip.or]
ch$EA==aml$EA
ch$OA==aml$OA


