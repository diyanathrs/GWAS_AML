#forest plots
library(readxl)
#devtools::install_github("NightingaleHealth/ggforestplot")
library(ggplot2)
library(gt)
library(tidyverse)
library(patchwork)
library(RColorBrewer)

dat <- read_xlsx('CHIP_JMA_AML_merged.xlsx', sheet = 5, .name_repair = "minimal")
dat$axis <- paste( dat$SNP, dat$GWAS)
glimpse(dat)
unique(dat$GWAS)

#flip HRs
aml  <- dat %>% dplyr::filter(GWAS %in% c( "Pan AML" ,"Normal AML")) 
ch <- dat %>% dplyr::filter(!GWAS %in% c( "Pan AML" ,"Normal AML")) 

ch$OR <- 1/ch$OR
ch$L95CL <- 1/ch$L95CL
ch$U95CL <- 1/ch$U95CL

#change Effect allele for aml
aml$EA <- aml$OA

dat <- rbind(aml,ch)
dat$P <- as.character(dat$P)
# add colors
unique(dat$SNP)
#rbPal <- colorRampPalette(c('blue','red'))
col.pal <- brewer.pal(name = 'Dark2',n = 3)

dat <- dat %>% mutate(Col=ifelse(SNP=='rs11212666',col.pal[1],ifelse(SNP=='rs2853677',col.pal[2],col.pal[3])))

## plot forest plot bars
p <- dat |>  ggplot(aes(y = axis)) + 
  theme_classic() 

p <- p + geom_point(aes(x=OR),col=dat$Col ,shape=15, size=3) +
  geom_linerange(aes(xmin=L95CL, xmax=U95CL), col=dat$Col) +  geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="Hazard Ratio") # facet_wrap(nrow = GWAS)
p

p_mid <- p + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank()) + guides(col="none")

p_mid

# wrangle results into pre-plotting table form
out <- dat |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(OR, L95CL, U95CL),
    ~ str_pad(
      round(.x, 2),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between HR estimate confidence intervals
  estimate_lab = paste0(OR, " (", L95CL, "-", U95CL, ")")) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      SNP = "SNP",
      estimate_lab = "HR (95% CI)",
      U95CL = "conf.high",
      L95CL = "conf.low",
      P = "GWAS p-value",
      GWAS = "Subtype",
      axis= "SNP",
      EA="E/A",
      Col="black"
    )
  ) #|>
#  mutate(SNP = fct_rev(fct_relevel(SNP, "rsid")))


glimpse(out)

## plot txt
head(out)

p_left <- out |>
  ggplot(aes(y=axis)) +
  geom_text(aes(x = 1, label = SNP),  hjust = 1, fontface = ifelse(out$SNP == "SNP", "bold", "plain"))+
  geom_text(aes(x = 0, label = GWAS), hjust = 0, fontface = ifelse(out$GWAS == "Subtype", "bold", "plain"))+
  theme_void() + coord_cartesian(xlim = c(0, 1))

p_left 
 

# p right
p_right <- out |> ggplot(aes(y=axis)) +  geom_text(aes(x = 0.2, label = estimate_lab), hjust = 0,
                     fontface = ifelse(out$estimate_lab == "HR (95% CI)", "bold", "plain"))+ theme_void() 

p_right <- p_right +
  geom_text(aes(x = 0, label = EA), hjust = 0.5, fontface = ifelse(out$EA == "E/A", "bold", "plain"))+
  geom_text(aes(x = 1, y = axis, label = P),
    hjust = 1,
    fontface = ifelse(out$P == "GWAS p-value", "bold", "plain")) +coord_cartesian(xlim = c(0, 1))

p_right

# plot all
layout <- c(
  area(t = 0, l = 0, b = 30, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 2.8, l = 5, b = 30, r = 10), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l =11, b = 30, r = 16) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
p_left + p_mid + p_right + plot_layout(design = layout)

#new + geom_segment(x= -20,y=20.5, color = "black", size=1, xend= 0,yend=20.5)
