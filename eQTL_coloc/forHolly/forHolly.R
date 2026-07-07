# collate hollys eqtl snps
snps <- c("rs28645857", "rs115692085" , "rs3791334" , "rs6450183" , "rs16829165")
#load files
files <- list.files(pattern = '*R4.datainplots.SuppTable2.txt')
eqtl <- lapply(files, read.table, header = T)
head(eqtl)
eqtl <- rbindlist(eqtl)

stopifnot(snps %in% eqtl$SNPrsID)
eqtlAML <- eqtl[eqtl$SNPrsID %in% snps]
names(eqtlAML)
eqtlAML <- eqtlAML[order(eqtlAML$eQTLPvalue)]
write.csv(eqtlAML, paste0('eQTLgenforHolly_',format(Sys.time(),'%m%b%y'),'.csv'), quote = F, row.names = F)
