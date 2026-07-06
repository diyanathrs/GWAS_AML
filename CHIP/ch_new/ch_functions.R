#ch functions 
# isolate top hits
aml.gwas <- lapply(aml.gwas, function(x){
  n.max <- max(x$N)
  x.sig <- x %>% filter(P < 5e-8 & N == n.max )
  return(x.sig) 
})

#annotate all
aml.gwas <- lapply(aml.gwas, function(x){
  x$rsid <- paste0(x$CHR,':',x$BP,'_',x$A2,'_', x$A1)
  x$snp.id <- paste0(x$CHR,':',x$BP)
  #merge
  x <- merge(x, hrc.rsid[,c(3,6)], by='rsid')
  return(x)
})
 