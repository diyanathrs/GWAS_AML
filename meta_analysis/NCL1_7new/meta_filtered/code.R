##subset N=max(N) from the meta outputs

require(data.table)
files <- list.files(path = '../',pattern = 'meta.gz', full.names = T)

for (i in seq_along(files)){
  dat <- setDT(read.table(files[i], header = T, check.names = F))
  nmax <- max(dat$N)
  dat <- dat[N==nmax]
  out.name <- gsub(x = files[i], pattern = '.meta.gz',replacement = '')
  out.name <- gsub(x = out.name, pattern = '..//', replacement = '')
  print(out.name)
  fwrite(dat, file = paste0(out.name,'_filtered','.meta') ,quote = F, sep = '\t' ,row.names = F)
}

print("Done")
