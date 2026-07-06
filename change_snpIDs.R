args=(commandArgs(TRUE)) 
if(length(args)==0){
    print("No arguments supplied.")
    
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    } 
}

gwas

#gwas <- 'NCL12' 
bim.lst <- list.files(path="./", pattern="*.bim")

bim.prefix <- gsub('_chr*_hardcalls','',bim.lst) 
	for (chrom in seq(1, 22)){ 
	dat <- read.delim(paste0(gwas,'_chr',chrom,'_hardcalls.bim'),header = F) 
	print(paste0('Loading ',gwas,'_chr',chrom,'_hardcalls.bim')) 
	snp.lst <- paste0(dat$V1,':',dat$V4,'_',dat$V6,'_',dat$V5) 
	dat$V2 <- snp.lst
	#snp.lst2 <- data.frame(dat$V2,snp.lst) 
	write.table(dat,file = paste0('updated_hardcalls/',gwas,'_chr',chrom,'_hardcalls','.bim'),quote = F,sep = '\t',col.names = F,row.names = F)
	}

print('All done !')
