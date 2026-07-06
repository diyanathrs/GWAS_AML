#!/usr/bin/Rscript

gwasnames <- c(paste0('NCL',seq(3,7))) 

prefix <- c('../NCL3NatComRevTrim/',
            '../NCL4newNatComRevTrim/',
            '../NCL5_ukb500kNatComRevTrim/',
            '../NCL6/',
            '../NCL7_mich/') 


for (chrom in seq(14,22)) {
	temp <- paste0(prefix,'hardcall/',gwasnames,'_chr',chrom,'_hardcalls')
  	#sink(file = paste0("chr",chrom,"mergeBAM.lst")) 
	temp2 <- data.frame(cbind(paste0(temp,".bim"),paste0(temp,".bed"),paste0(temp,".fam")))
  	write.table(temp, file = paste0('chr',chrom,'mergeBAM.lst'),quote = F,row.names = F,col.names = F,sep = '\t')

system(paste0("plink",
              " --bfile ", paste0("../NCL1_2NatComRevTrim/hardcall/NCL12_chr",chrom,"_hardcalls"),
              " --memory 32000 ",
              " --merge-list ", paste0("chr",chrom,"mergeBAM.lst"),
              " --make-bed",
	     # " --biallelic-only",
	     # " --set-missing-var-ids @:#_$r_$a",
              " --out ", paste0("merged_plink/","HRC_chr.",chrom) ) )
}
