args=(commandArgs(TRUE))
if(length(args)==0){
     print("No arguments supplied.")

 } else {
     for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
     }
 }


chrom

gwasnames <- c(paste0('NCL',seq(3,7)))

prefix <- c('../NCL3NatComRevTrim/impvRes/hardcalls/',
             '../NCL4newNatComRevTrim/impvRes/hardcalls/',
             '../NCL5_ukb500kNatComRevTrim/impvRes/hardcalls/',
             '../NCL6/impvRes/hardcalls/',
             '../NCL7_mich/impvRes2/hardcalls/')


#for (chrom in seq(1,15)) {
       #temp <- paste0(prefix,'hardcalls/',gwasnames,'_chr',chrom,'_hardcalls')
       #sink(file = paste0("chr",chrom,"mergeBAM.lst"))
       #temp2 <- data.frame(cbind(paste0(temp,".bim"),paste0(temp,".bed"),paste0(temp,".fam")))
       #write.table(temp, file = paste0('chr',chrom,'mergeBAM.lst'),quote = F,row.names = F,col.names = F,sep = '\t')

cmd <- paste0("qctool",
               " -g ", paste0("../NCL1_2NatComRevTrim/impvRes/hardcalls/NCL12_chr",chrom,"_hardcalls.bgen"),
	       " -s ", paste0("../NCL1_2NatComRevTrim/impvRes/hardcalls/NCL12_chr",chrom,"_hardcalls.sample"),
               " -g ", paste0(prefix[1], gwasnames[1],"_chr", chrom,"_hardcalls.bgen"),
	       " -s ", paste0(prefix[1], gwasnames[1],"_chr", chrom,"_hardcalls.sample"),
               " -g ", paste0(prefix[2], gwasnames[2],"_chr", chrom,"_hardcalls.bgen"),
	       " -s ", paste0(prefix[2], gwasnames[2],"_chr", chrom,"_hardcalls.sample"),
               " -g ", paste0(prefix[3], gwasnames[3],"_chr", chrom,"_hardcalls.bgen"),
  	       " -s ", paste0(prefix[3], gwasnames[3],"_chr", chrom,"_hardcalls.sample"),
               " -g ", paste0(prefix[4], gwasnames[4],"_chr", chrom,"_hardcalls.bgen"),
	       " -s ", paste0(prefix[4], gwasnames[4],"_chr", chrom,"_hardcalls.sample"),
               " -g ", paste0(prefix[5], gwasnames[5],"_chr", chrom,"_hardcalls.bgen"),
	       " -s ", paste0(prefix[5], gwasnames[5],"_chr", chrom,"_hardcalls.sample"),

               " -og ", paste0("merged_bgen_new/chr", chrom, "_merged.bgen"),
  	       " -os ", paste0("merged_bgen_new/chr", chrom, "_merged.sample"),
               " -compare-variants-by position "
)

print(cmd)

system(cmd)
