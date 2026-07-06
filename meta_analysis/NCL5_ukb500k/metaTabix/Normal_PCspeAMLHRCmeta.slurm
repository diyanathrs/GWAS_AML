#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL


module load R/3.3.1-intel-2017.03-GCC-6.3


R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'Normal'" ' metaRes2tabix.R

module load HTSlib/1.4.1-intel-2017.03-GCC-6.3

head -n 1 Normal_PCspeAMLHRC_assocLZ.res | awk '{for (i=1;i<=NF;i++) print $i}' | paste -s | awk -F"\t" '{print "#"$0}' > Normal_PCspeAMLHRC.header
tail -n+2 Normal_PCspeAMLHRC_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  Normal_PCspeAMLHRCTMP.gz
tabix -r Normal_PCspeAMLHRC.header Normal_PCspeAMLHRCTMP.gz > Normal_PCspeAMLHRCassoc_HRC.gz
tabix -s 1 -b 2 -e 2 Normal_PCspeAMLHRCassoc_HRC.gz

rm Normal_PCspeAMLHRCTMP.gz Normal_PCspeAMLHRC.header Normal_PCspeAMLHRC_assocLZ.res

