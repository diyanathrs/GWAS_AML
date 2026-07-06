#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk


module load R/4.2.1-foss-2022a


R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'t15_17'" ' metaRes2tabix.R

module load HTSlib/1.15.1-GCC-11.3.0

head -n 1 t15_17_PCspeAMLHRC_assocLZ.res | awk '{for (i=1;i<=NF;i++) print $i}' | paste -s | awk -F"\t" '{print "#"$0}' > t15_17_PCspeAMLHRC.header
tail -n+2 t15_17_PCspeAMLHRC_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  t15_17_PCspeAMLHRCTMP.gz
tabix -r t15_17_PCspeAMLHRC.header t15_17_PCspeAMLHRCTMP.gz > t15_17_PCspeAMLHRCassoc_HRC.gz
tabix -s 1 -b 2 -e 2 t15_17_PCspeAMLHRCassoc_HRC.gz

rm t15_17_PCspeAMLHRCTMP.gz t15_17_PCspeAMLHRC.header t15_17_PCspeAMLHRC_assocLZ.res

