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


R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'del5_7'" ' metaRes2tabix.R

module load HTSlib/1.15.1-GCC-11.3.0

head -n 1 del5_7_PCspeAMLHRC_assocLZ.res | awk '{for (i=1;i<=NF;i++) print $i}' | paste -s | awk -F"\t" '{print "#"$0}' > del5_7_PCspeAMLHRC.header
tail -n+2 del5_7_PCspeAMLHRC_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  del5_7_PCspeAMLHRCTMP.gz
tabix -r del5_7_PCspeAMLHRC.header del5_7_PCspeAMLHRCTMP.gz > del5_7_PCspeAMLHRCassoc_HRC.gz
tabix -s 1 -b 2 -e 2 del5_7_PCspeAMLHRCassoc_HRC.gz

rm del5_7_PCspeAMLHRCTMP.gz del5_7_PCspeAMLHRC.header del5_7_PCspeAMLHRC_assocLZ.res

