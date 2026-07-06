#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk


module load R/3.3.1-intel-2017.03-GCC-6.3


R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'status'" ' metaRes2tabix.R

module load HTSlib/1.4.1-intel-2017.03-GCC-6.3

head -n 1 status_PCspeAMLHRC_assocLZ.res | awk '{for (i=1;i<=NF;i++) print $i}' | paste -s | awk -F"\t" '{print "#"$0}' > status_PCspeAMLHRC.header
tail -n+2 status_PCspeAMLHRC_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  status_PCspeAMLHRCTMP.gz
tabix -r status_PCspeAMLHRC.header status_PCspeAMLHRCTMP.gz > status_PCspeAMLHRCassoc_HRC.gz
tabix -s 1 -b 2 -e 2 status_PCspeAMLHRCassoc_HRC.gz

rm status_PCspeAMLHRCTMP.gz status_PCspeAMLHRC.header status_PCspeAMLHRC_assocLZ.res

