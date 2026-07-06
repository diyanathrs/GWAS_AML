#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=ln03,sb017
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'status'" folderlist="'NCLpcspeData1_5.lst'"'  meta_summaryHits.R 

R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'status'" '  forestplot_addMAF.R

