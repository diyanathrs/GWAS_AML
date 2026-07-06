#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=ln03,sb017
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load R/4.2.1-foss-2022a

#R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'t15_17'" folderlist="'../NCLpcspeData.lst'"'  ../meta_summaryHits.R 

# added new code to include NCL6 & 7 in grepl in forestplot_addMAF.R
R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'t15_17'" '  ../forestplot_addMAF.R

