#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=ln03,sb017
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load R/4.2.1-foss-2022a

R CMD BATCH --vanilla --no-timing '--args dataprefix="'PCspeAMLHRC'" phenoprefix="'CBF'" '  meta_summaryQQMH.R

# different maf by dean
#R CMD BATCH --vanilla --no-timing '--args dataprefix="''" phenoprefix="'CBF'" '  meta_summaryQQMH_maf.R



