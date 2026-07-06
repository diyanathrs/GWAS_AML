#!/bin/bash
#SBATCH -A jamgaml                                                                                      
#SBATCH --partition=defq,long,bigmem                                                                    
#SBATCH --time=48:00:00                                                                                 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --exclude=ln03,sb017
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --workdir=/nobackup/proj/jamgaml/HRCimpvData/NCL6/ResultSummary

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'NCL6'" phenoprefix="'Any.Monosomy'" '  /nobackup/proj/jamgaml/HRCimpvData/res_summaryQQMH.R

