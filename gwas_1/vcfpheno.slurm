#!/bin/bash                                                                                             
#SBATCH -A jamgaml                                                                                      
#SBATCH --partition=defq,long,bigmem,short                                                                
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                                                                                      
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G                                                                                
#SBATCH --exclude=ln03,sb017
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#####SBATCH --array=1-22
#SBATCH --chdir=/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/association


module load R/3.3.1-intel-2017.03-GCC-6.3

###R CMD BATCH --vanilla --no-timing '--args impv.path="'/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/impvRes'" prefix="'NCL6'" chrom="'${SLURM_ARRAY_TASK_ID}'"  pheno="'/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/NCL6_phenotypeR1_12Feb2019.txt'" cov="'NCL6_sampleCleanPC.cov'" pcnum="''" '  /nobackup/proj/jamgaml/dean_AML/HRCimpvData/phenoPrepVCF.R

R CMD BATCH --vanilla --no-timing '--args impv.path="'/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/impvRes'" prefix="'NCL6'" pheno="'/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/NCL6_phenotypeR1_12Feb2019.txt'" cov="'NCL6_sampleCleanPC.cov'"  phenoname="'Trisomies'"' /nobackup/proj/jamgaml/dean_AML/HRCimpvData/phenoPrepVCFnew_dean.R


