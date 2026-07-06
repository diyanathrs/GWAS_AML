#!/bin/bash                                                                                               
#SBATCH -A jamgaml                                                                                        
#SBATCH --partition=defq,long,bigmem                                                                      
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                                                                                      
#SBATCH --cpus-per-task=1                                               
#SBATCH --mem-per-cpu=16G
#SBATCH --exclude=sb017                                                                         
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --chdir=/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/association

snptest_v2.5.2 \
-data /nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6/impvRes/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz TrisomiesNCL6_chr${SLURM_ARRAY_TASK_ID}_VCFphenocov.sample \
-genotype_field GP \
-frequentist 1 -method expected \
-cov_all -hwe \
-pheno Trisomies \
-o Trisomies_NCL6_chr${SLURM_ARRAY_TASK_ID}_HRCRes.out.gz

# change to ResultSummary folder

cd ../ResultSummary/

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args root.dir="'/nobackup/proj/jamgaml/dean_AML/HRCimpvData/NCL6'" dataprefix="'NCL6'" chrom="'${SLURM_ARRAY_TASK_ID}'" phenoprefix="'Trisomies'" '  /nobackup/proj/jamgaml/dean_AML/HRCimpvData/res_summary.R


