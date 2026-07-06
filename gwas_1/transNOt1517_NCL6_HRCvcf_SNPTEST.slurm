#!/bin/bash                                                                                               
#SBATCH -A jamgaml                                                                                        
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                                                                                      
#SBATCH --cpus-per-task=1                                               
#SBATCH --mem-per-cpu=16G
#SBATCH --exclude=sb017                                                                         
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --workdir=/nobackup/proj/jamgaml/HRCimpvData/NCL6/association

snptest_v2.5.2 \
-data /nobackup/proj/jamgaml/HRCimpvData/NCL6/impvRes/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz transNOt1517NCL6_chr${SLURM_ARRAY_TASK_ID}_VCFphenocov.sample \
-genotype_field GP \
-frequentist 1 -method expected \
-cov_all -hwe \
-pheno transNOt1517 \
-o transNOt1517_NCL6_chr${SLURM_ARRAY_TASK_ID}_HRCRes.out.gz

# change to ResultSummary folder

cd ../ResultSummary/

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args root.dir="'/nobackup/proj/jamgaml/HRCimpvData/NCL6'" dataprefix="'NCL6'" chrom="'${SLURM_ARRAY_TASK_ID}'" phenoprefix="'transNOt1517'" '  /nobackup/proj/jamgaml/HRCimpvData/res_summary.R


