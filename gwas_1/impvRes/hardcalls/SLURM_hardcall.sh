#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=08:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=16G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
#SBATCH --array=1-15                                                                                      
##SBATCH --chdir=${rootpath}/association


rootpath=$(pwd)

#impvpath="${rootpath}/../../impvRes"
gwas="NCL6"

module purge
#module load QCTOOL/v2.0-rc8-foss-2019a
#module load BCFtools/1.3-foss-2017b
                                                              
#qctool --vcf ${impvpath}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz dosage=GP --import-dosage-certainty 0.8 --make-bed --out ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls

#qctool -g ${impvpath}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -threshold 0.8 -og ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.bgen

#qctool -g ${impvpath}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -threshold 0.8 -snp-stats -osnp ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.txt

#bcftools view --force-samples --samples-file AML_cons_final.txt ../chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -O z -o ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.vcf.gz

module load QCTOOL/v2.0-rc8-foss-2019a

qctool -g ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.vcf.gz -threshold 0.8 -og ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.bgen -os ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls.sample
