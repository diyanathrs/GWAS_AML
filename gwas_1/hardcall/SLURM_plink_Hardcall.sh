#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=05:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=16G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
#SBATCH --array=1-22                                                                                      
##SBATCH --chdir=${rootpath}/association


rootpath=$(pwd)

impvpath="${rootpath}/../impvRes"
gwas="NCL6"

module purge
module load PLINK/2.00a3.7-foss-2022a 
                                                              
plink2 --vcf ${impvpath}/chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz dosage=GP --import-dosage-certainty 0.8 --record vcf --out ${gwas}_chr${SLURM_ARRAY_TASK_ID}_hardcalls

