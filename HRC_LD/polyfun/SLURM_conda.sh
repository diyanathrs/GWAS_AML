#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=10:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=32G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
##SBATCH --array=1-12                                                                                      
##SBATCH --chdir=${rootpath}/association


rootpath=$(pwd)

#impvpath="${rootpath}/../impvRes"
#gwas="NCL12"


# keep only controls not cases
#plink2 --bfile HRC_chr${SLURM_ARRAY_TASK_ID} --export bgen-1.2 --out merged_bgen/HRC_merged_chr${SLURM_ARRAY_TASK_ID} 

#R CMD BATCH --vanilla --no-timing '--args chrom="'${SLURM_ARRAY_TASK_ID}'" ' qctool_merge.R                                                              

conda env create -f polyfun.yml
