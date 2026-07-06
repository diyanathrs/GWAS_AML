#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=03:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=16G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
#SBATCH --array=1-1                                                                                      
##SBATCH --chdir=${rootpath}/association


rootpath=$(pwd)

#impvpath="${rootpath}/../impvRes"
#gwas="NCL12"

module purge
module load PLINK/2.00a3.7-foss-2022a 

#mkdir merged_bgen

# keep only controls not cases
plink2 --bfile HRC_chr${SLURM_ARRAY_TASK_ID} --export bgen-1.2 --out merged_bgen/HRC_merged_chr${SLURM_ARRAY_TASK_ID} 

cd merged_bgen
../../../finemapping/bgen_v1.1.4-CentOS6.8-x86_64/bgenix -g HRC_merged_chr${SLURM_ARRAY_TASK_ID}.bgen -index                                                              
