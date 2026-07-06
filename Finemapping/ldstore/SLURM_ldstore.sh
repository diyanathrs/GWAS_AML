#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=4:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=64G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
##SBATCH --array=11-14                                                                                      
##SBATCH --chdir=${rootpath}/association


rootpath=$(pwd)

#impvpath="${rootpath}/../impvRes"
#gwas="NCL12"
ldstore.exec="/mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64"

module purge
module load R 

#mkdir merged_bgen

# keep only controls not cases
#plink2 --bfile HRC_chr${SLURM_ARRAY_TASK_ID} --export bgen-1.2 --out merged_bgen/HRC_merged_chr${SLURM_ARRAY_TASK_ID} 

#Rscript Run_ldstore.R chrom=${SLURM_ARRAY_TASK_ID}                                                             

# run ldstore directly
/mnt/storage/nobackup/proj/jamgaml/dean_AML/finemapping/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64  --in-files HRC_all.master  --read-only-bgen  --write-bcor 
