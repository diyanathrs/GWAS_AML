#!/bin/bash
                            
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=02:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=16G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
#SBATCH --array=1-22                                                                                      
##SBATCH --chdir=${rootpath}/association


module purge
module load PLINK/2.00a3.7-foss-2022a 

#get the snp list from titania and use it to extract snps

plink --bfile ../merged_plink/HRC_chr.${SLURM_ARRAY_TASK_ID} --extract GWAS_r2filtered.snplist --make-bed --out chr${SLURM_ARRAY_TASK_ID}_meta_filtered

#mkdir merged_bgen
