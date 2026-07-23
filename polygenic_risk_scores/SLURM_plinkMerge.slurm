#!/bin/bash
#SBATCH --account=rockhpc_jamgaml
#SBATCH --partition=default_free                                                                                               
#SBATCH --time=10:00:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1                                                                                 
#SBATCH --mem-per-cpu=32G                                                                                  
#SBATCH --export=ALL                                                                                      
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                                                                                                                 


rootpath=$(pwd)

#impvpath="${rootpath}/../impvRes"
#gwas="NCL12"

module purge
module load R
#module load plink/2.0.0
module load PLINK/1.9b_6.21-x86_64

Rscript plink_merge.R                                                              
rm *mergeBAM.lst
