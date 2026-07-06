#!/bin/bash

#SBATCH --partition=defq,long                                                                  
#SBATCH --time=01:00:00                                                                        
#SBATCH --nodes=1                                                                              
#SBATCH --ntasks=1                                                                             
#SBATCH --cpus-per-task=1                                                                      
#SBATCH --mem-per-cpu=16G                                                                      
#SBATCH --export=ALL                                                                           
#SBATCH --mail-type=FAIL                                                                       
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                            
##SBATCH --array=1-22                                                                           
##SBATCH --chdir=${rootpath}/association


gwas="NCL6"

mkdir updated_hardcalls
module load R/4.2.1-foss-2022a
R CMD BATCH --vanilla --no-timing '--args gwas="'${gwas}'" ' ../../change_snpIDs.R

mv ../updated_hardcalls/*.bim .
