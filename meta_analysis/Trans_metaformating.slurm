#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem                                                                                  
#SBATCH --time=00:30:00                                                                                   
#SBATCH --nodes=1                                                                                         
#SBATCH --ntasks=1                                                                                        
#SBATCH --cpus-per-task=1
#SBATCH --exclude=ln03,sb017                                                                    
#SBATCH --mem-per-cpu=12G                                                                                 
#SBATCH --mail-type=FAIL                                                                                  
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk                                                            
#SBATCH --array=1-6
#SBATCH --chdir=/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes

module load R/4.2.1-foss-2022a
#print folderlist
echo NCLpcspeData_no6.lst

# removing N<50 GWAS for phenos
resdir=$(cat NCLpcspeData_no6.lst | awk -v kkk=${SLURM_ARRAY_TASK_ID} 'NR==kkk{print $1}')
echo ${resdir}
prefix=$(cat NCLpcspeData_no6.lst | awk -v kkk=${SLURM_ARRAY_TASK_ID} 'NR==kkk{print $2}')
echo ${prefix}


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'Trans'" res.dir="'${resdir}'" maf="'0.02'" ' metaNCL_gwas_maf.R


