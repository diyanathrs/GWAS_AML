#!/bin/bash
rootpath=`pwd`

subFile1="unzipping.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash                                                                     
#SBATCH -A jamgaml                                                              
#SBATCH --partition=short                                                        
#SBATCH --time=00:10:00                                                         
#SBATCH --nodes=1                                                               
#SBATCH --ntasks=1                                                              
#SBATCH --cpus-per-task=1                                                       
#SBATCH --mem-per-cpu=4G                                                       
#SBATCH --mail-type=FAIL                                                        
#SBATCH --array=1-22                                                            
#SBATCH --workdir=${rootpath}

unzip -P 'n56SmtqBDbRI6Y' chr_\${SLURM_ARRAY_TASK_ID}.zip

###7za x chr_\${SGE_TASK_ID}.zip -p'SUFzGdR0>m5tKx'

EOF

sbatch ${subFile1}

exit



subFile1="checking.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash                                                                     
#SBATCH -A jamgaml                                                              
#SBATCH --partition=short                                                        
#SBATCH --time=00:10:00                                                         
#SBATCH --nodes=1                                                               
#SBATCH --ntasks=1                                                              
#SBATCH --cpus-per-task=1                                                       
#SBATCH --mem-per-cpu=4G                                                       
#SBATCH --mail-type=FAIL                                                        
#SBATCH --array=1-22                                                            
#SBATCH --workdir=${rootpath}

md5sum chr_\${SLURM_ARRAY_TASK_ID}.zip >> mycheckmd5sum



EOF

sbatch ${subFile1}

