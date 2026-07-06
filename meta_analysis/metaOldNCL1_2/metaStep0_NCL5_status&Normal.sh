#!/bin/bash

folderlist=NCLpcspeData.lst
nodata=(`cat ${folderlist} | awk '{print $1}'`)
echo ${#nodata[@]}

rootpath=$(pwd)

# no. phenos
##nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
##declare -a nophenos=('panAMLnoAPL')
declare -a nophenos=('status' 'Normal')

#declare -a nophenos=( 'Trans' 'Complex' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')


echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_metaformating.slurm"
    cat > ${subFile} << EOF
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
#SBATCH --array=1-${#nodata[@]}
#SBATCH --chdir=${rootpath}

module load R/4.2.1-foss-2022a

resdir=\$(cat ${folderlist} | awk -v kkk=\${SLURM_ARRAY_TASK_ID} 'NR==kkk{print \$1}')
echo \${resdir}
prefix=\$(cat ${folderlist} | awk -v kkk=\${SLURM_ARRAY_TASK_ID} 'NR==kkk{print \$2}')
echo \${prefix}


R CMD BATCH --vanilla --no-timing '--args dataprefix="'\${prefix}'" phenoprefix="'${nophenos[i]}'" res.dir="'\${resdir}'" maf="'0.02'" ' metaNCL_gwas_maf.R


EOF

sbatch ${subFile}
sleep 1

rm ${subFile}

done

