#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 studyName"; exit; }

folderlist=CLL.lst
#folderlist=CLL_ELN2noSCT.lst
nodata=(`cat ${folderlist} | awk '{print $1}'`)
echo ${#nodata[@]}

rootpath=$(pwd)

# no. phenos
#nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
#declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')    

## declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')
## declare -a nophenos=('OS_Dx_Death_LFU_Status')

#declare -a nophenos=('OSstatus' 'RFSstatus')
declare -a nophenos=('OSstatus' )

nophenos=("${nophenos[@]/#/$1_}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}

for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_resformat.slurm"
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk
#SBATCH --array=1-${#nodata[@]}
#SBATCH --chdir=${rootpath}

module load R/3.4.3-foss-2017b-X11-20171023

resdir=\$(cat ${folderlist} | awk -v kkk=\${SLURM_ARRAY_TASK_ID} 'NR==kkk{print \$1}')
echo \${resdir}
prefix=\$(cat ${folderlist} | awk -v kkk=\${SLURM_ARRAY_TASK_ID} 'NR==kkk{print \$2}')
echo \${prefix}


R CMD BATCH --vanilla --no-timing '--args dataprefix="'\${prefix}'" phenoprefix="'${nophenos[i]}'" res.dir="'\${resdir}'" ' Surassoc2tabix.R


EOF

command1="sbatch ${subFile}"
jid1=$($command1 | awk ' { print $4 }')


subFile1="${nophenos[i]}_assoc2tabix.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk
#SBATCH --array=1-${#nodata[@]}
#SBATCH --chdir=${rootpath}

module load HTSlib/1.4.1-intel-2017.03-GCC-6.3

prefix=\$(cat ${folderlist} | awk -v kkk=\${SLURM_ARRAY_TASK_ID} 'NR==kkk{print \$2}')
echo \${prefix}

head -n 1 ${nophenos[i]}_\${prefix}_assocLZ.res | awk -F"\t" '{print "#"\$0}' > ${nophenos[i]}_\${prefix}.header
tail -n+2 ${nophenos[i]}_\${prefix}_assocLZ.res | sort -T ./ -k2,2n -k3,3n | bgzip >  ${nophenos[i]}_\${prefix}TMP.gz
tabix -r ${nophenos[i]}_\${prefix}.header ${nophenos[i]}_\${prefix}TMP.gz > ${nophenos[i]}_\${prefix}assoc_HRC.gz
tabix -s 2 -b 3 -e 3 ${nophenos[i]}_\${prefix}assoc_HRC.gz

rm ${nophenos[i]}_\${prefix}TMP.gz ${nophenos[i]}_\${prefix}.header ${nophenos[i]}_\${prefix}_assocLZ.res
EOF


command2="sbatch --dependency=afterok:${jid1} ${subFile1}"
jid2=$($command2 | awk ' { print $4 }')


rm ${subFile} ${subFile1}


done
