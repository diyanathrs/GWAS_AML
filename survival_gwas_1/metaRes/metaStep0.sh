#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 studyName"; exit; }

# folderlist=NCL1_4_folder.lst
# nodata=(`cat ${folderlist} | awk '{print $1}'`)
# echo ${#nodata[@]}
# no. phenos
# nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype') 

declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')

nophenos=("${nophenos[@]/#/$1_}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}

prefix="CLL_OEEfinalQced4imp"
resdir="/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE/ResultSummary"

for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_metaformating.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=8G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 

##$ -t 1-${#nodata[@]}
## module load apps/R/3.2.3

# resdir=\$(cat ${folderlist} | awk -v kkk=\${SGE_TASK_ID} 'NR==kkk{print \$1}')
# echo \${resdir}
# prefix=\$(cat ${folderlist} | awk -v kkk=\${SGE_TASK_ID} 'NR==kkk{print \$2}')
# echo \${prefix}


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" res.dir="'${resdir}'" ' metaNCL_gwas.R


EOF

qsub ${subFile}
sleep 1

rm ${subFile}

done




# qsub ${subFile}