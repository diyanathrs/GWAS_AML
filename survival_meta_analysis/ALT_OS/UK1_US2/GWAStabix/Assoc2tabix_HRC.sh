#!/bin/bash

rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLLbystudyFinalInfo90"
prefix="CLL_finalQced4imp"

#folderlist=CLL.lst
#nodata=(`cat ${folderlist} | awk '{print $1}'`)
#echo ${#nodata[@]}

condit="MAF25INFO90meta"

# no. phenos

declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')
## declare -a nophenos=('OS_Dx_Death_LFU_Status')
nophenos=("${nophenos[@]/#/$1}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_resformat.sge"
    cat > ${subFile} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=8G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a

###$ -t 1-${#nodata[@]}
###resdir=\$(cat ${folderlist} | awk -v kkk=\${SGE_TASK_ID} 'NR==kkk{print \$1}')
###echo \${resdir}
###prefix=\$(cat ${folderlist} | awk -v kkk=\${SGE_TASK_ID} 'NR==kkk{print \$2}')
###echo \${prefix}


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" root.dir="'${rootpath}'" ' plinkmeta_summary.R


EOF

jobname="${nophenos[i]}_resOUT"
qsub -N ${jobname} ${subFile}

subFile1="${nophenos[i]}_assoc2tabix.sge"
cat > ${subFile1} <<EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=6G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -m a
#$ -pe smp 2

##$ -t 1-${#nodata[@]}
##
##prefix=\$(cat ${folderlist} | awk -v kkk=\${SGE_TASK_ID} 'NR==kkk{print \$2}')
##echo \${prefix}

head -n 1 ${nophenos[i]}_${prefix}_assocLZ.res | awk -F"\t" '{print "#"\$0}' > ${nophenos[i]}_${prefix}.header
tail -n+2 ${nophenos[i]}_${prefix}_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  ${nophenos[i]}_${prefix}TMP.gz
tabix -r ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}TMP.gz > ${nophenos[i]}_${prefix}${condit}assoc_HRC.gz
tabix -s 1 -b 2 -e 2 ${nophenos[i]}_${prefix}${condit}assoc_HRC.gz

rm ${nophenos[i]}_${prefix}TMP.gz ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}_assocLZ.res
EOF

qsub -hold_jid ${jobname} ${subFile1}
sleep 1

rm ${subFile} ${subFile1}


done
