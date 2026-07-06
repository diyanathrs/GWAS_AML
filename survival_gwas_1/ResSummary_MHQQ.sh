#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 studyName"; exit; }

newdir=ResultSummary
[ -d ${newdir} ] || mkdir ${newdir}

rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE"
#impvpath="${rootpath}/impvRes"
prefix="CLL_OEEfinalQced4imp"


# no. phenos
#declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')

declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')

#declare -a nophenos=('OS_Dx_Death_LFU_Status')


nophenos=("${nophenos[@]/#/$1_}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""

echo ${nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    # tidying up the results based on maf and info
    subFile="${nophenos[i]}_${prefix}_HRCtidyup.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=10G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -t 1-22
#$ -wd ${rootpath}/ResultSummary


R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SGE_TASK_ID}'"  phenoprefix="'${nophenos[i]}'" '  ${rootpath}/res_summarySur.R

EOF

    jobname="${nophenos[i]}_${prefix}tidyup"
    qsub -N ${jobname} ${subFile}

    # QQ and MH plots
    subFile1="${nophenos[i]}_${prefix}_QQMH.sge"
    cat > ${subFile1} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -wd ${rootpath}/ResultSummary

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${rootpath}/res_summaryQQMHSur.R

EOF

qsub -hold_jid ${jobname} ${subFile1}
sleep 1

done

