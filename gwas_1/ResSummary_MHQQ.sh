#!/bin/bash

newdir=ResultSummary
[ -d ${newdir} ] || mkdir ${newdir}

rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/NCL5_ukb500k"
#impvpath="${rootpath}/impvRes"
prefix="NCL5_finalQced4imp"

# no. phenos
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""

for ((i=0; i<${#nophenos[@]}; i++))
do
    # tidying up the results based on maf and info
    subFile="${nophenos[i]}_${prefix}_HRCtidyup.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'  
#$ -m a
#$ -t 1-22
#$ -wd ${rootpath}/ResultSummary

## module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SGE_TASK_ID}'"  phenoprefix="'${nophenos[i]}'" '  /home/nwl15/WORKING_DATA/HRCimpvData/res_summary.R

EOF

    jobname="${nophenos[i]}_${prefix}tidyup"
    qsub -N ${jobname} ${subFile}

    # QQ and MH plots
    subFile1="${nophenos[i]}_${prefix}_QQMH.sge"
    cat > ${subFile1} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'
#$ -m a
#$ -wd ${rootpath}/ResultSummary

## module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  /home/nwl15/WORKING_DATA/HRCimpvData/res_summaryQQMH.R

EOF

qsub -hold_jid ${jobname} ${subFile1}
sleep 1

done

