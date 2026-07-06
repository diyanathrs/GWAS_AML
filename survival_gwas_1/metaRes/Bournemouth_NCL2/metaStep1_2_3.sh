#!/bin/bash

newdir=Hitsmeta
[ -d ${newdir} ] || mkdir ${newdir}

newdir=MH_QQplots
[ -d ${newdir} ] || mkdir ${newdir}


ncl4folder="CLL_finalQced4imp"

# no. phenos
declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')   
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${ncl4folder}meta.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=16G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a

plink --meta-analysis *${nophenos[i]}_${ncl4folder}_assoc.res + logscale report-all \\
--out ${nophenos[i]}_NCL_${ncl4folder}


EOF

jobname="meta_${nophenos[i]}_NCL_${ncl4folder}"
qsub -N ${jobname} ${subFile}

# qq AND MH plots

subFile1="${nophenos[i]}_${ncl4folder}meta_mhqq.sge"
    
    cat > ${subFile1} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=10G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a

module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryQQMH.R


EOF

qsub -hold_jid ${jobname} ${subFile1}

#sleep 1

# hits & forest plots

subFile2="${nophenos[i]}_${ncl4folder}meta_hits.sge"
    
    cat > ${subFile2} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=16G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a

##module load apps/R/3.2.3


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryHits.R 

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  forestplot_addMAF.R

EOF


qsub -hold_jid ${jobname} ${subFile2}

sleep 1


done





