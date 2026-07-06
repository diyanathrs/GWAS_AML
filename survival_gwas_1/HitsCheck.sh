#!/bin/bash

newdir=Hits_checkAll
[ -d ${newdir} ] || mkdir ${newdir}

rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLL"
prefix="CLL_finalQced4imp"

# no. phenos
#declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')

declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""

# hits screening
for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${prefix}_hitCheck.sge"
    cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 
#$ -wd ${rootpath}/Hits_checkAll


R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${rootpath}/hits_check.R

EOF

    jobname="hitcheck_${nophenos[i]}_${prefix}"
    qsub -N ${jobname} ${subFile}

done

# hits output
# 
for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile1="${nophenos[i]}_${prefix}_hitall.sge"
    cat > ${subFile1} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 
#$ -wd ${rootpath}/Hits_checkAll

R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${rootpath}/hits_all.R

EOF

    jobname1="hitsall_${nophenos[i]}_${prefix}"
    qsub -hold_jid "hitcheck_*" -N ${jobname1} ${subFile1}  

done

# subFile2=${prefix}_hitsPlots.sge
# cat > ${subFile2} << EOF
# #!/bin/bash
# #$ -V
# #$ -l h_vmem=16G
# #$ -R y
# #$ -j yes
# #$ -l h_rt=168:00:00   
# #$ -m a
# #$ -M Wei-Yu.Lin@newcastle.ac.uk 
# #$ -wd ${rootpath}/Hits_checkAll

# module load apps/R/3.2.3

# R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" '  /home/nwl15/WORKING_DATA/HRCimpvData/hits_all_summaryHRC.R

# EOF

# qsub -hold_jid "hitsall_*" ${subFile2}
