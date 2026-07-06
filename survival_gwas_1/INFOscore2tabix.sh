#!/bin/bash

prefix="OxfCLL4imp"
rootpath=$(pwd)
##rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE"

declare -a nophenos=('HRCinfo')
#declare -a nophenos=('status')    
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
#$ -l h_vmem=32G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" root.dir="'${rootpath}'" ' info2tabix.R


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
#$ -pe smp 4


head -n 1 ${nophenos[i]}_${prefix}_assocLZ.res | awk -F"\t" '{print "#"\$0}' > ${nophenos[i]}_${prefix}.header
tail -n+2 ${nophenos[i]}_${prefix}_assocLZ.res | sort -T ./ -k2,2n -k3,3n | bgzip >  ${nophenos[i]}_${prefix}TMP.gz
tabix -r ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}TMP.gz > ${nophenos[i]}_${prefix}assoc_HRC.gz
tabix -s 2 -b 3 -e 3 ${nophenos[i]}_${prefix}assoc_HRC.gz

rm ${nophenos[i]}_${prefix}TMP.gz ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}_assocLZ.res
EOF

qsub -hold_jid ${jobname} ${subFile1}

sleep 1

rm ${subFile} ${subFile1}


done
