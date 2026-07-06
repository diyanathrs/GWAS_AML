#!/bin/bash

rootpath=/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE

# change this path to where your gen and sample files
# datadir="/home/nwl15/WORKING_DATA/imputeTest/TRsurv"
datadir=${rootpath}/impvRes/mytestvcfs

# change file name of your master accordingly
# phenofile=CLL_phenoUpdate_02May2017.txt

subFile="sangerExtract.sge"
cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=16G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -t 6,10
#$ -wd ${rootpath}/SangerSNPs


qctool -g ${datadir}/CLL_OEEfinalQced4imp_chr\${SGE_TASK_ID}subset.gen.gz \\
-s ${datadir}/CLL_OEEfinalQced4imp_chr\${SGE_TASK_ID}subset.samples \\
-og sangerChr\${SGE_TASK_ID}subset.gen -os sangerChr\${SGE_TASK_ID}subset.sample \\
-incl-rsids sangerSNPid.lst 


plink --gen ${datadir}/CLL_OEEfinalQced4imp_chr\${SGE_TASK_ID}subset.gen.gz \\
--memory 16000 \\
--sample ${datadir}/CLL_OEEfinalQced4imp_chr\${SGE_TASK_ID}subset.samples \\
--extract sangerSNPid.lst \\
--recode A \\
--out chr\${SGE_TASK_ID}gt


EOF

qsub $subFile

rm ${subFile}