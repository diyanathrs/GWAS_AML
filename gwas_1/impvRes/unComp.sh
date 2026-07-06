#!/bin/bash

rootpath=`pwd`
## impvpath="${rootpath}/impvRes"
## prefix=CLL_finalQced4imp

#pheno="${rootpath}/NCL5_phenoQced_1821.txt"
#cov="${rootpath}/NCL5_1821.cov"
#pcnum="3"
## echo "$rootpath"
## exit

subFile1="unzipping.sge"
cat > ${subFile1} <<EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=10G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -l h='!compute3-11.clusterlan'
#$ -m a
#$ -t 1-22
#$ -wd ${rootpath}


##unzip -P YtZfgr0Gkx4pUO chr_\${SGE_TASK_ID}.zip
7za x chr_\${SGE_TASK_ID}.zip -p'SUFzGdR0>m5tKx'

EOF

jobname="unzipping"
qsub -N ${jobname} ${subFile1}


