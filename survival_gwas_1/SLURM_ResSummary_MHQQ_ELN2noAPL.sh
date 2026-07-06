#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 studyName"; exit; }

subset=ELN2noAPL
newdir=ResultSummary_${subset}
[ -d ${newdir} ] || mkdir ${newdir}

##rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE"
rootpath=$(pwd)

#impvpath="${rootpath}/impvRes"
prefix=Finland_${subset}
assoc=association_${subset}

#declare -a nophenos=('OSstatus' 'RFSstatus' )

declare -a nophenos=('OSstatus')

nophenos=("${nophenos[@]/#/$1_}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""

echo ${nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    # tidying up the results based on maf and info
    subFile="${nophenos[i]}_${prefix}_HRCtidyup.slurm"
    
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --export=ALL
##SBATCH --mail-type=FAIL
#SBATCH --array=3-3
#SBATCH --chdir=${rootpath}/${newdir}


module load R/3.4.3-foss-2017b-X11-20171023

R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" assoc.fol="'${assoc}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'"  phenoprefix="'${nophenos[i]}'" '  ${rootpath}/res_summarySur.R

#R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'"  phenoprefix="'${nophenos[i]}'" '  ${rootpath}/res_summarySur.R

EOF

    jobname="${nophenos[i]}_${prefix}tidyup"
    command1="sbatch -J ${jobname} ${subFile}"
    jid1=$($command1 | awk ' { print $4 }')
    

    # QQ and MH plots
    subFile1="${nophenos[i]}_${prefix}_QQMH.slurm"
    cat > ${subFile1} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --export=ALL
##SBATCH --mail-type=FAIL
#SBATCH --chdir=${rootpath}/${newdir}

module load R/3.4.3-foss-2017b-X11-20171023


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${rootpath}/res_summaryQQMHSur.R

EOF


command2="sbatch --dependency=afterok:${jid1} ${subFile1}"
jid2=$($command2 | awk ' { print $4 }')

done

