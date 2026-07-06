#!/bin/bash

(($# >= 3)) || { echo -e "\nUsage: bash $0 studyName resdir prefix"; exit; }


### declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')

declare -a nophenos=('OS2status')

nophenos=("${nophenos[@]/#/$1_}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}

rootpath=$(pwd)
##prefix="OxfCLL4imp"
##resdir="../ResultSummary"
resdir="$2"
prefix="$3"


echo ${resdir}
echo ${prefix}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_metaformating.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --workdir=${rootpath}

module load R/3.4.3-foss-2017b-X11-20171023

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" res.dir="'${resdir}'" ' metaNCL_gwasALT.R


EOF

sbatch ${subFile}
sleep 1

rm ${subFile}

done




# qsub ${subFile}
