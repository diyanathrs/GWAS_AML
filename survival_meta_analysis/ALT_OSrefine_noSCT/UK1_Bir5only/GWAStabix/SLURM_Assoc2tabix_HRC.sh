#!/bin/bash

prefix="AMLsur_noSCTnoAPL"

condit=""

# no. phenos

##declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')
## declare -a nophenos=('OS_Dx_Death_LFU_Status')

declare -a nophenos=('OSstatus')

nophenos=("${nophenos[@]/#/$1}")
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_resformat.slurm"
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk

module purge
module load R/3.4.3-foss-2017b-X11-20171023

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" ' plinkmeta_summary.R


EOF

command1="sbatch ${subFile}"
jid1=$($command1 | awk ' { print $4 }')

#jobname="${nophenos[i]}_resOUT"
#qsub -N ${jobname} ${subFile}

subFile1="${nophenos[i]}_assoc2tabix.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
##SBATCH --mail-type=FAIL
##SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk

module purge
module load HTSlib/1.9-foss-2017b

head -n 1 ${nophenos[i]}_${prefix}_assocLZ.res | awk -F"\t" '{print "#"\$0}' > ${nophenos[i]}_${prefix}.header
tail -n+2 ${nophenos[i]}_${prefix}_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  ${nophenos[i]}_${prefix}TMP.gz
tabix -r ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}TMP.gz > ${nophenos[i]}_${prefix}${condit}assoc_HRC.gz
tabix -s 1 -b 2 -e 2 ${nophenos[i]}_${prefix}${condit}assoc_HRC.gz

rm ${nophenos[i]}_${prefix}TMP.gz ${nophenos[i]}_${prefix}.header ${nophenos[i]}_${prefix}_assocLZ.res
EOF


command2="sbatch --dependency=afterok:${jid1} ${subFile1}"
jid2=$($command2 | awk ' { print $4 }')

rm ${subFile} ${subFile1}


done
