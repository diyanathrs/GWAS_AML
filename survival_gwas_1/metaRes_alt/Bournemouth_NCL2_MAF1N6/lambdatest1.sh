#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 infoindex"; exit; }

newdir=Hitsmeta
[ -d ${newdir} ] || mkdir ${newdir}

newdir=MH_QQplots
[ -d ${newdir} ] || mkdir ${newdir}


ncl4folder1="CLL_finalQced4imp"

ncl4folder="CLL_finalQced4imp_studyMAF1"

# no. phenos
# nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype') 
declare -a nophenos=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')   
echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${ncl4folder}meta$1.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=16G
#$ -l h='!compute3-11.clusterlan'
#$ -l h='!compute3-4.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 

module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" infoindex="'$1'" '  meta_Lamdatest1.R


EOF

jobname="meta_${nophenos[i]}_lambda1_$1"
#qsub -N ${jobname} ${subFile}

sleep 1

#rm ${subFile}

done





