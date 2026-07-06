#!/bin/bash

newdir=Hitsmeta
[ -d ${newdir} ] || mkdir ${newdir}

newdir=MH_QQplots
[ -d ${newdir} ] || mkdir ${newdir}


ncl4folder="gt_NCL1_5_ukb500kHRC"

# no. phenos
# nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')
# declare -a nophenos=('status' 'Normal')
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
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 

plink --meta-analysis ${nophenos[i]}*_assoc.res + logscale report-all \\
--out ${nophenos[i]}_NCL_${ncl4folder}


EOF

jobname="meta_${nophenos[i]}_NCL_${ncl4folder}"
qsub -N ${jobname} ${subFile}

# qq AND MH plots

subFile1="${nophenos[i]}_${ncl4folder}meta_mhqq.sge"
    
    cat > ${subFile1} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan' 
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 

##module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryQQMH.R


EOF

qsub -hold_jid ${jobname} ${subFile1}

# hits & forest plots

subFile2="${nophenos[i]}_${ncl4folder}meta_hits.sge"
    
    cat > ${subFile2} << EOF
#!/bin/bash
#$ -cwd -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'  
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 

##module load apps/R/3.2.3


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryHits.R 

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  forestplot_addMAF.R

EOF


qsub -hold_jid ${jobname} ${subFile2}

sleep 0.5


done

