#!/bin/bash

newdir=Hitsmeta
[ -d ${newdir} ] || mkdir ${newdir}

newdir=MH_QQplots
[ -d ${newdir} ] || mkdir ${newdir}


rootpath=$(pwd)

ncl4folder="AMLsur"

# no. phenos
declare -a nophenos=('RFSstatus')

echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${ncl4folder}meta.slurm"
    
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=FAIL
#SBATCH --workdir=${rootpath}

##plink --meta-analysis *${nophenos[i]}_*_assoc.res + logscale report-all \\
##--out ${nophenos[i]}_NCL_${ncl4folder}

##keep N==1

plinkb64 --meta-analysis *${nophenos[i]}_*_assoc.res + logscale report-all \\
--out ${nophenos[i]}_NCL_${ncl4folder}


EOF

#jobname="meta_${nophenos[i]}_NCL_${ncl4folder}"
#qsub -N ${jobname} ${subFile}

command1="sbatch ${subFile}"
jid1=$($command1 | awk ' { print $4 }')

# qq AND MH plots

subFile1="${nophenos[i]}_${ncl4folder}meta_mhqq.slurm"
    
    cat > ${subFile1} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --workdir=${rootpath}

module load R/3.4.3-foss-2017b-X11-20171023

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryQQMH.R


EOF


command2="sbatch --dependency=afterok:${jid1} ${subFile1}"
jid2=$($command2 | awk ' { print $4 }')


# hits & forest plots

subFile2="${nophenos[i]}_${ncl4folder}meta_hits.slurm"
    
    cat > ${subFile2} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --workdir=${rootpath}

module load R/3.4.3-foss-2017b-X11-20171023

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryHits.R 

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  forestplot_addMAF.R

EOF


command3="sbatch --dependency=afterok:${jid1} ${subFile2}"
jid3=$($command3 | awk ' { print $4 }')

sleep 1

done






