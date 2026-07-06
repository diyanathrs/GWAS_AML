#!/bin/bash

newdir=Hitsmeta
[ -d ${newdir} ] || mkdir ${newdir}

newdir=MH_QQplots
[ -d ${newdir} ] || mkdir ${newdir}


ncl4folder="PCspeAMLHRC"
#ncl4folder2="PCspeAMLHRC_maf"

folderlist="NCLpcspeData.lst"

# no. phenos
# nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)

#declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype' )

# keep NCL6 for only status and Normal
declare -a nophenos=('status' 'Normal')

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
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load PLINK/1.9b_6.21-x86_64-GCC-11.2.0


plink --meta-analysis ../${nophenos[i]}_MAF0.0095*assoc.res + logscale report-all \\
--out ${nophenos[i]}_NCL_${ncl4folder}

# maf version by dean
#plink --meta-analysis ${nophenos[i]}*_MAF_assoc.res + logscale report-all \\
#--out ${nophenos[i]}_NCL_${ncl4folder2}

EOF

command1="sbatch ${subFile}"
jid1=$($command1 | awk ' { print $4 }')

# qq AND MH plots

subFile1="${nophenos[i]}_${ncl4folder}meta_mhqq.slurm"
    
    cat > ${subFile1} << EOF
#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclude=ln03,sb017
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load R/4.2.1-foss-2022a

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryQQMH.R

# different maf by dean
#R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder2}'" phenoprefix="'${nophenos[i]}'" '  meta_summaryQQMH_maf.R



EOF

command2="sbatch --dependency=afterok:${jid1} ${subFile1}"
jid2=$($command2 | awk ' { print $4 }')


# hits & forest plots

subFile2="${nophenos[i]}_${ncl4folder}meta_hits.slurm"
    
    cat > ${subFile2} << EOF
#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=1:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=ln03,sb017
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load R/4.2.1-foss-2022a

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" folderlist="'${folderlist}'"'  meta_summaryHits.R 

# added new code to include NCL6 & 7 in grepl in forestplot_addMAF.R
R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" '  forestplot_addMAF.R

EOF


sbatch --dependency=afterok:${jid2} ${subFile2}

sleep 0.5


done

