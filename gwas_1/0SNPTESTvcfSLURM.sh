#!/bin/bash

(($# == 4)) || { echo -e "\nUsage: bash $0 prefix pheno cov pcnum"; exit; }

##rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/NCL5_ukb500k"
scriptpath=/nobackup/proj/jamgaml/HRCimpvData
rootpath=$(pwd)

impvpath="${rootpath}/impvRes"

prefix="$1"
pheno="${rootpath}/$2"
cov="${rootpath}/$3"
pcnum="$4"

echo $prefix
echo $pheno
echo $cov
echo $pcnum


# prefix="NCL5_finalQced4imp"
# pheno="${rootpath}/NCL5_ukb500kQcedPheno.txt"
# cov="${rootpath}/NCL5_sampleCleanPC.cov"
# pcnum="3"

subFile1="vcfpheno.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --workdir=${rootpath}/association


module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args impv.path="'${impvpath}'" prefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'"  pheno="'${pheno}'" cov="'${cov}'" pcnum="'${pcnum}'" '  ${scriptpath}/phenoPrepVCF.R


EOF

command1="sbatch ${subFile1}"
jid1=$($command1 | awk ' { print $4 }')

# jobname="vcfpheno"
# qsub -N ${jobname} ${subFile1}


# no. phenos
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')
#declare -a nophenos=('Normal') 
echo ""
echo "no. phenos for chr ${chrom} are ${#nophenos[@]}"
echo ""

for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${prefix}_HRCvcf_SNPTEST.slurm"
    
    cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --workdir=${rootpath}/association

snptest_v2.5.2 \\
-data ${impvpath}/chr\${SLURM_ARRAY_TASK_ID}.dose.vcf.gz ${prefix}_chr\${SLURM_ARRAY_TASK_ID}_VCFphenocov.sample \\
-genotype_field GP \\
-frequentist 1 -method expected \\
-cov_all -hwe \\
-pheno ${nophenos[i]} \\
-o ${nophenos[i]}_${prefix}_chr\${SLURM_ARRAY_TASK_ID}_HRCRes.out.gz


EOF


command2="sbatch --dependency=afterok:${jid1} ${subFile}"
jid2=$($command2 | awk ' { print $4 }')

# summary results
subFile2="${nophenos[i]}_${prefix}_HRCtidyup.slurm"
cat > ${subFile2} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --workdir=${rootpath}/ResultSummary

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'"  phenoprefix="'${nophenos[i]}'" '  ${scriptpath}/res_summary.R


EOF


command3="sbatch --dependency=afterok:${jid2} ${subFile2}"
jid3=$($command3 | awk ' { print $4 }')


# QQ MH 
subFile3="${nophenos[i]}_${prefix}_QQMH.slurm"
cat > ${subFile3} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --workdir=${rootpath}/ResultSummary

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${scriptpath}/res_summaryQQMH.R


EOF


command4="sbatch --dependency=afterok:${jid3} ${subFile3}"
jid4=$($command4 | awk ' { print $4 }')


sleep 1


done


