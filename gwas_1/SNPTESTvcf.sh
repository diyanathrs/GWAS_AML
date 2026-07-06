#!/bin/bash
rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/NCL5_ukb500k"
impvpath="${rootpath}/impvRes"
prefix="NCL5_finalQced4imp"

pheno="${rootpath}/NCL5_ukb500kQcedPheno.txt"
cov="${rootpath}/NCL5_sampleCleanPC.cov"
pcnum="3"

subFile1="vcfpheno.sge"
cat > ${subFile1} <<EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 
#$ -t 1-22
#$ -wd ${rootpath}/association


##module load apps/R/3.2.3

R CMD BATCH --vanilla --no-timing '--args impv.path="'${impvpath}'" prefix="'${prefix}'" chrom="'\${SGE_TASK_ID}'"  pheno="'${pheno}'" cov="'${cov}'" pcnum="'${pcnum}'" '  /home/nwl15/WORKING_DATA/HRCimpvData/phenoPrepVCF.R

EOF

jobname="vcfpheno"
qsub -N ${jobname} ${subFile1}


# no. phenos
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')
#declare -a nophenos=('status' 'Normal') 
echo ""
echo "no. phenos for chr ${chrom} are ${#nophenos[@]}"
echo ""

for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${prefix}_HRCvcf_SNPTEST.sge"
    
    cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=12G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00
#$ -l h='!compute3-11.clusterlan'
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 
#$ -t 1-22
#$ -wd ${rootpath}/association

snptest_v2.5.2 \\
-data ${impvpath}/chr\${SGE_TASK_ID}.dose.vcf.gz ${prefix}_chr\${SGE_TASK_ID}_VCFphenocov.sample \\
-genotype_field GP \\
-frequentist 1 -method expected \\
-cov_all -hwe \\
-pheno ${nophenos[i]} \\
-o ${nophenos[i]}_${prefix}_chr\${SGE_TASK_ID}_HRCRes.out.gz


EOF

qsub -hold_jid ${jobname} ${subFile}
sleep 1
#rm ${subFile}

done


