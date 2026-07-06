#!/bin/bash

rootpath="/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE"
impvpath="${rootpath}/impvRes"
prefix=CLL_OEEfinalQced4imp

#pheno="${rootpath}/NCL5_phenoQced_1821.txt"
#cov="${rootpath}/NCL5_1821.cov"
#pcnum="3"

subFile1="vcf2gen.sge"
cat > ${subFile1} <<EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=10G
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -t 1-22
#$ -wd ${rootpath}/impvSubset


R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SGE_TASK_ID}'" '  ${rootpath}/00_step1.R


qctool -g ${impvpath}/chr\${SGE_TASK_ID}.dose.vcf.gz -incl-rsids ${prefix}_chr\${SGE_TASK_ID}_MarkersPassed.lst -vcf-genotype-field GP -og ${prefix}_chr\${SGE_TASK_ID}subset.gen -os ${prefix}_chr\${SGE_TASK_ID}subset.sample


EOF

jobname="vcf2gen"
qsub -N ${jobname} ${subFile1}


# # no. phenos
# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype')
# #declare -a nophenos=('status' 'Normal') 
# echo ""
# echo "no. phenos for chr ${chrom} are ${#nophenos[@]}"
# echo ""

# for ((i=0; i<${#nophenos[@]}; i++))
# do
#     subFile="${nophenos[i]}_${prefix}_HRCvcf_SNPTEST.sge"
    
#     cat > ${subFile} << EOF
# #!/bin/bash
# #$ -V
# #$ -l h_vmem=12G
# #$ -R y
# #$ -j yes
# #$ -l h_rt=168:00:00   
# #$ -m a
# #$ -t 1-22
# #$ -wd ${rootpath}/association

# snptest_v2.5.2 \\
# -data ${impvpath}/chr\${SGE_TASK_ID}.dose.vcf.gz ${prefix}_chr\${SGE_TASK_ID}_VCFphenocov.sample \\
# -genotype_field GP \\
# -frequentist 1 -method expected \\
# -cov_all -hwe \\
# -pheno ${nophenos[i]} \\
# -o ${nophenos[i]}_${prefix}_chr\${SGE_TASK_ID}_HRCRes.out.gz


# EOF

# qsub -hold_jid ${jobname} ${subFile}
# sleep 1
# #rm ${subFile}

# done


