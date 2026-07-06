#!/bin/bash

rootpath=$(pwd)

impvpath="${rootpath}/impvRes"
script_path="/nobackup/proj/jamgaml/dean_AML/HRCimpvData"

prefix="NCL6"
pheno="${rootpath}/NCL6_phenotypeR1_12Feb2019.txt"
## cov="${rootpath}/ncl1_2_3836.cov"
## pcnum="3"
cov="NCL6_sampleCleanPC.cov"


# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype' 'abnormality' 'anyTriMon' 'transNOt1517' 'panAMLnoAPL')
## declare -a nophenos=('status' 'Normal')
## declare -a nophenos=('panAMLnoAPL')

declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies')


for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile1="vcfpheno.slurm"
    cat > ${subFile1} <<EOF
#!/bin/bash                                                                                             
#SBATCH -A jamgaml                                                                                      
#SBATCH --partition=defq,long,bigmem,short                                                                
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                                                                                      
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G                                                                                
#SBATCH --exclude=ln03,sb017
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#####SBATCH --array=1-22
#SBATCH --chdir=${rootpath}/association


module load R/3.3.1-intel-2017.03-GCC-6.3

###R CMD BATCH --vanilla --no-timing '--args impv.path="'${impvpath}'" prefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'"  pheno="'${pheno}'" cov="'${cov}'" pcnum="'${pcnum}'" '  ${script_path}/phenoPrepVCF.R

R CMD BATCH --vanilla --no-timing '--args impv.path="'${impvpath}'" prefix="'${prefix}'" pheno="'${pheno}'" cov="'${cov}'"  phenoname="'${nophenos[i]}'"' ${script_path}/phenoPrepVCFnew.R


EOF

#command1="sbatch ${subFile1}"
#jid1=$($command1 | awk ' { print $4 }')

echo ${jid1}

#jid2+=(${jid1})


done

#depid1=$(IFS=:; echo "${jid2[*]}")

#echo ${depid1}



## jobname="vcfpheno"
## qsub -N ${jobname} ${subFile1}
## command1="sbatch -J ${jobname} ${subFile1}"
## jid1=$($command1 | awk ' { print $4 }')
# no. phenos
# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype' 'abnormality' 'anyTriMon' 'transNOt1517' 'panAMLnoAPL')
## declare -a nophenos=('status' 'Normal')
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies')
#declare -a nophenos=('panAMLnoAPL')
#echo ""
#echo "no. phenos for chr ${chrom} are ${#nophenos[@]}"
#echo ""

for ((i=0; i<${#nophenos[@]}; i++))
do
    subFile="${nophenos[i]}_${prefix}_HRCvcf_SNPTEST.slurm"
    
    cat > ${subFile} << EOF
#!/bin/bash                                                                                               
#SBATCH -A jamgaml                                                                                        
#SBATCH --partition=defq,long                                                                      
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1                                                                                      
#SBATCH --cpus-per-task=1                                               
#SBATCH --mem-per-cpu=16G
#SBATCH --exclude=sb017                                                                         
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --chdir=${rootpath}/association

#snptest_v2.5.2 \\
#-data ${impvpath}/chr\${SLURM_ARRAY_TASK_ID}.dose.vcf.gz ${nophenos[i]}${prefix}_chr\${SLURM_ARRAY_TASK_ID}_VCFphenocov.sample \\
#-genotype_field GP \\
#-frequentist 1 -method expected \\
#-cov_all -hwe \\
#-pheno ${nophenos[i]} \\
#-o ${nophenos[i]}_${prefix}_chr\${SLURM_ARRAY_TASK_ID}_HRCRes.out.gz

# change to ResultSummary folder

cd ../ResultSummary/

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'" phenoprefix="'${nophenos[i]}'" '  ${script_path}/res_summary.R


EOF

# command2="sbatch --partition=long --array=1-2 --dependency=afterok:${jid1} ${subFile}"
# jid2=$($command2 | awk ' { print $4 }')

# echo ${jid2[@]}

# command22="sbatch --partition=defq,long,bigmem --array=3-22 --dependency=afterok:${jid1} ${subFile}"
# jid22=$($command22 | awk ' { print $4 }')

# echo ${jid22[@]}

# command3="sbatch --dependency=afterok:${jid2} ${subFile2}"
# jid3=$($command2 | awk ' { print $4 }')


#command2="sbatch --dependency=afterok:${depid1} ${subFile}"
#jid22=$($command2 | awk ' { print $4 }')

echo "job ids for snptest run are ${jid2[@]}"


##jid4+=(${jid2[@]} ${jid22[@]})
#jid4+=(${jid22[@]})

#echo ${jid4[@]}

#depids=$(IFS=:; echo "${jid4[*]}")

#echo $depids


subFile3="${nophenos[i]}_${prefix}_qqmh.slurm"
cat > ${subFile3} <<EOF
#!/bin/bash
#SBATCH -A jamgaml                                                                                      
#SBATCH --partition=defq,long,bigmem                                                                    
#SBATCH --time=48:00:00                                                                                 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --exclude=ln03,sb017
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --chdir=${rootpath}/ResultSummary

module load R/3.3.1-intel-2017.03-GCC-6.3

R CMD BATCH --vanilla --no-timing '--args dataprefix="'${prefix}'" phenoprefix="'${nophenos[i]}'" '  ${script_path}/res_summaryQQMH.R

EOF

#command4="sbatch --dependency=afterok:${depids} ${subFile3}"
#jid4=$($command4 | awk ' { print $4 }')



done


## qsub -hold_jid ${jobname} ${subFile}
## sleep 1
## rm ${subFile}

# subFile2="${nophenos[i]}_${prefix}_ResSummary.slurm"
# cat > ${subFile2} <<EOF
# #!/bin/bash                                                                                               
# #SBATCH -A jamgaml                                                                                        
# #SBATCH --partition=defq,long,bigmem                                                                      
# #SBATCH --time=48:00:00                                                                                   
# #SBATCH --nodes=1                                                                                         
# #SBATCH --ntasks=1                                                                                        
# #SBATCH --cpus-per-task=1                                                                                 
# #SBATCH --mem-per-cpu=16G                                                                                  
# #SBATCH --export=ALL                                                                                      
# #SBATCH --mail-type=FAIL                                                                                  
# #SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk                                                            
# #SBATCH --array=1-22                                                                                      
# #SBATCH --chdir=${rootpath}/ResultSummary


# module load R/3.3.1-intel-2017.03-GCC-6.3

# R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'" phenoprefix="'${nophenos[i]}'" '  ${script_path}/res_summary.R


# EOF
