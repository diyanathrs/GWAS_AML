#!/bin/bash

ncl4folder="PCspeAMLHRC"

# no. phenos
# nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
# declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype' 'abnormality' 'anyTriMon' 'transNOt1517')

declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies')
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
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk


module load R/3.3.1-intel-2017.03-GCC-6.3


R CMD BATCH --vanilla --no-timing '--args dataprefix="'${ncl4folder}'" phenoprefix="'${nophenos[i]}'" ' metaRes2tabix.R

module load HTSlib/1.4.1-intel-2017.03-GCC-6.3

head -n 1 ${nophenos[i]}_${ncl4folder}_assocLZ.res | awk '{for (i=1;i<=NF;i++) print \$i}' | paste -s | awk -F"\t" '{print "#"\$0}' > ${nophenos[i]}_${ncl4folder}.header
tail -n+2 ${nophenos[i]}_${ncl4folder}_assocLZ.res | sort -T ./ -k1,1n -k2,2n | bgzip >  ${nophenos[i]}_${ncl4folder}TMP.gz
tabix -r ${nophenos[i]}_${ncl4folder}.header ${nophenos[i]}_${ncl4folder}TMP.gz > ${nophenos[i]}_${ncl4folder}assoc_HRC.gz
tabix -s 1 -b 2 -e 2 ${nophenos[i]}_${ncl4folder}assoc_HRC.gz

rm ${nophenos[i]}_${ncl4folder}TMP.gz ${nophenos[i]}_${ncl4folder}.header ${nophenos[i]}_${ncl4folder}_assocLZ.res

EOF

sbatch ${subFile}

done

