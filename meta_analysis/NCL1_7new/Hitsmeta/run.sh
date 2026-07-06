#!bin/bash

module purge
module load R/4.2.1-foss-2022a

declare -a nophenos=('status' 'Normal' 'CBF' 't15_17'  'Complex' 'del5_7' 'Trisomies' 'Any.Monosomy' 'anyTriMon' 'abnormality' 'monosomal.karyotype' 'Trans')

#declare -a nophenos=('status' 'Normal')
#Any.Monosomy Trans
#monosomal.karyotype

echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}


for ((i=0; i<${#nophenos[@]}; i++))
do
echo ${nophenos[i]}
# added new code to include NCL6 & 7 in grepl in forestplot_addMAF.R
Rscript --verbose  forestplot_addMAF.R ${nophenos[i]} > forestplot_addMAF.Rout

done
exit

