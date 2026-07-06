#!/bin/bash

folderlist=NCLpcspeData.lst
nodata=(`cat ${folderlist} | awk '{print $1}'`)


# no. phenos
##nophenos=(`cat ${listdir}/${prefix}_chr1pheno.lst | awk '{print $1}' `)
declare -a nophenos=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype'  'anyTriMon' 'abnormality')


echo ""
echo "no. phenos are ${#nophenos[@]}"
echo ""
echo ${nophenos[@]}
echo ${#nophenos[@]}


for proj in ${nodata[@]};
do
echo "For $proj.."
for pheno in ${nophenos[@]}
do
count=(`ls $proj/$pheno*.rds | wc -l`)
echo "Rds files in $pheno = $count"
if (($count < 22))
then
echo "============<<<<<<<<"
fi
done
done



  
                                            
                                                                     

