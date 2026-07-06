#!/bin/bash

prefix="NCL6"
declare -a phenolist=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies' 'Any.Monosomy' 'monosomal.karyotype' 'abnormality' 'anyTriMon')
#declare -a phenolist=('status' 'Normal' 'CBF' 'Trans' 'Complex' 't15_17' 'del5_7' 'Trisomies')
#declare -a phenolist=('status')

for ((i=0; i<${#phenolist[@]}; i++))
do
    # 
## chunks=(`find ./ -type f -name "${phenolist[i]}_${prefix}_HRCvcf_SNPTEST.sge.o*"`)

    chunks=(`find ./ -type f -name "slurm*.out" | xargs grep -l "${phenolist[i]} ~ (baseline)"`)
    awk '/^Phenotype\ summary/ {printline = 1}/^Data\ Summaries/ {printline=0} printline' ${chunks[0]} | grep -v ^$ | cat -n  > ${phenolist[i]}_${prefix}_model.summary
    for ((j=1; j<${#chunks[@]}; j++))
    do
	awk '/^Phenotype\ summary/ {printline = 1}/^Data\ Summaries/ {printline=0} printline' ${chunks[j]} | grep -v ^$ | cat -n >> ${phenolist[i]}_${prefix}_model.summary
    done
    cat ${phenolist[i]}_${prefix}_model.summary | sort | uniq > ${phenolist[i]}_${prefix}_model.summary1


done




rm *.summary
more *.summary1 > ${prefix}_modelHRC.summary
rm *.summary1
