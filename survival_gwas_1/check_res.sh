#!/bin/bash

prefix=("ELN2noAPL"  "ELN2noSCT"  "noAPL" "Normal"  "noSCTnoAPL")

for index in "${!prefix[@]}";
	do
	echo "Number of res files in ${prefix[$index]}"
	ls ResultSummary_${prefix[$index]}/*.rds | wc -l
	done

#check GWAStabix are done
for index in "${!prefix[@]}";
        do
        echo "Checking GWAS Tabix is done for ${prefix[$index]}"
	if ! find GWAStabix -name "*_${prefix[$index]}assoc_HRC.gz" | grep -q .; then
   	 echo "NO Tabix done for ${prefix[$index]} <<<<<<<<<<"
	fi
	done
	echo "Done!"
