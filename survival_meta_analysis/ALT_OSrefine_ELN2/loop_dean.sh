#!/bin/bash

while IFS= read -r line; do 
	bash SLURM_metaStep0.sh $line; 
	done < mystudy_eln2.lst
