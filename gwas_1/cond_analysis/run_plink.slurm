#!/bin/bash

#SBATCH --partition=defq,long
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk


input="../association/statusNCL6_chr11_VCFphenocov.sample"
bfile="../hardcall/NCL6_chr11_hardcalls"

module purge
module load R
module load PLINK/2.00a3.7-foss-2022a

R CMD BATCH --vanilla --no-timing '--args input="'${input}'" bfile="'${bfile}'" ' plink_condAnalysis.R
