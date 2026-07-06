#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load PLINK/1.9b_6.21-x86_64-GCC-11.2.0


plink --meta-analysis ../assoc_res_maf2_50cases/t15_17_MAF0.02*assoc.res + logscale report-all \
--out t15_17_NCL_PCspeAMLHRC


