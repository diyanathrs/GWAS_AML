#!/bin/bash
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk

plink --meta-analysis Complex*_assoc.res + logscale report-all \
--out Complex_NCL_PCspeAMLHRC


