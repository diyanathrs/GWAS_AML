#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=ln03,sb017
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

##module load R/4.2.1-foss-2022a
module purge
module load R/3.3.1-intel-2017.03-GCC-6.3
module load PLINK/1.9b_6.21-x86_64-GCC-11.2.0

Rscript --vanilla /nobackup/proj/jamgaml/dean_AML/LZtest/getlocuszoom.R Trisomies PCspeAMLHRC rs145717902 5 32323770



