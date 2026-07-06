#!/bin/bash
#SBATCH --partition=default_free
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --account=rockhpc_jamgaml
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
##SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk

module load PLINK/1.9b_6.21-x86_64-GCC-11.2.0


plink --meta-analysis ../status_MAF0.02*assoc.res + logscale report-all \
--out status_NCL_PCspeAMLHRC

# maf version by dean
#plink --meta-analysis status*_MAF_assoc.res + logscale report-all \
#--out status_NCL_

