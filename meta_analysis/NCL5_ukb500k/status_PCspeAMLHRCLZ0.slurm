#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --chdir=/nobackup/proj/jamgaml/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL5_ukb500k/Hitsmeta

module load R/3.3.1-intel-2017.03-GCC-6.3
module load HTSlib/1.4.1-intel-2017.03-GCC-6.3

Rscript --vanilla /nobackup/proj/jamgaml/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL5_ukb500k/lz_hitsmetaHRC.R "status" "PCspeAMLHRC"

