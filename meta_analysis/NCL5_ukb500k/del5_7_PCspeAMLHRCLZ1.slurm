#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,short,long,bigmem
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --exclude=ln03,sb017
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --chdir=/nobackup/proj/jamgaml/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL5_ukb500k/Hitsmeta

module load R/3.3.1-intel-2017.03-GCC-6.3

pheno=($(cat del5_7*lzname.txt | awk '{print $1}'))
study=($(cat del5_7*lzname.txt | awk '{print $2}'))
snps=($(cat del5_7*lzname.txt | awk '{print $3}'))
chr=($(cat del5_7*lzname.txt | awk '{print $4}'))
pos=($(cat del5_7*lzname.txt | awk '{print $5}'))

if [ ${#snps[@]} -eq 0 ]
then

exit

fi

for ((i=0; i<${#snps[@]}; i++))
do

bash /nobackup/proj/jamgaml/LZtest/lzsubmit.sh ${pheno[i]} ${study[i]} ${snps[i]} ${chr[i]} ${pos[i]}


done


