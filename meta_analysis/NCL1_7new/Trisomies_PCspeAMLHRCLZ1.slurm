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
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk
#SBATCH --chdir=/nobackup/proj/jamgaml/dean_AML/HRCimpvData/metaNCL1_5NatComRevTrimSubtypes/NCL1_7new/Hitsmeta

module load R/4.2.1-foss-2022a

pheno=($(cat Trisomies*lzname.txt | awk '{print $1}'))
study=($(cat Trisomies*lzname.txt | awk '{print $2}'))
snps=($(cat Trisomies*lzname.txt | awk '{print $3}'))
chr=($(cat Trisomies*lzname.txt | awk '{print $4}'))
pos=($(cat Trisomies*lzname.txt | awk '{print $5}'))

if [ ${#snps[@]} -eq 0 ]
then

exit

fi

for ((i=0; i<${#snps[@]}; i++))
do

bash /nobackup/proj/jamgaml/dean_AML/LZtest/lzsubmit.sh ${pheno[i]} ${study[i]} ${snps[i]} ${chr[i]} ${pos[i]}


done


