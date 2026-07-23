#!/bin/bash

#SBATCH --partition=defq,long
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --export=ALL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk
#SBATCH --array=2


# for complex do cojo for chr2 only

./../gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --bfile ../../HRC_LD/merged_bgen_new/HRC_chr.${SLURM_ARRAY_TASK_ID} --chr ${SLURM_ARRAY_TASK_ID} --cojo-file Complex_cond_in.ma --maf 0.02 --cojo-slct --cojo-p 5e-5 --out cojo_slct_complex_chr${SLURM_ARRAY_TASK_ID}

#./../gcta-1.94.4-linux-kernel-3-x86_64/gcta64 --bfile ../../HRC_LD/merged_plink/HRC_chr.11  --maf 0.02 --chr 11  --cojo-file status_cond_in.ma --cojo-cond cond.snplist --out cond_rs4930561
