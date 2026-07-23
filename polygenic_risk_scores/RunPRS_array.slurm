#!/bin/bash
#SBATCH --account=rockhpc_jamgaml
#SBATCH --partition=default_free
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --export=ALL
#SBATCH --array=0-5
#SBATCH --output=logs/prs_%A_%a.out
#SBATCH --error=logs/prs_%A_%a.err

module load plink/2.0.0

# Score file from COJO results
SCORE=PanAML_score.txt

# Prefixes of the merged PLINK datasets
DATASETS=(
"../NCL1_2NatComRevTrim/hardcall/NCL12_allchr"
"../NCL3NatComRevTrim/hardcall/NCL3_allchr"
"../NCL4newNatComRevTrim/hardcall/NCL4_allchr"
"../NCL5_ukb500kNatComRevTrim/hardcall/NCL5_allchr"
"../NCL6/hardcall/NCL6_allchr"
"../NCL7_mich/hardcall/NCL7_allchr"
)

BFILE=${DATASETS[$SLURM_ARRAY_TASK_ID]}
STUDY=$(basename ${BFILE})

echo "Running PRS for ${STUDY}"

plink2 \
    --bfile ${BFILE} \
    --score ${SCORE} 1 2 3 header cols=+scoresums \
    --out ${STUDY}

echo "Finished ${STUDY}"
