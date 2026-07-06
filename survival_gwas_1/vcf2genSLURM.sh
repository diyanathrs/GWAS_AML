#!/bin/bash

(($# >= 1)) || { echo -e "\nUsage: bash $0 sampleKeep_list"; exit; }


rootpath=$(pwd)
impvpath="${rootpath}/impvRes"
prefix=Finland
samplelist="$1"

echo "${rootpath}/impvSubset/${samplelist}"

if [ -e "${rootpath}/impvSubset/${samplelist}" ]
then
    echo "file exists, Yes"
else 
    echo "not found"
    exit
fi

subFile1="vcf2gen.slurm"
cat > ${subFile1} <<EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long,bigmem
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=Wei-Yu.Lin@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --workdir=${rootpath}/impvSubset


##module load R/3.3.1-intel-2017.03-GCC-6.3
##R CMD BATCH --vanilla --no-timing '--args root.dir="'${rootpath}'" dataprefix="'${prefix}'" chrom="'\${SLURM_ARRAY_TASK_ID}'" '  ${rootpath}/00_step1.R
##qctool -g ${impvpath}/chr\${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -incl-rsids ${prefix}_chr\${SLURM_ARRAY_TASK_ID}_MarkersPassed.lst -vcf-genotype-field GP -og ${prefix}_chr\${SLURM_ARRAY_TASK_ID}subset.gen -os ${prefix}_chr\${SLURM_ARRAY_TASK_ID}subset.sample

module load BCFtools/1.3-intel-2017.03-GCC-6.3

bcftools view -i 'MAF >= 0.0095  && R2 >= 0.3' --force-samples -S ${samplelist} -Oz -o ${prefix}.chr\${SLURM_ARRAY_TASK_ID}.vcf.gz --threads \${SLURM_NTASKS_PER_NODE} ${impvpath}/chr\${SLURM_ARRAY_TASK_ID}.dose.vcf.gz

qctool_v2.0.1 -g ${prefix}.chr\${SLURM_ARRAY_TASK_ID}.vcf.gz -vcf-genotype-field GP -og ${prefix}_chr\${SLURM_ARRAY_TASK_ID}subset.gen -os ${prefix}_chr\${SLURM_ARRAY_TASK_ID}subsetMOD.sample




EOF

sbatch ${subFile1}



