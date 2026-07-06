#!/bin/bash

rootpath=$(pwd)

# change this path to where your gen and sample files
# datadir="/home/nwl15/WORKING_DATA/imputeTest/TRsurv"
datadir=${rootpath}/impvSubset_eln2

# change file name of your master accordingly
phenofile="AMLSurpheno_uk12GerHun_13Mar2019.txt"

prefix="Finland"

# 
# declare -a study=('Bournemouth' 'Cardiff' 'Hull1' 'Hull2' 'Newcastle1' 'Newcastle2')

declare -a study=('Finland')

# 3 time-dependant variables, xxx-Status, and time
# change if your column names are not the same below
# declare -a cencol=('OS_Dx_Death_LFU_Status' 'OS_Rx_Death_LFU_status' 'TTFT_Dx_Rx_Status')
# declare -a timecol=('OS_Dx_Death_LFU' 'OS_Rx_Death_LFU' 'TTFT_Dx_Rx')

declare -a cencol=('OSstatus' 'RFSstatus' 'OS2status')
declare -a timecol=('OSdays' 'RFSdays' 'OS2days')

# declare -a cencol=('OS_Dx_Death_LFU_Status')
# declare -a timecol=('OS_Dx_Death_LFU')

# for each phenotype, generating the sge script and submitting to the cluster
# for each phenotype, requesting 8 cpus to fit cox models (-pe smp 8), changing this if necesary
# by default, chr1- 22 runs will be submitted as job arrays, if you want to test one chromosome only, change (-t 1-22)
# for example, -t 22-22, will only run chromosome 22 data
# change -M to your email address
for ((j=0; j<${#study[@]}; j++))
##for ((j=0; j<1; j++))
do

    for ((i=0; i<${#cencol[@]}; i++))
    do 
	subFile="${study[j]}_${cencol[i]}_Surv.slurm"
	cat > ${subFile} << EOF
#!/bin/bash
#SBATCH -A jamgaml
#SBATCH --partition=defq,long
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6000
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=d.s.ranasinghe2@newcastle.ac.uk
#SBATCH --array=1-22
#SBATCH --chdir=${rootpath}/association

## SBATCH --workdir=\${TMPDIR}
## SBATCH --exclusive
## SCRATCH=/scratch/\${USER}/\${SLURM_JOB_ID}
## SCRATCH=\${TMPDIR}
## mkdir -p \${SCRATCH} || exit $?
##cp /nobackup/proj/jamgaml/HRCimpvData/CLLbystudy/association/CLL_phenoUpdate_13July2017.txt \${SCRATCH}/ || exit $?
# cp ${rootpath}/association/${phenofile} \${SCRATCH}/ || exit $?
# cp ${datadir}/OxfCLL4imp_chr\${SLURM_ARRAY_TASK_ID}subset*.gen \${SCRATCH}/ || exit $?
# cp ${datadir}/OxfCLL4imp_chr\${SLURM_ARRAY_TASK_ID}subset*.sample \${SCRATCH}/ || exit $?
# cp ${rootpath}/HRbystudyFMSrun_rscript.R \${SCRATCH}/ || exit $?
# cd \${SCRATCH}
## module load R/3.3.1-intel-2017.03-GCC-6.3
# R CMD BATCH --vanilla --no-timing '--args data.dir="'${datadir}'" chrom="'\${SGE_TASK_ID}'" phenofile="'${phenofile}'" cen.col="'${cencol[i]}'" time.col="'${timecol[i]}'" study="'${study[j]}'" ' ${rootpath}/HRbystudyFMSrun_rscript.R ${study[j]}_${cencol[i]}_chr\${SGE_TASK_ID}.Rout
# Rscript --vanilla ${rootpath}/HRbystudyFMSrun_rscript.R "${phenofile}" "${cencol[i]}" "${timecol[i]}" "${datadir}" "${study[j]}"
# Rscript --vanilla ${rootpath}/HRbystudyFMSrun_rscript.R "${phenofile}" "${cencol[i]}" "${timecol[i]}" "\${SCRATCH}" "${study[j]}"
# cp \${SCRATCH}/${study[j]}_${cencol[i]}_chr\${SLURM_ARRAY_TASK_ID}_HRCSurRes.out \${SLURM_SUBMIT_DIR}/association/ || exit $?

module load R

Rscript --vanilla ${rootpath}/HRbystudyFMSrun_rscript.R "${phenofile}" "${cencol[i]}" "${timecol[i]}" "${datadir}" "${study[j]}" "${prefix}"


EOF

	sbatch ${subFile}
	sleep 0.5


    done
    
    sleep 0.5

done


