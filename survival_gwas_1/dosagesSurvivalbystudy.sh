#!/bin/bash

rootpath=/home/nwl15/WORKING_DATA/HRCimpvData/CLL_OEE

# change this path to where your gen and sample files
# datadir="/home/nwl15/WORKING_DATA/imputeTest/TRsurv"
datadir=${rootpath}/impvSubset

# change file name of your master accordingly
phenofile=CLL_Hull3pheno.txt


# 
## declare -a study=('Bournemouth' 'Cardiff' 'Hull1' 'Hull2' 'Newcastle1' 'Newcastle2')
declare -a study=('Hull3')

# 3 time-dependant variables, xxx-Status, and time
# change if your column names are not the same below
# declare -a cencol=('OS_Dx_Death_LFU_Status' 'OS_Rx_Death_LFU_status' 'TTFT_Dx_Rx_Status')
# declare -a timecol=('OS_Dx_Death_LFU' 'OS_Rx_Death_LFU' 'TTFT_Dx_Rx')

declare -a cencol=('OS_Dx_Death_LFU_Status' 'TTFT_DX_RX_Status' 'OS_RX_to_Death_LFU_status')
declare -a timecol=('OS_Dx_Death_LFU' 'TTFT_DX_RX' 'OS_RX_to_Death_LFU')



# for each phenotype, generating the sge script and submitting to the cluster
# for each phenotype, requesting 8 cpus to fit cox models (-pe smp 8), changing this if necesary
# by default, chr1- 22 runs will be submitted as job arrays, if you want to test one chromosome only, change (-t 1-22)
# for example, -t 22-22, will only run chromosome 22 data
# change -M to your email address
for ((j=0; j<${#study[@]}; j++))
#for ((j=0; j<2; j++))
do

    for ((i=0; i<${#cencol[@]}; i++))
    do 
	subFile="${study[j]}_${cencol[i]}_Surv.sge"
	cat > ${subFile} << EOF
#!/bin/bash
#$ -V
#$ -l h_vmem=6G
#$ -l h='!compute3-11.clusterlan'
#$ -R y
#$ -j yes
#$ -l h_rt=168:00:00   
#$ -m a
#$ -M Wei-Yu.Lin@newcastle.ac.uk 
#$ -t 1-22
#$ -pe smp 8
#$ -wd ${rootpath}/association

R CMD BATCH --vanilla --no-timing '--args data.dir="'${datadir}'" chrom="'\${SGE_TASK_ID}'" phenofile="'${phenofile}'" cen.col="'${cencol[i]}'" time.col="'${timecol[i]}'" study="'${study[j]}'" ' ${rootpath}/HRbystudyFMSrun_rscript.R ${study[j]}_${cencol[i]}_chr\${SGE_TASK_ID}.Rout


EOF


	qsub -N ${study[j]}_${cencol[i]} ${subFile}
	sleep 1

    done
    
    sleep 1

done

