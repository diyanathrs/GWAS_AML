#!bin/bash

module purge
module load R/4.2.1-foss-2022a

# added new code to include NCL6 & 7 in grepl in forestplot_addMAF.R
Rscript --verbose  forestplot_addMAF.R Normal > forestplot_addMAF.Rout

exit

