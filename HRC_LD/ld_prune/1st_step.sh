#!/bin/bash

#get the snp list from titania and use it to extract snps

plink --bfile ../merged_plink/HRC_chr.2 --extract HRCsnp_extract.txt --make-bed --out chr2_meta_filtered




