#!/usr/bin/env bash

(($# >= 2)) || { echo -e "\nUsage: bash $0 meta summary statistics file"; exit; }

# for each meta file


Rscript prep_meta.R $1

python ../munge_polyfun_sumstats.py \
  --sumstats sumstat.meta \
  --n 12938 \
  --out $2_munged.parquet 


mkdir -p output

python ../polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix output/$2 \
    --sumstats $2_munged.parquet \
    --ref-ld-chr baselineLF2.2.UKB/baselineLF2.2.UKB. \
    --w-ld-chr baselineLF2.2.UKB/weights.UKB. \
    --allow-missing



