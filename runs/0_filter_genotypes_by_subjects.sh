#!/bin/bash

cd $(git rev-parse --show-toplevel)

IDS=data/datasets/ids_list/CMR_GBR_fid_iid.txt
DIR=data/transforms/GenomicPCA/genotypes

for CHR in `seq 1 22`; do 
  plink -bfile $DIR/ukb_cal_chr${CHR}_v2_GBR_indiv --keep $IDS --make-bed --out $DIR/only_cmr/ukb_cal_chr${CHR}_v2_CMR_GBR_indiv
done

