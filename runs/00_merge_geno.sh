#!/bin/bash

BFILE=ukb_cal_chr1_v2_CMR_GBR_indiv
FILES_TO_MERGE="files_to_merge.txt"
OUT=ukb_cal_all_chrs_v2_CMR_GBR_indiv

plink --bfile $BFILE --merge-list $FILES_TO_MERGE --make-bed --out $OUT
