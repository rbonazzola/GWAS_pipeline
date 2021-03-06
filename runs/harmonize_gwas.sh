#!/bin/bash

# This script performs harmonization with liftover (from hg19 to hg38)
# of GWAS summary statistics produced with the BGENIE tool

# TODO: add summary-gwas-imputation as submodule
# TODO: read paths from another file
GENETICS_FOLDER="${GENETICS_FOLDER}"
EXPERIMENT="2020-09-11_02-13-41"

for NUMBER in `seq 1 7`; do
  python ${GENETICS_FOLDER}/summary-gwas-imputation/src/gwas_parsing.py \
    -gwas_file /home/home01/scrb/src/GWAS_pipeline/output/coma/2020-09-11_02-13-41/BGEN/GWAS__z${NUMBER}__std_covariates__GBR__BGEN__qc.tsv \
    -output_column_map SNP variant_id \
    -output_column_map a_1  effect_allele \
    -output_column_map a_0  non_effect_allele \
    -output_column_map BP position \
    -output_column_map CHR chromosome --chromosome_format \
    -output_column_map af frequency \
    -output_column_map P pvalue  \
    -output_column_map BETA effect_size \
    -output_column_map SE standard_error \
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error  \
    -snp_reference_metadata ${GENETICS_FOLDER}/gwas_ref_data/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
    -liftover ${GENETICS_FOLDER}/gwas_ref_data/data/liftover/hg19ToHg38.over.chain.gz \
    -output /home/home01/scrb/src/GWAS_pipeline/output/coma/${EXPERIMENT}/BGEN/GWAS__z${NUMBER}__std_covariates__GBR__BGEN__qc__harmonized.txt.gz
done
