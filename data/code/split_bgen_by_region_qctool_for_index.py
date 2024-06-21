
# bgenix \
# -g data/datasets/imputed/ukb22828_c4_b0_v3.bgen \
# -i data/datasets/imputed/ukb_imp_chr4_v3.bgen.bgi \
# -incl-range 04:${REGION} > data/transforms/genotypes_by_region/ukb_chr04_${REGION}_all.bgen
# 
# qctool \
# -g data/transforms/genotypes_by_region/ukb_${CHR}_${REGION}_all.bgen \
# -og data/transforms/genotypes_by_region/ukb_${CHR}_${REGION}.bgen \
# -s data/datasets/imputed/ukb22828_c1_b0_v3_s487202.sample \
# -incl-rsids data/transforms/snps_files/by_region/ukb_chr4_${REGION}__maf_gt_0.005_info_gt_0.3.txt \ 
# -incl-samples data/ids_list/all_cmr_ids.txt
# 
# rm data/transforms/genotypes_by_region/ukb_${CHR}_${REGION}_all.bgen

import os, sys
import pandas as pd
import argparse
import shlex
from subprocess import call, check_output
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip())

from PATHS import *
from split_bgen_by_region_qctool import *

if __name__ == "__main__":

    df = get_regions_df()

    parser = argparse.ArgumentParser()
    parser.add_argument("--region_index", type=int)
    args = parser.parse_args()

    r = args.region_index
    
    region = df.iloc[r]
    region = get_region_data(region)

    input_genotype_file = INPUT_GENOTYPE_FILE_PATTERN.format(chromosome=region['chromosome'])
    index_file = INDEX_FILE_PATTERN.format(chromosome=region['chromosome'])
    incl_rsid_file = INCLUDE_SNP_FILES[0].format(chromosome=region['chromosome'])
    full_df = pd.read_csv(incl_rsid_file, sep=" ", header=None)    

    reduced_incl_snp_file = generate_incl_snps_file(full_df, region)
    incl_range = "{chromosome_adjusted}:{start_pos}-{end_pos}".format(**region)
    temp_genotype_file = TEMP_GENOTYPE_FILE_PATTERN.format(**region)
    output_genotype_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region)

    bgenix_command = build_bgenix_command(
        input_genotype_file, 
        index_file, 
        temp_genotype_file, 
        incl_range
    )

    qctool_command = build_qctool_command(
        temp_genotype_file, 
        output_genotype_file, 
        incl_rsid_files=[reduced_incl_snp_file], 
        incl_range=None, 
        incl_samples=INCLUDE_SAMPLES, 
        samples=SAMPLES_FILE
    )
    
    print(bgenix_command)
    print(qctool_command)