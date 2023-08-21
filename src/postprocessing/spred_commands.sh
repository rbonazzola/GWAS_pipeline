#!/bin/bash
/home/rodrigo/PhD/software/MetaXcan/software/SPrediXcan.py --model_db_path /media/rodrigo/bd2fea02-a909-49bc-9c6c-e0750941ca71/data/GTEx/GTEx_v8/prediction_models/mashr_eqtl/eqtl/mashr/mashr_Heart_Left_Ventricle.db --covariance /media/rodrigo/bd2fea02-a909-49bc-9c6c-e0750941ca71/data/GTEx/GTEx_v8/prediction_models/mashr_eqtl/eqtl/mashr/mashr_Heart_Left_Ventricle.txt.gz --gwas_folder /media/rodrigo/bd2fea02-a909-49bc-9c6c-e0750941ca71/data/GWAS/output/coma/2020-09-11_02-13-41/test_std_covariates_PC__GBR_unrelated__qc/BGEN --gwas_file_pattern GWAS__z5__harmonized__test_std_covariates_PC__GBR_unrelated__qc.tsv --output_file /media/rodrigo/bd2fea02-a909-49bc-9c6c-e0750941ca71/data/GWAS/output/coma/2020-09-11_02-13-41/test_std_covariates_PC__GBR_unrelated__qc/SPrediXcan/2020-09-11_02-13-41__z5__Heart_Left_Ventricle__GTEx_v8__mashr.csv --snp_column variant_id --effect_allele_column effect_allele --non_effect_allele_column non_effect_allele --beta_column effect_size --pvalue_column pvalue