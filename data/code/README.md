# Genetic data pre-processing
## Imputed genotypes
The following instructions allow to filter the set of imputed genotype data (in BGEN format) for a set of SNPs and subjects, and to split these files into regions to ease the parallelization of the subsequent tasks.
1. Generate a preliminary set of SNPs to filter for, so that subsequent tasks are faster, based on MAF and imputation INFO score.
2. Filter for SNPs generated in step 1, for subjects and split the genome into regions.
3. Compute SNP-wise statistics (MAF, missingness) and further filter the BGEN files.

The following paths must be specified. Placeholders for chromosome number can be `chromosome` or `chromosome_adjusted`, according to whether the number is padded to the left with zero for one-digit numbers (`1` or `01`).

| Path                           | Description                                                                                                                        |
|--------------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| `REGIONS_FILE`                 | Path to a tab-separated file containing the columns `chr`, `start`, and `stop`. The chromosomes are in the format `chrN`.          |
| `INPUT_GENOTYPE_FILE_PATTERN`  | Path to the original BGEN files, containing a replaceable `{chromosome}` placeholder, e.g. `data/datasets/imputed/ukb22828_c{chromosome}_b0_v3.bgen`. |
| `INDEX_FILE_PATTERN`           | Path to the BGI index files for the above BGEN files.                                                                              |
| `INCLUDE_SAMPLES`              | Path to a file with one subject ID per line. These subjects are used to calculate the SNP statistics and determine which SNPs to include in the final files. |
| `TEMP_GENOTYPE_FILE_PATTERN`   | File pattern to an intermediate set of BGEN files. These contain the genotype information for all the subjects, for those SNPs passing the QC criteria according to the `INCLUDE_SNP_FILES`. |
| `INCLUDE_SNP_FILES`            | List of SNPs, one file per chromosome, with SNPs passing the QC criteria.                                                          |
| `REDUCED_SNP_FILE_PATTERN`     | Same as above, but per region.                                                                                                     |
| `OUTPUT_GENOTYPE_FILE_PATTERN` | Same as `TEMP_GENOTYPE_FILE_PATTERN`, but also filtering for subjects.                                                             |
| `SNPS_STATS_FILE_PATTERN`      | SNP statistics computed within the `INCLUDE_SAMPLES` file.                                                                         |
| `SNPS_LIST_FILE_PATTERN`       | List of SNPs, one file per region, with SNPs passing the QC criteria.                                                              |
| `FINAL_BGEN_FILE_PATTERN`      | Final BGEN files, one per region, for the subjects in the `INCLUDE_SAMPLES` file, with SNPs from the `SNPS_LIST_FILE_PATTERN` files.|


### Step 1
Generate a file with the SNPs to be kept in a preliminary filtering stage:
```
bash generate_include_snps_list.sh
```

### Step 2
Filter BGENs as per the list of SNPs generated in 1, a list of subjects, and split into small regions:
```
python split_bgen_by_region_qctool.py
```

### Step 3
Compute SNP-wise statistics and further filter the BGENs from step 2:
```
python filter_bgens_for_qc_criteria.py
```
