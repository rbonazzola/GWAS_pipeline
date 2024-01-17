# GWAS pipeline

This repository contains code to execute GWAS on UK Biobank data, using the Plink and BGENIE tools, and perform downstream analysis on the results.

The folder `src/` contains scripts to perform pre-processing of the data, GWAS execution, and post-processing of the GWAS summary statistics ((like Manhattan plots or Q-Q plot).
The folder `download_data/` contains scripts to download data from the UK Biobank and filter the genotype files.

## Requirements
### Software environment
A Dockerfile is provided for running most of the code in this repository. You can also download the [corresponding Docker image](https://hub.docker.com/r/rbonazzola/gwas_pipeline) from DockerHub.
Alternatively, you can use the Dockerfile as a guide to build your own environment without Docker.

The image is based on Ubuntu 22.04 and contains the following tools:
- R 4.3.1
- Python 3.11
- qctool 2.2.0
- bgenie 1.3
- plink 1.9
- GreedyRelated (to remove related subjects)

If your are on a platform on which Singularity is supported (but Docker isn't), you should be able to convert the Docker image into a Singularity SIF file directly from DockerHub, by using the following command: 

```bash
singularity build <SING_IMAGE_NAME>.sif docker-daemon://<DOCKER_IMAGE_NAME>:<TAG>
```

For computing genetic PCs, a different Docker image is used, namely the one provided by the author of `flashpca`: see https://github.com/gabraham/flashpca/blob/master/docker.md.

For LocusZoom plots, another Docker image will be provided soon.

## Overview
The pipeline consists of scripts for:

1) **Fetching the data**: download genotype and covariate data from the UK Biobank.
2) **Data pre-processing**: 
	- 2a. Genetic data: filtering out subjects and genetic variants and, optionally, splitting the genome into small regions to ease with subsequent parallel processing.
	- 2b. Generate genetic PCs.
	- 2c. Phenotypic data: adjusting the phenotypes for a set of covariates and performing inverse-normalization on the phenotypic scores (custom R script)
	- 2d. Filtering for related subjects and other characteristics: scripts to produce the final files that will be the input to the GWAS.
3) **Executing GWAS**: self-explanatory. Currently, we support Plink and BGENIE.
4) **Compile results and generate figures**: compile the GWAS results and generate Manhattan plots and Q-Q plots. You can also generate LocusZoom plots.
5) **Downstream analysis**: Integrate data from other sources in order to interpret the results. We currently support: proximity analysis using `biomaRt`, gene ontology term enrichment with `g:Profiler`, transcriptome-wide association studies with `S-PrediXcan` and pleiotropy analysis using the IEU GWAS Open Project.

Steps 2 to 4 rely on a single `yaml` configuration file.

## Usage

#### Fetching data
Instructions on how to download each kind of genetic data can be found [in this link](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/ukbgene_instruct.html).

#### Pre-processing data
The script that performs this task is `src/preprocess_files_for_GWAS.R`.
Example of usage (GWAS on left-ventricular end-diastolic volume, run on unrelated British subjects using Plink):

```
Rscript src/preprocess_files_for_GWAS \
  --phenotype_file data/phenotypes/cardiac_phenotypes/lvedv.csv
  --phenotypes LVEDV
  --columns_to_exclude id 
  --samples_to_include data/ids_list/british_subjects.txt
  --samples_to_exclude data/ids_list/related_british_subjects.txt
  --covariates_config_yaml config_files/standard_covariates.yml
  --output_file output/lvedv_adjusted_british.tsv
  --gwas_software plink
  --overwrite_output
```

#### Executing GWAS
One needs to execute the command

`python main.py -c <YAML_FILE>`

The complexity is located in the YAML configuration file.

#### Post-processing

This script takes a set of BGENIE output files, one for each region, possibly containing the summary statistics of many different phenotypes, and compiles them into one file per phenotype (which spans the whole genome):
```bash
Rscript src/postprocessing/concatenate_gwas.R \
--input_results_dir gwas_output/by_region \
--output_filename_pattern output/GWAS__{phenotype} \
--phenotypes LVEDV
```

The following sorts the GWAS summary statistics by chromosome and position. This can be convenient for some applications (for example LocusZoom plots, or if the file wants be indexed):
```bash
bash src/postprocessing/sort_by_chr_and_pos.sh gwas_output/GWAS_LVEDV.tsv
```

The following command generates Manhattan plots, Q-Q plots and per-region summaries (best association for each ~2Mb region):
```bash
Rscript src/postprocessing/gwas_analysis.R \
--gwas_folder gwas_output/cardiac_indices \
--gwas_pattern GWAS_{phenotype} \
--phenotypes LVEDV LVESV LVM LVEF LVSV \
--cache_rds
```

You can set the `--debug` flag if you want additional output.

## Configuration file

    chromosomes: 20-22  
    data_dir: <>  
    output_dir: <>  
    individuals: <>  
    filename_patterns: {
        genotype: {
            bed: <>,
            bim: <>, 
            fam: <>
        },
        phenotype: {
            phenotype_file: <>,
            phenotype_file_tmp: <>,
            phenotype_list: <>,   
            covariates: <>
        },
        gwas: "gwas_output/{phenotype}/gwas__{phenotype}"
    }

- `chromosomes` (optional): it consists of comma-separated ranges, where a range is a single number (e.g. `1`) or a proper range (e.g. `3-7`).
- `data_dir` (optional): directory where the input data is stored. If provided and the genotype and phenotype file paths are relative paths, they are interpreted as hanging from `data_dir`.
- `output_dir` (optional): directory where the output data will be stored.
- `individuals` (optional): path to a file containing the subset of individuals on which to run GWAS, one ID per line.
- `filename_patterns`: rules that determine the paths of the input and output files.
  - `genotype`:
    - `bed` (required): file name pattern containing that may contain the substring `{chromosome}` in which case it is replaced by the actual chromosome number. 
    - `bim` (required): idem previous.
    - `fam` (required): idem previous.
  - `phenotype`:
    - `phenotype_file` (required): file containing the phenotypes to run GWAS on.
    - `phenotype_file_delim` (optional, default: `,`): file delimiter.  
    - `phenotype_tmp_file` (optional): name of a temporary file that will be used as input to the GWAS software if the previous does not match the expected format.
    - `covariates` (not used yet).
  - `gwas` (required): name pattern for output files. It can contain the fields `{phenotype}` and `{suffix}`.
  - `plink` (optional, default: `plink`): path to the BGENIE executable.
