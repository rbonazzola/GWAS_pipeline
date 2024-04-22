library(dplyr)
library(ggplot2)
library(glue)
library(argparse)
options(error = recover)
##################################################################################################################

main <- function(args) {
  
  print(getwd())
  
  long_range_LD_file <- glue(args$long_range_ld_file)
  
  bed_file_pattern <- args$bed_file_pattern
  bim_file_pattern <- args$bim_file_pattern
  fam_file_pattern <- args$fam_file_pattern
  
  ids_file <- args$ids_file
  
  extract_variants_in_range <- function(chromosome, start, end) {
    bim_file <- glue(bim_file_pattern)
    snp_df <- read.table(bim_file, header=FALSE)
    colnames(snp_df) <- c("chromosome", "rsid", "pos_cm", "position", "allele_1", "allele_2")
    snp_df <- snp_df %>% filter(position > start & position < end)
    snp_df$rsid
  }
  
  extract_variants_ambiguous_strand <- function(chromosome) {
    bim_file <- glue(bim_file_pattern)
    snp_df <- read.table(bim_file, header=FALSE)
    colnames(snp_df) <- c("chromosome", "rsid", "pos_cm", "position", "allele_1", "allele_2")
    snp_df <- snp_df %>% filter(allele_1 == "C" & allele_2 == "G" | allele_1 == "G" & allele_2 == "C" | allele_1 == "A" & allele_2 == "T" | allele_1 == "T" & allele_2 == "A")
    as.character(snp_df$rsid)
  }
  
  get_snps_below_maf_thr <- function(chromosome, maf_threshold=0.025) {
    maf_file <- glue(maf_file_pattern)
    print(maf_file)
    maf_df <- read.table(maf_file, header=TRUE)
    maf_df <- maf_df %>% filter(MAF < maf_threshold)
    as.character(maf_df$SNP)
  }
  
  get_snps_above_miss_thr <- function(chromosome, missingness_threshold=0.015) {
    missingness_file <- glue(missingness_file_pattern)
    missingness_df <- read.table(missingness_file, header=TRUE)
    missingness_df <- missingness_df %>% filter(F_MISS > missingness_threshold)
    as.character(missingness_df$SNP)
  }
  
  get_snps_below_hwe_p_thr <- function(chromosome, hwe_p_threshold=1e-5) {
    hwe_file <- glue(hwe_file_pattern)
    hwe_df <- read.table(hwe_file, header=TRUE)
    hwe_df <- hwe_df %>% filter(P < hwe_p_threshold)
    as.character(hwe_df$SNP)
  }
  
  long_range_ld_df <- read.table(long_range_LD_file, sep = "\t", header=TRUE)
  long_range_ld.lst <- long_range_ld_df %>% select(1:3) %>% split(., seq(nrow(.))) # dataframe to list by rows
  snps_in_ld <- unlist(lapply(long_range_ld.lst, function(x) {x <- as.integer(x); extract_variants_in_range(x[1], x[2], x[3])}))
  
  ambiguous_strand_snps <- unlist(lapply(1:22, extract_variants_ambiguous_strand))
  
  MAF_COMMAND_PATTERN <- "plink --keep {ids_file} --freq {bedbimfam} --out {out_bfile}"
  MISSING_COMMAND_PATTERN <- "plink --keep {ids_file} --missing {bedbimfam} --out {out_bfile}"
  HWE_COMMAND_PATTERN <- "plink --keep {ids_file} --hardy {bedbimfam} --out {out_bfile}"
  PLINK_LD_PRUNE_COMMAND_PATTERN <- "plink --exclude {snps_to_exclude_file} --indep-pairwise 100 10 0.1 {bedbimfam} --make-bed --out {out_bfile}"
  
  maf_file_pattern <- paste0(out_bfile_pattern, ".frq")
  missingness_file_pattern <- paste0(out_bfile_pattern, ".lmiss")
  prunein_file_pattern <- paste0(out_bfile_pattern, ".prune.in")
  hwe_file_pattern <- paste0(out_bfile_pattern, ".hwe")
  prunein_snps_file <- glue("{transforms_dir}/snps_for_pca_prune.in")
  
  for (chromosome in 1:22) {

    bedfile <- glue(bed_file_pattern)
    bimfile <- glue(bim_file_pattern)
    famfile <- glue(fam_file_pattern)
    bedbimfam <- glue("--bed {bedfile} --bim {bimfile} --fam {famfile}")
    
    out_bfile <- glue(out_bfile_pattern)
    plink_command_maf <- glue(MAF_COMMAND_PATTERN)
    plink_command_missing <- glue(MISSING_COMMAND_PATTERN)
    plink_command_hwe <- glue(HWE_COMMAND_PATTERN)
    
    print(plink_command_maf)
    system(plink_command_maf)
    
    print(plink_command_missing)
    system(plink_command_missing)
    
    print(plink_command_hwe)
    system(plink_command_hwe)
  }
  
  snps_below_maf_thr <- unlist(lapply(1:22, get_snps_below_maf_thr))
  snps_above_miss_thr <- unlist(lapply(1:22, get_snps_above_miss_thr))
  snps_below_hwe_p_thr <- unlist(lapply(1:22, get_snps_below_hwe_p_thr))
  
  logging::loginfo(glue("SNPs with ambiguous strand: {length(ambiguous_strand_snps)}"))
  logging::loginfo(glue("SNPs belonging to long-range LD regions found: {length(snps_in_ld)}"))
  logging::loginfo(glue("SNPs with MAF > 2.5% found: {length(snps_below_maf_thr)}"))
  logging::loginfo(glue("SNPs with missingness > 1.5% found: {length(snps_above_miss_thr)}"))
  logging::loginfo(glue("SNPs with HWE p-value below 1e-5: {length(snps_below_hwe_p_thr)}"))
  
  snps_to_exclude <- unique(c(ambiguous_strand_snps, snps_in_ld, snps_below_maf_thr, snps_above_miss_thr, snps_below_hwe_p_thr))
  snps_exclude_df <- as.data.frame(snps_to_exclude)
  logging::loginfo(glue("Total number of SNPs to exclude: {length(snps_to_exclude)}"))
  
  write.csv(snps_exclude_df, snps_to_exclude_file, quote = FALSE, row.names = FALSE)
  
  for(chromosome in 1:22) {
    bedfile <- glue(bed_file_pattern)
    bimfile <- glue(bim_file_pattern)
    famfile <- glue(fam_file_pattern)
    bedbimfam <- glue("--bed {bedfile} --bim {bimfile} --fam {famfile}")

    out_bfile <- glue(out_bfile_pattern)
    PLINK_LD_PRUNE_COMMAND <- glue(PLINK_LD_PRUNE_COMMAND_PATTERN)
    print(PLINK_LD_PRUNE_COMMAND)
    system(PLINK_LD_PRUNE_COMMAND)
  }
  
  prunein_snps <- unlist(lapply(1:22, function(chromosome) read.table(glue(prunein_file_pattern))[,1]))
  prunein_snps_df <- as.data.frame(prunein_snps)
  write.csv(prunein_snps_df, prunein_snps_file, quote = FALSE, row.names = FALSE)
  
  if (!file.exists(files_to_merge)) {
    for (chromosome in 1:22) {
      file <- glue(out_bfile_pattern)
      print(glue("echo {file} >> {files_to_merge}"))
      system(glue("echo {file} >> {files_to_merge}"))
    }
  }
  
  merge_command <- "plink --bfile {geno_file_for_pca_chr1} --merge-list {files_to_merge} --make-bed --out {geno_file_for_pca}"
  
  if (!file.exists(geno_file_for_pca)) {
    print(glue(merge_command))
    system(glue(merge_command))
  }
  
  extract_id <- function(x) strsplit(x, ":")[[1]][1]
  
  if (!file.exists(genomic_pca_file)) {
    genomic_pca <- flashpcaR::flashpca(geno_file_for_pca, ndim=10)
    pca_proj <- genomic_pca$projection
    pca_proj_df <- as.data.frame(pca_proj)
    colnames(pca_proj_df) <- paste0("PC", 1:ncol(pca_proj_df))
    pca_proj_df$ID <- rownames(pca_proj_df)
    pca_proj_df <- cbind(pca_proj_df %>% select(ID), pca_proj_df %>% select(-ID)) # reorder columns
    write.table(pca_proj_df, file=genomic_pca_file , quote=FALSE, row.names=FALSE, sep = "\t")
  } else {
    pca_proj_df <- read.table(genomic_pca_file, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  pca_proj_df$ID <- unname(sapply(pca_proj_df$ID, extract_id))
  rownames(pca_proj_df) <- pca_proj_df$ID
  
}

dataset_dir <- "data/datasets"
intermediate_dir <- "data/intermediate"
transforms_dir <- "data/transforms"

# Create argument parser
parser <- ArgumentParser(description = "Perform genetic PCA on a set of subjects using BED/BIM/FAM files.")

parser$add_argument("--bed_file_pattern", 
                    help = "Path to the .bed file", default="{dataset_dir}/calls/ukb22418_c{chromosome}_b0_v2.bed")

parser$add_argument("--bim_file_pattern", 
                    help = "Path to the .bim file", default="{dataset_dir}/calls/ukb_snp_chr{chromosome}_v2.bim")

parser$add_argument("--fam_file_pattern", 
                    help = "Path to the .fam file", default="{dataset_dir}/calls/ukb22418_c{chromosome}_b0_v2_s488170.fam")

parser$add_argument("--long_range_ld_file", 
                    help = "Path to the long-range LD SNPs file", default="{dataset_dir}/long_range_LD_regions.bed")

parser$add_argument("--ids_file", 
                    help = "Path to the file including IDs of subjects to perform PCA for (in PLINK format).", 
                    default="{dataset_dir}/ids_list/all_OCT_GBR_plink.txt")

parser$add_argument("--output_pc_file", 
                    help = "Output path of the genetic PC file.", 
                    default="{transforms_dir}/genomicPCs_GBR.tsv")

args <- parser$parse_args()

files_to_merge <- glue("{transforms_dir}/GenomicPCA_OCT/files_to_merge.txt")
snps_to_exclude_file <- glue("{transforms_dir}/GenomicPCA_OCT/exclude_snps_for_pca.txt")

out_bfile_pattern <- "{transforms_dir}/GenomicPCA_OCT/genotypes/ukb_cal_chr{chromosome}_v2_GBR_indiv"

geno_file_for_pca_for_chr <- "{transforms_dir}/GenomicPCA_OCT/genotypes/ukb_cal_chr{chromosome}_v2_GBR_indiv"
geno_file_for_pca_chr1 <- { chromosome=1; glue(geno_file_for_pca_for_chr) }
geno_file_for_pca <- glue("{transforms_dir}/GenomicPCA_OCT/ukb_cal_all_chrs_v2_GBR_indiv")

long_range_LD_file <- glue(args$long_range_ld_file)
genomic_pca_file <- glue(args$output_pc_file)

main(args)

##################################################################################################################




