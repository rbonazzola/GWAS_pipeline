library(tidyverse)
library(logging)
library(argparse)

parser <- ArgumentParser()

# phenotype_file = "/home/rodrigo/01_repos/CardiacMotionRL/mlruns/2/8c1ffa20cacc4b6c88e18159e01867b4/artifacts/output/latent_vector_adj-all_GBR.csv"
# Phenotype file
parser$add_argument("-p", "--phenotype_file", required=TRUE, help="Path to the (original) phenotype file. This and all paths must be either absolute paths or paths relative to this repo's root directory")
parser$add_argument("--phenotypes", default=NULL, nargs="+")
parser$add_argument("--relatedness_file", default="/home/rodrigo/01_repos/GWAS_pipeline/data/datasets/ukb11350_rel_s488226.dat")

# Samples
# parser$add_argument("--samples_file", default="/home/rodrigo/tmp/example.sample")
parser$add_argument("--bgen_sample_file", default="~/01_repos/GWAS_pipeline/data/datasets/ids_list/ukb22828_c1_b0_v3_s487159.sample", help="Your sample file, it should be the one that the BGEN files are linked to")
parser$add_argument("--output_dir", default="~/01_repos/CardiacGWAS/results")
parser$add_argument("--tmpdir", default="tmp/")
parser$add_argument("--seed", default=42, nargs="+")
parser$add_argument("--gzip_output", default=FALSE, action="store_true", help="Whether to gzip the output files.")
parser$add_argument("--na_code", default="-999", help="Code representing NA.")

# Output
parser$add_argument("-o", "--output_file_prefix", required=TRUE)
parser$add_argument("--overwrite_output", default=FALSE, action="store_true", help="Flag indicating if this script should be re-run upon finding a previously generated file with the same name as --output_file.")

args <- parser$parse_args()

g <- glue::glue

SEEDS = args$seed
# SAMPLES_FILE <- "~/01_repos/GWAS_pipeline/data/datasets/ids_list/ukb22828_c1_b0_v3_s487159.sample"
# SAMPLES_FILE <- "~/01_repos/GWAS_pipeline/data/datasets/ids_list/all_cmr_ids_bgenie.txt"
samples_file <- args$bgen_sample_file
samples_df <- readr::read_table(samples_file, skip = 2, col_names = c("id_1", "id_2", "missing"))
samples_df <- as.data.frame(samples_df)
names(samples_df)[1] <- "ID"
samples_df <- mutate(samples_df, ID=as.integer(ID))
all_ids <- samples_df$ID

rel_file <- args$relatedness_file
rel_df <- read_delim(rel_file)

fracs_replication <- c(0.0, 0.10) # , 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5)

# GreedyExec <- "~/doctorado/software/GreedyRelated/bin/GreedyRelated"
GreedyExec <- "GreedyRelated"

partition_suffix_pat <- "{100*frac_replication}n{100*(1-frac_replication)}" 

if (!file.exists(args$tmpdir)) {
  logging::loginfo(g("Folder {args$tmpdir} does not exist. Creating..."))
  dir.create(args$tmpdir)
}

for (SEED in SEEDS) {
  
  disc_sizes <- numeric()
  repl_sizes <- numeric()
  
  dir.create(g("{args$output_dir}/seed_{SEED}/"))
  OUTPUT_FP_DISC = "{args$output_dir}/seed_{SEED}/{args$output_file_prefix}-discovery_{partition_suffix}__extended_to_{sample_file_n}.csv"  
  OUTPUT_FP_REPL = "{args$output_dir}/seed_{SEED}/{args$output_file_prefix}-replication_{partition_suffix}__extended_to_{sample_file_n}.csv"  
  
  logging::loginfo("Using seed {SEED} for discovery/replication split" %>% g)
  
  for (frac_replication in fracs_replication) {
    
    set.seed(SEED)
    
    ids_discovery <- samples_df %>% sample_frac(size = 1-frac_replication) %>% .$ID
    
    dir.create(g("{args$tmpdir}/GreedyRelated/"))
    subset_rel1_df <- rel_df %>% filter(ID1 %in% ids_discovery) %>% filter(ID2 %in% ids_discovery)
    subset_rel1_df <- subset_rel1_df %>% dplyr::select(-HetHet, -IBS0)
    intermediate_rel1 <- "{args$tmpdir}/GreedyRelated/greedy_frac_{1-frac_replication}.dat" %>% g
    write_delim(subset_rel1_df, intermediate_rel1, delim="\t")
    
    exclude_subjects_file <- "{args$tmpdir}/GreedyRelated/exclude_subjects_{1-frac_replication}.txt" %>% g
    
    greedy_related_cmd1 <- "{GreedyExec} -r {intermediate_rel1} --seed 1 --id1 ID1 --id2 ID2 -f Kinship -o {exclude_subjects_file}" %>% g
    system(greedy_related_cmd1)
    
    # The file does not have a header (hence col_names=FALSE)
    exclude_subjects <- read_delim(exclude_subjects_file, col_names=FALSE)$X1
    
    ids_discovery <- setdiff(ids_discovery, exclude_subjects)
    ids_replication <- setdiff(all_ids, ids_discovery)
    
    subset_rel2_df <- rel_df %>% filter(ID1 %in% ids_replication) %>% filter(ID2 %in% ids_replication)
    subset_rel2_df <- subset_rel2_df %>% dplyr::select(-HetHet, -IBS0)
    intermediate_rel2 <- "{args$tmpdir}/GreedyRelated/greedy_frac_{1-frac_replication}.dat" %>% g
    write_delim(subset_rel2_df, intermediate_rel2, delim="\t")
    
    exclude_subjects_file2 <- "{args$tmpdir}/GreedyRelated/exclude_subjects_{1-frac_replication}.txt" %>% g
    
    greedy_related_cmd2 <- "{GreedyExec} -r {intermediate_rel2} --seed 1 --id1 ID1 --id2 ID2 -f Kinship -o {exclude_subjects_file2}" %>% g
    system(greedy_related_cmd2)
    
    exclude_subjects <- read_delim(exclude_subjects_file2, col_names=FALSE)$X1
    ids_replication <- setdiff(ids_replication, exclude_subjects)
    
    disc_sizes <- c(disc_sizes, length(ids_discovery))
    repl_sizes <- c(repl_sizes, length(ids_replication))
    
    partition_suffix <- partition_suffix_pat %>% g
    
    write_csv(
      data.frame(ids_discovery),
      "{args$tmpdir}/GreedyRelated/ids_discovery_{partition_suffix}.txt" %>% g
    )
    
    write_csv(
      data.frame(ids_replication),
      "{args$tmpdir}/GreedyRelated/ids_replication_{partition_suffix}.txt" %>% g
    )
    
  }
  
  # print(disc_sizes); print(repl_sizes)
  
  ################################################################################
    
  phenotype_file = args$phenotype_file
  phenotype_df = read.delim(phenotype_file, sep = "\t")
  print(phenotype_df)
  phenotype_df[duplicated(phenotype_df$ID),] = args$na_code
  phenotype_df[duplicated(phenotype_df$ID), "ID"] = -100 * (1:sum(duplicated(phenotype_df$ID)))
  phenotype_df$ID <- as.character(phenotype_df$ID)
  rownames(phenotype_df) <- phenotype_df$ID
  
  print(length(unique(rownames(phenotype_df))))
  
  SAMPLES_FILES <- c(
    # "63000"="~/01_repos/GWAS_pipeline/data/datasets/ids_list/all_cmr_ids_bgenie.txt"
    nrows(readr::read_table(args$bgen_sample_file, skip = 2))=args$bgen_sample_file
    # "487000"="~/01_repos/GWAS_pipeline/data/datasets/ids_list/ukb22828_c1_b0_v3_s487159.sample"
  )
  
  for (sample_file_n in names(SAMPLES_FILES)) {
    
    sample_file <- SAMPLES_FILES[sample_file_n]
    samples_df <- readr::read_table(sample_file, skip = 2, col_names = c("id_1", "id_2", "missing"))
    samples_df <- as.data.frame(samples_df)
    names(samples_df)[1] <- "ID"
    samples_df <- mutate(samples_df, ID=as.character(ID))
    
    for (frac_replication in fracs_replication) {
      
      partition_suffix <- partition_suffix_pat %>% g
      ids_discovery <- read.delim("{args$tmpdir}/GreedyRelated/ids_discovery_{partition_suffix}.txt" %>% g)[,1]
      ids_replication <- read.delim("{args$tmpdir}/GreedyRelated/ids_replication_{partition_suffix}.txt" %>% g)[,1]
      
      print("{frac_replication}: {length(ids_replication)}+{length(ids_discovery)}" %>% g)
      
      ids_discovery <- as.character(ids_discovery)
      ids_replication <- as.character(ids_replication)
      
      df_discovery = phenotype_df[ids_discovery,]
      df_replication = phenotype_df[ids_replication,]
      
      pheno_df <- samples_df %>% left_join(df_discovery, by = "ID") # %>% select("ID", all_of(pheno_names)) # select(.dots = pheno_names)
      pheno_df <- pheno_df %>% dplyr::select(-id_2, -missing)
      pheno_df <- cbind(pheno_df$ID, pheno_df %>% dplyr::select(-ID) %>% round(digits=4))
      
      ofile_discovery <- OUTPUT_FP_DISC %>% g
      print(head(pheno_df))
      readr::write_delim(
          pheno_df, 
          ofile_discovery, 
          col_names = TRUE, delim = "\t", na = args$na_code
      )
      
      if (args$gzip_output) {
        system(g("gzip {ofile_discovery}"))
      }
      
      pheno_df <- samples_df %>% left_join(df_replication, by = "ID") # %>% select("ID", all_of(pheno_names)) # select(.dots = pheno_names)
      pheno_df <- pheno_df %>% dplyr::select(-id_2, -missing)
      pheno_df <- cbind(pheno_df$ID, pheno_df %>% dplyr::select(-ID) %>% round(digits=4))
      
      ofile_replication <- OUTPUT_FP_REPL %>% g
      readr::write_delim(
          pheno_df, 
          ofile_replication, 
          col_names = TRUE, delim = "\t", na = args$na_code, 
      )
      
      if (args$gzip_output) {
        system(g("gzip {ofile_replication}"))
      }
      
    }
  }
}
