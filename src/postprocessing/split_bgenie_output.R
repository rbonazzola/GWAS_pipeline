suppressPackageStartupMessages({
  library(tidyverse)
})

# Changes the output of BGENIE into a more standard format
modify_one_pheno_df <- function(df) {
  
  phenotype <- colnames(df) %>% strsplit(split="_|\\.") %>% .[[1]] %>% .[c(1,2)] %>% paste0(collapse="_")
  
  beta_column <- paste0(phenotype, "_beta")
  se_column <- paste0(phenotype, "_se")
  t_column <- paste0(phenotype, "_t")
  log10p_column <- paste0(phenotype, ".log10p")
  
  col_names <- c(beta_column, se_column, t_column, log10p_column)
  colnames(df) <- c("BETA", "SE", "T", "LOG10P")
  df <- mutate(df, P=10^(-LOG10P))
  df$LOG10P <- NULL
  df
}

create_file <- function(dfs, phenotype, region) {
  
  print(phenotype)
  df <- dfs[[phenotype]]
  folder <- glue::glue("{results_dir}/{phenotype}")
  
  if (!file.exists(folder)) {
    # If the folder doesn't exist, create it
    dir.create(folder)
  }
  
  filename <- glue::glue("{folder}/{region}.tsv")
  write_tsv(df, filename)
}

# for (filename in list.files(GWAS_BY_REGION_DIR) %>% .[grepl("txt.gz", .)]) {

main <- function(results_dir, filename) {
    
    filename = glue::glue("{results_dir}/{filename}")

    region <- basename(filename) %>% gsub(".txt.gz", "", .)
    print(region)
    
    df = read.table(gzfile(filename), header=TRUE, sep=" ")
    
    df_first <- df %>% rename(CHR=chr, SNP=rsid, BP=pos, AF=af, INFO=info)
    
    snp_data_cols <- c("CHR", "SNP", "BP", "AF", "a_0", "a_1", "INFO")
    snps_data <- df_first %>% select(all_of(snp_data_cols))
    
    gwas_sumstats <- df_first %>% select(-any_of(snp_data_cols))
    
    # Each phenotype has 4 columns
    n_phenos <- ncol(gwas_sumstats) / 4 
    col_indices <- lapply(1:n_phenos, function(i) 4*(i-1)+1:4)
    
    dfs <- lapply(col_indices, function(indices) gwas_sumstats[,indices])
    
    pheno_names <- sapply(dfs, function(df) { colnames(df) %>% strsplit(split="_|\\.") %>% .[[1]] %>% .[c(1,2)] %>% paste0(collapse="_")}) 
    names(dfs) <- pheno_names
    
    pheno_names <- pheno_names[2:length(pheno_names)]
    dfs <- dfs[pheno_names]
    dfs <- lapply(dfs, modify_one_pheno_df)
    
    snps_folder <- glue::glue("{results_dir}/snps_info")
    if (!file.exists(snps_folder)) {
      # If the folder doesn't exist, create it
      dir.create(snps_folder)
    }
    write_tsv(snps_data, glue::glue("{snps_folder}/{region}__snps_data.tsv"))
    
    for (phenotype in names(dfs)) {
      folder <- glue::glue("{results_dir}/{phenotype}")
      if ( file.exists(glue::glue("{folder}/{region}.tsv")) ) {
          logging::loginfo(glue::glue("{region} has already been processed"))
          break
      }
      create_file(dfs, phenotype, region)
    }
}


parser <- argparse::ArgumentParser()
parser$add_argument("--input_results_dir")
parser$add_argument("--files", nargs="+")
args <- parser$parse_args()

results_dir = args$input_results_dir
for (filename in args$files) {
  main(args$input_results_dir, filename)
}
