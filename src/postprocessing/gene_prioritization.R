suppressPackageStartupMessages({
  #install.packages("gprofiler2")
  # library(gprofiler2)
  library(tidyverse)
  library(biomaRt)
  library(logging)
  # install.packages("argparse")
  library(argparse)
})

get_regions_df <- function(path=LD_INDEP_REGIONS) {
  regions <- read.delim(path, sep = "\t")
  regions <- regions %>% group_by(chr) %>% mutate(id = paste0(row_number())) %>% as.data.frame
  regions$chr <-  sub("\\s+$", "", regions$chr)
  regions <- regions %>% mutate(id=paste(chr, id, sep = "_")) %>% ungroup() 
  regions$chr <- regions$chr %>% gsub("chr", "", .) %>% as.numeric()
  regions
}


get_gene_positions <- function() {
  
  selected_attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype", "strand")
  ensembl <- biomaRt::useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", verbose=TRUE)
  dataset <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  pp <- biomaRt::getBM(attributes = selected_attributes, mart = ensembl)
  
  pp <- pp %>% dplyr::filter(chromosome_name %in% as.character(1:22))
  
  # pp <- pp %>% filter(gene_biotype == "protein_coding")
  pp
  
}


get_genes_around_snps <- function(genes_df, snps_df, window_size) {
  genes_around_gwas_hits <- merge(genes_df, snps_df, by.x = "chromosome_name", by.y = "CHR")
  genes_around_gwas_hits <- genes_around_gwas_hits %>% mutate(distance_from_tss=abs(BP-start_position))
  genes_around_gwas_hits[genes_around_gwas_hits$strand==-1, "distance_from_tss"] <- genes_around_gwas_hits %>% filter(strand==-1) %>% {abs(.$BP-.$end_position)}
  genes_around_gwas_hits <- genes_around_gwas_hits %>% filter(abs(BP-start_position) < window_size/2)
  genes_around_gwas_hits
}


get_genes_for_term <- function(term_name) {
  ensgs <- gprofiler2::gost(term_name)$meta$genes_metadata$query$query_1$ensgs
  genes_df %>% filter(ensembl_gene_id %in% ensgs)
}

main <- function(args) {
  
  snps_df <- read.csv(args$snps_file)
  genes_df <- get_gene_positions()
  
  logging::loginfo(glue::glue("Window of size: {args$window_size}"))
  genes_around_gwas_hits <- get_genes_around_snps(genes_df, snps_df, args$window_size)
  
  logging::loginfo(glue::glue("Writing nearby genes to {args$output_nearby_genes}"))
  if (args$only_protein_coding) {
    genes_around_gwas_hits <- genes_around_gwas_hits %>% filter(gene_biotype=="protein_coding")
  }
  genes_around_gwas_hits <- genes_around_gwas_hits %>% filter(!duplicated(.))
  write_csv(genes_around_gwas_hits, file=args$output_nearby_genes)

  snps_inside_genes <- genes_around_gwas_hits %>% arrange(chromosome_name, start_position) %>% 
      mutate("is_inside"=(sign(BP-start_position)*sign(BP-end_position) < 0)) %>%
      filter(is_inside) 
  
  logging::loginfo(glue::glue("Writing SNPs that are inside genes to {args$snps_in_genes_file}"))
  write_csv(snps_inside_genes, file=args$snps_in_genes_file)
}


###################################################################################################################

# sheet_name <- "loci_mapping"
# spreadsheet_url <- glue::glue("https://docs.google.com/spreadsheets/d/1XvVDFZSvcWWyVaLaQuTpglOqrCGB6Kdf6c78JJxymYw/gviz/tq?tqx=out:csv&sheet={sheet_name}")
# loci_data <- googlesheets4::read_sheet(spreadsheet_url)
# loci_data <- loci_data %>% filter(is.na(duplicated))

parser <- ArgumentParser(description = "A simple R script with argparse")
parser$add_argument("--snps_file", help = "Input SNP file, with columns \"CHR\", \"SNP\", \"BP\", \"region\". ", default="~/tmp/OCT_snps.txt")
parser$add_argument("--window_size", "--window-size", help = "Window size", default=1e6, type="numeric")
parser$add_argument("--only_protein_coding", "--only-protein-coding", action="store_true", help = "Whether to include only protein-coding gene.")
parser$add_argument("--output_nearby_genes", help = "Path to output file with nearby genes.", default="~/tmp/OCT_all_snps/genes_arounds_gwas_hits.csv")
parser$add_argument("--snps_in_genes_file", help = "Path to file.", default="~/tmp/OCT_all_snps/gwas_hits_inside_genes.csv")

args <- parser$parse_args()

main(args)


# TERM_REGEX <- "cardiac|heart|muscle|sarco|calcium|contractile|myo|Z disc|I band|atri|cardi"
# 
# go_results <- gprofiler2::gost(genes_around_gwas_hits$hgnc_symbol)
# 
# relevant_term_ids <- go_results$result %>% .[grepl(TERM_REGEX, .$term_name), "term_id"]

################################################################################################################

# genes_for_terms <- lapply(relevant_term_ids, get_genes_for_term)
# names(genes_for_terms) <- relevant_term_ids

# get_genes_for_term(term_name) %>% filter(ensembl_gene_id %in% genes_around_gwas_hits$ensembl_gene_id)

# genes_in_term <- sapply(
#   genes_for_terms, 
#   function(x) {
#     intersect(x$ensembl_gene_id, genes_around_gwas_hits$ensembl_gene_id)
#   }
# )
# 
# candidate_genes <- genes_df %>% 
#   filter(ensembl_gene_id %in% unlist(genes_in_term)) %>% 
#   mutate(chromosome_name=as.numeric(chromosome_name)) 
# 
# # include distance to SNP
# merge(
#   candidate_genes %>% dplyr::select(-gene_biotype), 
#   genes_around_gwas_hits %>% dplyr::select(ensembl_gene_id, SNP, BP, region, distance_from_tss), 
#   by="ensembl_gene_id"
# ) %>% 
#   arrange(chromosome_name, start_position) %>% 
#   mutate("is_inside"=(sign(BP-start_position)*sign(BP-end_position) < 0)) %>%
#   filter(is_inside)
# 
# heart_development_go_term <- "GO:0007507"
# genes_around_gwas_hits %>% filter(ensembl_gene_id %in% get_genes_for_term(heart_development_go_term)$ensembl_gene_id) %>% 
#   mutate("is_inside"=(sign(BP-start_position)*sign(BP-end_position) < 0)) %>%
#   filter(is_inside) %>% 
#   arrange(as.numeric(chromosome_name), start_position)
# 
# 
# #  Export table to LaTex
# go_df = go_results$result %>% dplyr::select(-query, -source, -parents, -significant, -precision, -recall, -source_order)
# print(xtable(go_df, type="latex"), file="go_term_enrichment.txt")
# 
# gconvert(genes_around_gwas_hits$ensembl_gene_id, target="HGNC") %>% { . = .; names(.) = .$name}# 