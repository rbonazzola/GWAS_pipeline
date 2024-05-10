library(tidyverse)

source("~/01_repos/GWAS_pipeline/src/postprocessing/retrieve_gene_positions.R")
g <- glue::glue

INTRON_DIR = "/home/rodrigo/01_repos/GWAS_pipeline/output/coma/finetuned_mlruns/SPrediXcan/best_intron_associations"
heart_assocs = "{INTRON_DIR}/joined_heart_intron_associations.csv" %>% g %>% read_csv
artery_assocs = "{INTRON_DIR}/joined_artery_intron_associations.csv" %>% g %>% read_csv

# genes_hg38_df <- get_gene_positions(reference_genome=38)
# genes_hg38_df <- genes_hg38_df %>% filter(gene_biotype == "protein_coding")

# parse_intron_name <- function(intron_name) {
#   kk = strsplit(intron_name, split = "_")[[1]]
#   data.frame(
#     chromosome=kk[2], 
#     start=as.integer(kk[3]), 
#     end=as.integer(kk[4])
#   )
# }

# This seems fine but sometimes I get more than one gene, whereas sometimes I get no gene
# Therefore, I use the mappings from GTEx files
# intron2gene <- function(intron) {
#   gene_row = genes_hg38_df %>% filter(intron['chromosome'] == chromosome_name & intron['start'] > start_position & intron['end'] < end_position)
#   print(intron)
#   intron_name <- "{intron['chromosome']}_{as.integer(intron['start'])}_{as.integer(intron['end'])}" %>% g
#   gene_row$hgnc_symbol
#   # names(kk) = intron_name
# }
# to_intron_name = function(row) paste(c("intron", row[1], row[2], row[3]), collapse="_")

#dd <- bind_rows(
  # lapply(unique(heart_assocs$gene_name), function(intron) heart_assocs[intron])
# )


# genes = apply(dd, 1, intron2gene)
# 
# intron_names = apply(dd, 1, to_intron_name) %>% gsub(" ", "", .)
# 
# intron2gene_dict <- as.list(genes)
# names(intron2gene_dict) <- intron_names
# 
# gene_names_in_order <- sapply(
#   heart_assocs$gene, 
#   function(intron) 
#     paste(intron2gene_dict[[gene]], collapse=", ")
# ) %>% gsub(", ,", ", ", .)

# $gene_name <- gene_names_in_order

SPREDIXCAN_DIR <- "/home/rodrigo/01_repos/GWAS_pipeline/output/coma/finetuned_mlruns/SPrediXcan/"

heart_sgene_mapping <- read_csv("{SPREDIXCAN_DIR}/intron2gene_mapping.csv" %>% g, col_names = c("intron", "gene"))
heart_sgenes <- heart_sgene_mapping$gene
names(heart_sgenes) <- gsub("chr", "intron_", heart_sgene_mapping$intron)

artery_sgene_mapping <- read_csv("{SPREDIXCAN_DIR}/intron2gene_mapping.csv" %>% g, col_names = c("intron", "gene"))
artery_sgenes <- artery_sgene_mapping$gene
names(artery_sgenes) <- gsub("chr", "intron_", artery_sgene_mapping$intron)

heart_spred_hits <- heart_assocs %>% mutate(gene_name=heart_sgenes[gene])  
heart_spred_hits <- heart_spred_hits %>% group_by(tissue, gene) %>% filter(pvalue == min(pvalue))
heart_spred_hits <- heart_spred_hits %>% filter(pvalue < 1e-8) %>% ungroup 
heart_spred_hits <- heart_spred_hits %>% dplyr::select(-pred_perf_r2, -pred_perf_pval, -pred_perf_qval) 
heart_spred_hits <- heart_spred_hits %>% .[order(.$pvalue),]

# gene_names_in_order <- sapply(artery_assocs$gene, function(gene) paste(intron2gene_dict[[gene]], collapse=", ")) %>% gsub(", ,", ", ", .)
# artery_assocs$gene_name <- gene_names_in_order

artery_spred_hits <- artery_assocs %>% mutate(gene_name=artery_sgenes[gene])   
artery_spred_hits <- artery_spred_hits %>% group_by(tissue, gene) %>% filter(pvalue == min(pvalue))
artery_spred_hits <- artery_spred_hits %>% filter(pvalue < 1e-8) %>% ungroup 
artery_spred_hits <- artery_spred_hits %>% dplyr::select(-pred_perf_r2, -pred_perf_pval, -pred_perf_qval) 
artery_spred_hits <- artery_spred_hits %>% .[order(.$pvalue),]

rbind(artery_spred_hits, heart_spred_hits) %>% 
  write.table(
    "{SPREDIXCAN_DIR}/best_assocs_HeartArtery_tissues.csv" %>% g,
    row.names = FALSE, quote = FALSE, sep = ","
  )