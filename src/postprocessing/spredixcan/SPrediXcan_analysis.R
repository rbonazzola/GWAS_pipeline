library(dplyr)

spred_eqtl <- read_delim("~/01_repos/CardiacMotionGWAS/results/SPrediXcan/eqtl_best.csv")
spred_eqtl <- spred_eqtl %>% .[order(.$pvalue),]
spred_eqtl <- spred_eqtl %>% filter(pvalue < 1e-10)
spred_eqtl %>% filter(tissue %in% c("HRTAA", "HRTLV", "ARAORT"))
# spred_eqtl$gene %>% unique

genes <- spred_eqtl$gene %>% unique %>% strsplit("\\.") %>% sapply(function(x) x[1])
# _name

genes_df %>% filter(ensembl_gene_id %in% genes) %>% .[order(as.integer(.$chromosome), .$start_position),]