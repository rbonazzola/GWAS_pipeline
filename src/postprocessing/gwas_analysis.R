suppressPackageStartupMessages({
  library(tidyverse)
  library(qqman)
  library(glue)
  library(argparse)
  library(logging)
})

# source("plots_and_others.R")

setwd(system("git rev-parse --show-toplevel", intern = TRUE))

fetch_args <- function() {

    parser <- ArgumentParser()
    parser$add_argument("--gwas_folder", nargs="+")
    parser$add_argument("--output_folder", default=NULL)
    parser$add_argument("--gwas_pattern", default="GWAS__{phenotype}", help="File pattern including the \"{phenotype}\" field, *without extensions*.")
    parser$add_argument("--phenotypes", nargs="+", default=sapply(0:7, function(x) paste0("z00", x)) )
    parser$add_argument("--title", default=FALSE, action="store_true")
    parser$add_argument("--cache_rds", action="store_true", default=FALSE)
    parser$add_argument("--overwrite_rds", action="store_true", default=FALSE)
    parser$add_argument("--no_manhattan", action="store_true", default=FALSE)
    parser$add_argument("--no_qqplot", action="store_true", default=FALSE)
    parser$add_argument("--no_summary", action="store_true", default=FALSE)
    parser$add_argument("--qqplot_pooled", action="store_true", default=FALSE)
    parser$add_argument("--color_odd_chr", default="mediumblue")
    parser$add_argument("--color_even_chr", default="lightblue")
    parser$add_argument("--figure_type", default="eps")
    parser$add_argument("--no_axis_labels",  action="store_true", default=FALSE)
    parser$add_argument("--debug", action="store_true", default=FALSE)
    
    args <- parser$parse_args()
    args
    
}
####################################################################################################

main <- function(args) {
  
    # TODO: put every file pattern into a configuration file
    plot_params <- list()
    plot_params$colors = c(args$color_odd_chr, args$color_even_chr)
    plot_params$no_axis_labels = args$no_axis_labels
    # color1 <- "hotpink4"
    # color2 <- "palevioletred2"
    
    which_plot <- list()
    which_plot$manhattan <- !args$no_manhattan
    which_plot$qqplot <- !args$no_qqplot
    which_plot$summary <- !args$no_summary
    
    if (args$debug) setLevel("DEBUG")

    gwas_folder <- args$gwas_folder
    output_dir <- ifelse(is.null(args$output_folder), args$gwas_folder, args$output_folder)  
    fp <- build_file_patterns(gwas_folder, output_dir)
    
    logging::logdebug(glue("Phenotypes to be processed: {args$phenotypes}"))
    logging::logdebug(fp)
    
    for (run_id in args$gwas_folder) {
      
        # TOFIX: run_id is not what I want here.
        logging::loginfo("Processing run with ID {run_id}" %>% glue)
      
        if ( !file.exists(glue::glue(fp$figs_dir))) { 
          dir.create(glue::glue(fp$figs_dir), showWarnings = TRUE)
        }
        
        pvals_pooled <- vector(length = 0)
        
        lapply(args$phenotypes, process_phenotype, file_patterns=fp, plot_params=plot_params, which_plot=which_plot)
        
        # POOLED PHENOTYPES' Q-Q PLOT
        if (args$qqplot_pooled) {
          logging::loginfo("Gathering all the associations ({length(pvals_pooled)}) into a single QQ-plot..." %>% glue)
          plot_qqplot(pvals_pooled, ofile=glue(qqplot_all_fp))
        }
    }
}

build_file_patterns <- function(gwas_folder, output_folder, gwas_filepattern="GWAS_{phenotype}") {

    # text files
    gwas_fp <- file.path(gwas_folder, paste0(gwas_filepattern, ".tsv"))
    # gwas_fp_rds <- file.path(gwas_folder, paste0(args$gwas_pattern, ".rds"))
    gwas_summary_fp <- file.path(output_folder, "summaries", paste0(args$gwas_pattern, "__regionwise_summary.tsv"))
    
    # figures
    figs_dir <- file.path(output_folder, "figures")
    
    if (args$figure_type == "png") { manhattan_fp <- file.path(figs_dir, paste0(args$gwas_pattern, "__manhattan.png")) } 
    else if (args$figure_type == "eps") { manhattan_fp <- file.path(figs_dir, paste0(args$gwas_pattern, "__manhattan.eps")) }
    
    qqplot_fp <- file.path(figs_dir, paste0(args$gwas_pattern, "__QQ-plot.png"))
    qqplot_all_fp <- file.path(figs_dir, "GWAS__all__QQ-plot.png")
    
    fp <- list()
    fp$manhattan_fp <- manhattan_fp
    fp$qqplot_fp <-qqplot_fp
    fp$qqplot_all_fp <- qqplot_all_fp
    fp$gwas_summary_fp <- gwas_summary_fp
    fp$gwas_fp <- gwas_fp
    fp$figs_dir <- figs_dir
    fp
}

# Process LD-independent genomic regions
load_ld_indep_regions <- function(ld_indep_regions_file="data/ld_indep_regions/fourier_ls-all_EUR_hg19.bed") {
  
    regions <- read.delim(ld_indep_regions_file, stringsAsFactors = F)
    regions <- regions %>% group_by(chr) %>% mutate(id = paste0(row_number()))
    regions$chr <-  sub("\\s+$", "", regions$chr)
    regions <- regions %>% mutate(id=paste(chr, id, sep = "_")) %>% ungroup()  
    regions
    
}

regionwise_summary <- function(gwas_df, regions, ofile) {
  
    gwas_list <- list()
    
    for (chr_number in 1:22) {
        reduced_regions <- regions %>% filter(chr==paste0("chr", chr_number))
        reduced_gwas <- gwas_df %>% filter(CHR==chr_number)
        reduced_gwas <- reduced_gwas %>% mutate(region=cut(BP, breaks=reduced_regions$start, labels=head(reduced_regions$id, -1)))
        gwas_list <- c(gwas_list, list(reduced_gwas))
    }
    
    best_p_per_region <- bind_rows(gwas_list) %>% filter(!is.na(region)) %>% group_by(region) %>% slice(which.min(P)) %>% ungroup()
    dir.create(dirname(ofile))
    write.csv(best_p_per_region, file = ofile, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

####################################################################################################

load_gwas_df <- function(file, overwrite_rds=FALSE, cache_rds=TRUE) {
  
    gwas_f <- file
    gwas_f_rds <- gsub("tsv", "rds", gwas_f)
    
    if (!file.exists(gwas_f) && !file.exists(gwas_f_rds) || (file.exists(gwas_f) && file.info(gwas_f)$size == 0)) {
        logging::logdebug(glue("File {file} does not exist or it is empty."))
        return(NULL)
    }
    
    rds_log_flag <- TRUE
    if (file.exists(gwas_f_rds) && !overwrite_rds) {
        if (rds_log_flag){
          logging::logwarn(glue("Found a previously cached RDS file: {gwas_f_rds})."))
          logging::logwarn("If you want the old file replaced, please run this script again with the --overwrite_rds flag")
          rds_log_flag <- FALSE
        }
        gwas_df <- readRDS(gwas_f_rds)
        logging::loginfo("RDS file read successfully.")
    } else {
        logging::loginfo("Reading the GWAS file in text format: {gwas_f}." %>% glue)
        gwas_df <- read_tsv(gwas_f, col_names = TRUE)
        if (cache_rds) {
          logging::loginfo(glue("Caching GWAS in RDS format."))
          saveRDS(gwas_df, gwas_f_rds)
        }
    }
    gwas_df %>% filter(!is.na(P) & P != 0)
}

plot_manhattan <- function(gwas_df, ofile, figure_type, colors, no_axis_labels=FALSE) {
  
  # ymax is closest multiple of 5 
  ymax <- 5 * (as.integer(max(-log10(gwas_df$P)) / 5) + 1)
  
  if (figure_type == "eps") {
    postscript(file=ofile, width = 12, height=4)
  } else if (figure_type == "png") {
    png(ofile, res=100, width = 3000, height = min(80 * ymax, 1200))
  }
  
  # gwas_df <- rbind(gwas_df %>% filter(P < 1e-3), gwas_df %>% filter(P < 1e-2 & P > 3e-3) %>% sample_frac(0.1), gwas_df %>% filter(P < 3e-2 & P > 1e-3) %>% sample_frac(0.5))
  # gwas_df
  
  if (no_axis_labels) {
  
      pp <- qqman::manhattan(
        gwas_df, 
        # main=plot_title, 
        cex = 0.75, cex.axis = 2, xlab="", ylab="",
        genomewideline = -log10(1.5e-10),
        suggestiveline = -log10(5e-8), 
        col = colors,
        ylim = c(0, ymax)
      )
      
  } else {
    
      pp <- qqman::manhattan(
        gwas_df, 
        # main=plot_title, 
        chrlabs = 1:22,
        cex = 0.75, cex.axis = 2,
        genomewideline = -log10(1.5e-10),
        suggestiveline = -log10(5e-8), 
        col = colors,
        ylim = c(0, ymax)
      )
      
  }
      
  dev.off()
}

plot_qqplot <- function(pvals, ofile) {
  
    if (length(pvals) == 0)
      logging::logerror("List of p-values is empty. Skipping Q-Q plot.")
    
    png(glue(ofile), res=100, width = 1000, height = 1000)
    
    pp <- qqman::qq(pvals, 
        # main=plot_title, 
        cex.axis=2, 
        col = "blue4"
    )
    dev.off()
}

########################################################################

process_phenotype <- function(phenotype, file_patterns, plot_params, which_plot) {
  
    logging::loginfo("Processing {phenotype}..." %>% glue)

    logging::loginfo("Loading GWAS summary statistics...")
    gwas_df <- load_gwas_df(file=glue(file_patterns$gwas_fp), overwrite_rds = args$overwrite_rds, cache_rds = args$cache_rds)
    if (is.null(gwas_df)) { return(NULL) }
    
    if (which_plot$manhattan) {
      man_f <- glue(file_patterns$manhattan_fp)
      logging::loginfo(glue("Creating Manhattan plot at {man_f}"))
      plot_manhattan(gwas_df, ofile=man_f, figure_type = args$figure_type, colors=plot_params$colors, no_axis_labels=plot_params$no_axis_labels)
    }
    
    if (which_plot$qqplot) {
      qq_f <- glue(file_patterns$qqplot_fp)
      logging::loginfo(glue("Creating Q-Q plot at {qq_f}"))
      plot_qqplot(gwas_df$P, ofile=qq_f)
    }
    
    if (which_plot$summary) {
      # Produce region-wise summary (best p-value per region)
      regions <- load_ld_indep_regions()
      regionwise_summary(gwas_df, regions, ofile=glue(fp$gwas_summary_fp))
    }
  
}

########################################################################

args <- fetch_args()
main(args)