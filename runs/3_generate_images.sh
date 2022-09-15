RUNID=1
EXPID=2c83b43cbb314a9795911d3981360067

Rscript src/postprocessing/gwas_analysis.R \
--output_folder . \
--gwas_folder output/CardiacCOMA/${RUNID}_${EXPID} \
--gwas_pattern GWAS__${RUNID}_${EXPID}_{phenotype} \
--qqplot_pooled \
--cache_rds
