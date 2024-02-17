import shlex
import re
import os
from pprint import pprint
import subprocess
import pandas as pd
import yaml
import logging

logging.basicConfig(level=logging.WARNING,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')

from easydict import EasyDict

HARMONIZATION = False

##############################################################################################################

GWAS_COLUMNS = {
  "--snp_column": "SNP",
  "--effect_allele_column": "a_1",
  "--non_effect_allele_column": "a_0",
  "--beta_column": "BETA",
  "--pvalue_column": "P"
}     

GWAS_COLUMNS = {
  "--snp_column": "panel_variant_id",
  "--effect_allele_column": "effect_allele",
  "--non_effect_allele_column": "non_effect_allele",
  "--beta_column": "effect_size",
  "--pvalue_column": "pvalue"
}     


PARSING_SCRIPT = "/home/home01/scrb/01_repos/GWAS_pipeline/utils/summary-gwas-imputation/src/gwas_parsing.py"
SPRED_EXEC = "/home/home01/scrb/01_repos/GWAS_pipeline/utils/MetaXcan/software/SPrediXcan.py"
VARIANT_METADATA = "/home/home01/scrb/01_repos/GWAS_pipeline/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz"

##############################################################################################################
################################### S-PrediXcan commands

def harmonization_command(gwas_file, gwas_harmon_file):
    
    #### HARMONIZATION

    REF_DATA_ARGS = [
        "-snp_reference_metadata", "%s METADATA" % VARIANT_METADATA, 
        "-liftover", config['CHAIN_FILE']
    ]
    
    HARMON_COLMAPPING_OPTIONS = [
        '-output_column_map', 'SNP', 'variant_id', 
        '-output_column_map', 'a_1', 'effect_allele',
        '-output_column_map', 'a_0', 'non_effect_allele',
        '-output_column_map', 'BP', 'position',
        '-output_column_map', 'CHR', 'chromosome', '--chromosome_format', 
        '-output_column_map', 'af', 'frequency', 
        '-output_column_map', 'P', 'pvalue',
        '-output_column_map', 'BETA', 'effect_size',
        '-output_column_map', 'SE', 'standard_error'    
    ]
    
    OUTPUT_ORDER = [ 
        '-output_order', 
        'variant_id', 
        'panel_variant_id', 
        'chromosome', 
        'position', 
        'effect_allele', 
        'non_effect_allele', 
        'frequency', 
        'pvalue', 
        'zscore', 
        'effect_size', 
        'standard_error' 
    ]
                
    # harmon_commands = {}
    harm_command = [ 'python', PARSING_SCRIPT ]
    harm_command.extend([ '-gwas_file', gwas_file, '-output', gwas_harmon_file ])                
    harm_command.extend(REF_DATA_ARGS)
    harm_command.extend(HARMON_COLMAPPING_OPTIONS)
    harm_command.extend(OUTPUT_ORDER)    
    
    return " ".join(harm_command)
    
    
def spredixcan_command(gwas_file, spredixcan_ofile, model_db, covariance_file):
    
    def to_command(t):
        t = list(t.items())
        t = [item for sublist in t for item in sublist]
        spredixcan_command = [SPRED_EXEC] + t
        spredixcan_command = " ".join(spredixcan_command)
        return spredixcan_command

    args = {
      "--model_db_path": model_db,
      "--covariance": covariance_file,
      "--gwas_folder": os.path.dirname(gwas_file),
      "--gwas_file_pattern": os.path.basename(gwas_file),
      "--output_file": spredixcan_ofile,
      "--model_db_snp_key": "varID",
      "--keep_non_rsid": ""
      # "--output_file": spredixcan_output_pattern.format(run, run, z, tissue),
    }
    
    # args = {
    #    "--model_db_path": models[tissue],
    #    "--covariance": covariances[tissue],
    #    "--gwas_folder": os.path.dirname(gwas_file),
    #    '--gwas_file_pattern': os.path.basename(gwas_file),
    #    "--output_file": spredixcan_ofile,
    #    # "--output_file": spredixcan_output_pattern.format(run, run, z, tissue),
    # }
    
    args.update(GWAS_COLUMNS)

    return to_command(args)
    
    # with open("spred_commands.sh", "w") as script:
    #     script.write("#!/bin/bash\n")
    #     
    #     for tissue in arguments:
    #         for phenotype in arguments[tissue]:
    #             spredixcan_command = to_command(arguments[tissue][phenotype])
    #             script.write(spredixcan_command + "\n")
    
####################################################################################
            
def main(config, args):
    
    # gwas_harmon_pattern = "GWAS__{}__harmonized__test_std_covariates_PC__GBR_unrelated__qc.tsv"
  
    import re
    gwas_file_regex = re.compile(config.gwas_file_regex)
    gwas_harmon_file_regex = re.compile(config.gwas_harmon_regex)
    phenotype = "1223dbdba9314a0e8427a519ff524ef9_z000"
    tissue = "Heart_Left_Ventricle"
    
    config.gwas_folder = "/home/home01/scrb/01_repos/GWAS_pipeline/output/coma/finetuned_mlruns" # /GWAS_1223dbdba9314a0e8427a519ff524ef9_z000.tsv
    
    # gwas_file = lambda phenotype: config.gwas_pattern.format(phenotype=phenotype)
    gwas_file = args.gwas_file
    gwas_harmon_file = args.gwas_harmon_file
    # phenotype = gwas_file_regex.match(gwas_file).group(1)
    phenotype = gwas_harmon_file_regex.match(gwas_harmon_file).group(1)
    print(phenotype)
    logging.info("Processing %s", phenotype)
    gwas_harmon_file = lambda phenotype: config.gwas_harmon_pattern.format(phenotype=phenotype)
    
    # gwas_ifile = os.path.join(config.gwas_folder, gwas_file)
    gwas_ofile = os.path.join(config.gwas_folder, gwas_harmon_file(phenotype))
    
    # harm_command = harmonization_command(gwas_ifile, gwas_ofile)
    
    with open(args.job_file, "wt") as ff:
        ff.write("module load anaconda; source activate imlabtools\n\n")
        if HARMONIZATION == True:
            ff.write(harm_command)
        ff.write("\n")
    # subprocess.call(harm_command)
    
    ####################################################################################################
  
    for tissue in args.tissues:

        model_db = f"{config.model_db_folder}/mashr_{tissue}.db"
        covariance_file = f"{config.model_db_folder}/mashr_{tissue}.txt.gz"
    
        # spredixcan_ofolder = ""
        spredixcan_ofile = config.spredixcan_output_pattern.format(phenotype, tissue)
    
        spredixcan_output_pattern = os.path.join(config.gwas_folder, config.spredixcan_output_pattern)
    
        spred_command = spredixcan_command(gwas_ofile, spredixcan_ofile, model_db, covariance_file)

        with open(args.job_file, "at") as ff:
            ff.write(spred_command)
            ff.write("\n")

        # print(spred_command)
        # subprocess.call([spred_command])

    logging.info("Written job to %s", args.job_file)
    

            
if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(description="Create batch jobs for harmonization and S-PrediXcan.")
    
    parser.add_argument('--config', "-c", type=str, default="config_postprocessing_gtex_v8_mashr.yaml")
    parser.add_argument('--gwas_folder', type=str)
    parser.add_argument('--gwas_file', type=str)
    parser.add_argument('--gwas_harmon_file', type=str)
    parser.add_argument('--job_file', type=str)
    parser.add_argument('--tissues', type=str, nargs="+")

    
    # Adding an argument to control the verbosity level
    parser.add_argument(
        '--verbose', action='store_true',
        help="Enable verbose logging (display INFO and DEBUG messages)."
    )
    
    args = parser.parse_args()
    
    # Update logging level based on the --verbose argument
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)
    
    ##################################################
    
    with open(args.config) as config_f:
        
        if not os.path.exists(args.config):
            raise ValueError(f"File {args.config} does not exist.")
            
        config = yaml.safe_load(config_f)
        config = EasyDict(config)
        
    ##################################################
        
    # latent_variables = ["z" + str(i) for i in range(8)]
    # phenotypes = [(run, z) for run in config['runs'] for z in latent_variables]    
    
    main(config, args)
    
