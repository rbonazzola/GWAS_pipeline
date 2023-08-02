# from IPython import embed
import shlex
from subprocess import call, check_output
import os, sys
os.chdir(check_output(shlex.split("git rev-parse --show-toplevel")).strip())

import logging

logging.basicConfig(
   level=logging.INFO,
   format="%(asctime)s [%(levelname)s] %(message)s",
   handlers=[
      logging.StreamHandler()
   ]
)

logger = logging.getLogger()
logger.info(f"Switched to folder: {os.getcwd()}")
sys.path.append(".")
from data.code.PATHS import *

import pandas as pd
from utils.ARC_helpers.SGE_utils import *


def get_regions_df():

    df = pd.read_csv(REGIONS_FILE, sep = "\t")
    df.columns = [x.strip() for x in df.columns]
    return df


def generate_incl_snps_file(full_df, region, replace_previous=False):

    '''
    Reads a file of SNP summary statistics for a full chromosome (full_df),
    and produces a file for a specific region in that chromosome.
    '''
    
    # df.columns = [x.strip() for x in df.columns]
    ofilename = REDUCED_SNP_FILE_PATTERN.format(chromosome=chromosome, start_pos=region["start_pos"], end_pos=region["end_pos"])

    if not os.path.exists(ofilename) or replace_previous:
        df_filtered = full_df[(region["start_pos"] < full_df.iloc[:,1]) & (region["end_pos"] > full_df.iloc[:,1])]
        included_snps_df = df_filtered.iloc[:,[0]]
        included_snps_df.to_csv(ofilename, index=False, header=False)
        
    return ofilename



def get_region_data(row):

    i, region = row 
    chromosome = region.chr[3:].strip()  # chrN ---> N

    chromosome_adjusted = "0" + chromosome if len(chromosome) == 1 else chromosome

    region = {
        "chromosome": chromosome_adjusted,
        "chromosome_adjusted": chromosome_adjusted,
        "start_pos": region.start,
        "end_pos": region.stop
    }

    return region


def build_bgenix_command(input_genotype_file, index_file, output_genotype_file, incl_range):

    command = ["bgenix", "-g", input_genotype_file, "-i", index_file, "-incl-range", incl_range, ">", output_genotype_file]
    command = " ".join(command)
    return command
      
 
def build_qctool_command(input_genotype_file, output_genotype_file, incl_rsid_files=None, excl_rsid_files=None, incl_range=None, incl_samples=None, samples=None):

    command = ["qctool", "-g", input_genotype_file, "-og", output_genotype_file]

    if samples is not None:
        command += ["-s", samples]

    if incl_rsid_files is not None:
        command += ["-incl-rsids"] + incl_rsid_files

    if excl_rsid_files is not None:
        command += ["-excl-rsids"] + excl_rsid_files

    if incl_range is not None:
        command += ["-incl-range", incl_range]

    if incl_samples is not None:
        command += ["-incl-samples", incl_samples]

    command = " ".join(command)
    return command


def get_QCed_snps(bgen_file, snps_stats_file=None, snps_list_file=None, bgen_ofile=None):

  if snps_stats_file is None:
      snps_stats_file = bgen_file[:-5] + "_snps_stats.txt"
      snps_stats_file = str(snps_stats_file)

  if snps_list_file is None:
      snps_list_file = bgen_file[:-5] + "_snps_passing_QC.txt"

  if bgen_ofile is None:
      bgen_ofile = bgen_file[:-5] + "_snps_passing_QC.bgen"

  qctool_stats_command = f"qctool -g {bgen_file} -snp-stats -osnp {snps_stats_file}"
  check_call(shlex.split(qctool_stats_command))
 
  snps_stats_df = pd.read_csv(snps_stats_file, sep=" ", comment="#")
  snps_stats_df = snps_stats_df[["rsid", "HW_exact_p_value", "impute_info", "missing_proportion", "minor_allele_frequency"]]                                                                                                                   
  
  qc_condition = (snps_stats_df.HW_exact_p_value > 1e-5) & (snps_stats_df.impute_info > 0.3) & (snps_stats_df.minor_allele_frequency > 0.01) 
  rsids = snps_stats_df[qc_condition].reset_index().rsid
  rsids.to_csv(snps_list_file, index=False)

  qctool_filter_bgen = build_qctool_command(input_genotype_file=bgen_file, output_genotype_file=bgen_ofile, incl_rsid_files=[snps_list_file]) 
  check_call(shlex.split(qctool_filter_bgen))



if __name__ == "__main__":

    import argparse

    df = get_regions_df()

    parser = argparse.ArgumentParser(description="Calculate areas of different geometric shapes.")

    # Additional options
    parser.add_argument(
        '--test', 
        action='store_true',
        help="Run a test"
    )

    # Adding an argument to control the verbosity level
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help="Enable verbose logging (display INFO and DEBUG messages)."
    )

    args = parser.parse_args()
    
    # Update logging level based on the --verbose argument
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    chromosome = 1

    if args.test:
        logging.info("Running a test on 5 regions")
        df_iter = df.sample(n=5).reset_index(drop=True).iterrows()
    else:
        logging.info(f"Submitting jobs for the {df.shape[0]} regions.")
        df_iter = df.iterrows()

    for i, region in df_iter:

        region = get_region_data(region)
        
        new_chromosome = region["chromosome"].lstrip("0")

        if chromosome != new_chromosome:

            try:
                logger.info(f"Created {count} batch job scripts for chromosome {chromosome}... ")
            except:
                pass
            count = 0
            chromosome = new_chromosome
            # print(f"Processing chromosome {chromosome}")
            input_genotype_file = INPUT_GENOTYPE_FILE_PATTERN.format(chromosome=chromosome)
            index_file = INDEX_FILE_PATTERN.format(chromosome=chromosome)

            if not os.path.exists(input_genotype_file):
                logger.warning(f"File {input_genotype_file} does not exist. Skipping...")
                continue
            incl_rsid_file = INCLUDE_SNP_FILES[0].format(chromosome=chromosome)

            full_df = pd.read_csv(incl_rsid_file, sep=" ", header=None)

        reduced_incl_snp_file = generate_incl_snps_file(full_df, region)
        incl_range = "{chromosome_adjusted}:{start_pos}-{end_pos}".format(**region)
        temp_genotype_file = TEMP_GENOTYPE_FILE_PATTERN.format(**region)
        output_genotype_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region)

        bgen_file = OUTPUT_GENOTYPE_FILE_PATTERN.format(**region) 
        snps_stats_file = SNPS_STATS_FILE_PATTERN.format(**region)
        snps_list_file = SNPS_LIST_FILE_PATTERN.format(**region)
        bgen_ofile = FINAL_BGEN_FILE_PATTERN.format(**region)


        path_extension = "PATH=/home/home01/scrb/bin:$PATH"
        
        bgenix_command = build_bgenix_command(
            input_genotype_file, 
            index_file, 
            temp_genotype_file, 
            incl_range
        )

        qctool_command = build_qctool_command(
            temp_genotype_file, 
            output_genotype_file, 
            incl_rsid_files=[reduced_incl_snp_file], 
            incl_range=None, 
            incl_samples=INCLUDE_SAMPLES, 
            samples=SAMPLES_FILE
        )        
        
        filter_bgen_py_command = f"python data/code/filter_bgens_for_qc_criteria.py --intermediate_bgen_file {bgen_file} --snps_stats_file {snps_stats_file} --snps_list_file {snps_list_file} --final_bgen_file {bgen_ofile}"

        commands = [
            path_extension, 
            "module unload openmp; module unload intel; module load gnu",
            # bgenix_command, 
            # qctool_command, 
            # f"rm {temp_genotype_file}",
            "module load anaconda", "source activate base",
            filter_bgen_py_command
        ]

        count += 1

        jobname = "chr{chromosome_adjusted}_{start_pos}-{end_pos}".format(**region)
        logs_file_prefix = f"{os.getcwd()}/runs/logs/preproc_bgen/{jobname}"
  
        submit_sge_job(
            jobname,
            commands=commands,
            memory_limit="8G", 
            walltime="01:00:00",
            stdout=f"{logs_file_prefix}.out",
            stderr=f"{logs_file_prefix}.err",
            dry_run=False
        )
        
        if i == 1:
            logger.info(f"Logging stdout and stderr to {logs_file_prefix}.out/err")
            logger.debug("\n".join(commands))


    logging.info(f"Created {count} batch job scripts for chromosome {chromosome}... ")
