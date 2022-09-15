#!/bin/bash

GOODRUNS="../CardiacCOMA/runs/good_runs_with_gwas_column.csv"
GOODRUNS="../CardiacCOMA/runs/good_runs.csv"
FILTER="-w False"
FILTER=""

#for RUNID in `grep $FILTER $GOODRUNS | cut -d, -f2`; do 
for RUNID in `cat ../CardiacCOMA/runs/good_runs.csv | cut -d, -f2`; do 
  python main_bgenie.py --run_id $RUNID --experiment_id 1 --gwas_software bgenie --steps_to_run 1 --dry-run;
done
