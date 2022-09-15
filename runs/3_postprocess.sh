#!/bin/bash

# for RUNID in `grep -w False ../CardiacCOMA/runs/good_runs_with_gwas_column.csv | cut -d, -f2`; do 

for RUNID in `cut -d, -f2 ../CardiacCOMA/runs/good_runs.csv`; do
  python main_bgenie.py --run_id $RUNID --experiment_id 1 --gwas_software bgenie --steps_to_run 3; 
done
