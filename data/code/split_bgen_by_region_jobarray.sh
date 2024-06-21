#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N split_bgen_by_region
#$ -t 1-1703
#$ -o logs/split_bgen/output_$TASK_ID.log
#$ -e logs/split_bgen/error_$TASK_ID.log

# Print the hostname of the machine the job is running on
echo "Running on $(hostname)"

# Print the task ID
echo "Task ID: $SGE_TASK_ID"

# Define the input file based on the task ID
INPUT_FILE="input_file_${SGE_TASK_ID}.txt"

# Run your command with the input file
python  $INPUT_FILE