#!/bin/sh
#SBATCH --job-name array
#SBATCH --partition commons
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH --output=log%a.txt
#SBATCH --array=61-100

python several_mutate.py $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID


