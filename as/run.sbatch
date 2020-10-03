#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --partition=commons
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --export=ALL
#SBATCH --output=log.txt
#SBATCH --gres=gpu:volta:1

cd $SLURM_SUBMIT_DIR

python -u run_md.py PDBNAME output.dcd

