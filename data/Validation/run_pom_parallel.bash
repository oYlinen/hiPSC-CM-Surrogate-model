#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --array=0-100
#SBATCH --partition=small
#SBATCH --mem=10000
#SBATCH --cpus-per-task=1



cd /scratch/svc_td_cbig/olli/pom/run/matlab_5k_cntr/
module load miniforge3/24.9.0
module load MATLAB/2024b
source activate olli

VALUES=({0..999})

srun python ../../POMtool.py --config=matlab_config.yaml --patch_count=100 --patch_idx=${VALUES[SLURM_ARRAY_TASK_ID]} --seed=52 > /scratch/svc_td_cbig/olli/pom/run/log/java_${SLURM_ARRAY_TASK_ID}.log 2>&1
