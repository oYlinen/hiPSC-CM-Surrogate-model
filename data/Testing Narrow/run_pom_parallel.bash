#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --array=101-200
#SBATCH --partition=small
#SBATCH --mem=10000
#SBATCH --cpus-per-task=1



cd /scratch/svc_td_cbig/olli/pom/run/matlab_10k_cntr_testing_range2/
module load miniforge3/24.9.0
module load MATLAB/2024b
source activate olli

VALUES=({0..999})

srun python ../../POMtool.py --config=matlab_config.yaml --patch_count=200 --patch_idx=${VALUES[SLURM_ARRAY_TASK_ID]} --seed=54 > /scratch/svc_td_cbig/olli/pom/run/log/java_${SLURM_ARRAY_TASK_ID}.log 2>&1
