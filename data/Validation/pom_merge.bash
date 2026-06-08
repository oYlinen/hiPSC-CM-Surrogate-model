#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=05:00:00
#SBATCH --partition=small
#SBATCH --mem=10000
#SBATCH --cpus-per-task=1



module load miniforge3/24.9.0
module load MATLAB/2024b
source activate olli

srun python ../../POMtool.py merge --config=matlab_config.yaml --patch_count=100

