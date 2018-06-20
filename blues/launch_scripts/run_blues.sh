#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node 1
#SBATCH --gres=gpu:p100:1
#SBATCH -t 24:00:00

# Echo commands
set -x

# Load stuff
module load anaconda5/5.0.0-3.6
module load cuda
source activate blues-nathan

# Change directory
cd $SCRATCH/blues/blues

# Run
python example_rotmove.py
