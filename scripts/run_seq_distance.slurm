#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=6G

#SBATCH --job-name=calc_LVD
#SBATCH --output=seq_distance_SLURM.out 
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


setpkgs -a anaconda3

source activate verBio

python seq_distance.py
