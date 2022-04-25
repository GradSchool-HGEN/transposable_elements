#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem=6G

#SBATCH --job-name=map_dict
#SBATCH --output=map_dict.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


setpkgs -a anaconda3

source activate verBio
./run_alignHMM.sh
