#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=0-5:00:00
#SBATCH --mem=16G


#SBATCH --job-name=store_path
#SBATCH --output=store_path.slurm.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


module load GCC/5.4.0-2.26 
setpkgs -a anaconda3

source activate verBio

python classify_enhancer_overlap.py

echo "DONE"