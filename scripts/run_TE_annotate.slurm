#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-23:00:00
#SBATCH --mem=1G


#SBATCH --job-name=TE_annotate_MER20
#SBATCH --output=TE_annotate_MER20.slurm.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


module load GCC/5.4.0-2.26 BEDTools/2.26.0
# setpkgs -a anaconda3

# source activate verBio

# python TE_annotate.py -h

# bedtools --help
# python TE_annotate.py "MER20" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100.bedGraph" -o "/dors/capra_lab/abraha1/projects/transposable_elements/results/"

python TE_annotate_analysis.py


