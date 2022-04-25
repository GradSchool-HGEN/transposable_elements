#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=3-00:00:00
#SBATCH --mem=40G


#SBATCH --job-name=fimo_repeatmasker
#SBATCH --output=fimo_repeatmasker.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


module load GCC/5.4.0-2.26 BEDTools/2.26.0
# setpkgs -a anaconda3

# source activate verBio

# python TE_annotate.py -h

# bedtools --help
# python TE_annotate.py "MER20" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100.bedGraph" -o "/dors/capra_lab/abraha1/projects/transposable_elements/results/"



cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/regulatory_potential

python run_fimo_repeatmasker.py