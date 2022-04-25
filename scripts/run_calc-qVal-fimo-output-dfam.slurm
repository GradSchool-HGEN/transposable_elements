#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=0-23:00:00
#SBATCH --mem=100G


#SBATCH --job-name=calc_qVal_fimo-dfam
#SBATCH --output=calc_qval_fimo-dfam-output
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


module load GCC/5.4.0-2.26 BEDTools/2.26.0
# setpkgs -a anaconda3

# source activate verBio

# python TE_annotate.py -h

# bedtools --help
# python TE_annotate.py "MER20" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100.bedGraph" -o "/dors/capra_lab/abraha1/projects/transposable_elements/results/"

cd /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/fimo-dfam

qvalue --column 7 --append fimo-output_hg38_dfam.nrph.hits_tidy > q-val_fimo-output_hg38_dfam.nrph.hits_tidy 
