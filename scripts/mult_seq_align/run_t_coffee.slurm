#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --mem=12G

#SBATCH --job-name=t_coffee_v1_all
#SBATCH --output=t_coffee_v1_all.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta"

setpkgs -a anaconda3
source activate verBio

cd /dors/capra_lab/abraha1/projects/transposable_elements/results/mult_align/t_coffee/


echo $inputfile
start=$SECONDS
echo "t_coffee, test using default, input file is full fasta, and 8 threads"

t_coffee $inputfile -type=dna 
duration=$(( SECONDS - start ))
echo $duration

#consider these settting 
# t_coffee <inputfile.fasta> -type=dna 
