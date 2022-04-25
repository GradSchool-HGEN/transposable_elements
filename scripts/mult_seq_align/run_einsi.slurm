#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --mem=24G

#SBATCH --job-name=einsi_v1
#SBATCH --output=einsi_v1_full.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta"

setpkgs -a anaconda3
source activate verBio

cd /dors/capra_lab/abraha1/projects/transposable_elements/results/mult_align

echo $inputfile
start=$SECONDS
echo "mafft, test using default, input file is full fasta, and 8 threads"
mafft --thread 8 --genafpair --maxiterate 1000 $inputfile > einsi-output_MER20
duration=$(( SECONDS - start ))
echo $duration
