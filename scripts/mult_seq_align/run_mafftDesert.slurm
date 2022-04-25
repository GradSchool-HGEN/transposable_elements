#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --mem=36G

#SBATCH --job-name=mafft_DesertRemoved
#SBATCH --output=mafft_DesertRemoved.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


#inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta"
inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20_nonDesert.fasta"

setpkgs -a anaconda3
source activate verBio

cd /dors/capra_lab/abraha1/projects/transposable_elements/results/mult_align

echo $inputfile
start=$SECONDS
echo "mafft, MER20 with > 20 bases in desert regions are removed, and 8 threads"
mafft --thread 8 $inputfile > mafft-output_DesertRemoved
duration=$(( SECONDS - start ))
echo $duration
