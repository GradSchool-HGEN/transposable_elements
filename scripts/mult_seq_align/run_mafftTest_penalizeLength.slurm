#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --mem=12G

#SBATCH --job-name=mafft_test_v3_penalizeLenth
#SBATCH --output=maff_test_v3_penalizeLength.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta"

setpkgs -a anaconda3
source activate verBio

cd /dors/capra_lab/abraha1/projects/transposable_elements/results/mult_align

echo $inputfile
start=$SECONDS
echo "mafft, peanlizing length, input file is full fasta, and 8 threads"

mafft --thread 8 --parttree --ep 0.123 $inputfile > mafft_default.
#mafft --thread 8 $inputfile > mafft_default.

duration=$(( SECONDS - start ))
echo $duration
