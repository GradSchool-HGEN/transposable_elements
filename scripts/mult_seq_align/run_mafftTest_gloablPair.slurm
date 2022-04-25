#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=20:00:00
#SBATCH --mem=80G

#SBATCH --job-name=mafft_test_v4_gloablpair
#SBATCH --output=maff_test_v4_globalpair.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


inputfile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta"

setpkgs -a anaconda3
source activate verBio

cd /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/multAlign/mafft

echo $inputfile
start=$SECONDS
echo "mafft, test using default, input file is full fasta, and 8 threads"

mafft --globalpair --thread 8 $inputfile > mafft-out_MER20_gloablpair
#mafft --thread 8 $inputfile > mafft_default.

duration=$(( SECONDS - start ))
echo $duration
