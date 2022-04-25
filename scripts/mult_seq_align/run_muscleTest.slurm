#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=10:00:00
#SBATCH --mem=6G

#SBATCH --job-name=faster_muscle_miters3
#SBATCH --output=faster_muscle_miters3.out     
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


cd /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam

input=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_MER20.fasta

echo "muscle, faster option using maxiter3"
muscle -in hg38_dfam.nrph.hits_tidy_MER20.fasta -out muscle-output_MER20_maxiter3 -maxiters 3 
