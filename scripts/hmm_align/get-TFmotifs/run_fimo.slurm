#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=0-23:00:00
#SBATCH --mem=20G


#SBATCH --job-name=fimo-MER20_cons
#SBATCH --output=fimo-MER20_dfam_hg38.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align/get-TFmotifs

./get-fimo.sh