#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=2-12:00:00
#SBATCH --mem=80G


#SBATCH --job-name=v2fimo_hg38
#SBATCH --output=v2fimo_hg38.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 


HG38_FASTA_FILE="/dors/capra_lab/data/dna/hg38.fasta"
FIMO_BACKGROUND_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_hg38/markov_background_hg38" 
FIMO_MOTIFDATABASE_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/HOCOMOCOv10_HUMAN_mono_meme_format.meme_filtered"
FIMO_OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_hg38/fimo_out_hg38_v2" 



fimo --o $FIMO_OUTPUT_DIR --qv-thresh --thresh 0.05 --bgfile $FIMO_BACKGROUND_FILE $FIMO_MOTIFDATABASE_FILE $HG38_FASTA_FILE 