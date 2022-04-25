#!/bin/sh

#this script will run FIMO (meme-suite) on dfam hits of transposable elements
#created 2017-12-15 16:34:15
#Abin Abraham
date

#cd to the directory for the batch fimo to be run on:
cd /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/split_dfamfasta/

#SET UP VARIABLES 
outputdir=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/split_dfamfasta/results/
markovFile=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_markovBackground
motifFile=/dors/capra_lab/data/TF_motifs/meme_motif_databases.12.15/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme

for i in *.fa; do 
echo "processing file .. "

fimo -verbosity 1 --o $outputdir -bgfile $markovFile $motifFile $i

done

