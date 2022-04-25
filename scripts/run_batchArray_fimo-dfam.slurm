#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --array=1-4
#SBATCH --time=1-01:00:00  
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=14G


#SBATCH --job-name=run-leftOver-fimo-dfamT_batchArray
#SBATCH --output=leftOver-fimo-dfamTE_batchArray.out 
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL

curDir=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/split_dfamfasta/toRun/
outputdirPATH=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/split_dfamfasta
markovFile=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy_markovBackground
motifFile=/dors/capra_lab/data/TF_motifs/meme_motif_databases.12.15/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme

cd $curDir
myfile=$(ls ${curDir} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}') 
echo $myfile

fimo -verbosity 1 --text -bgfile $markovFile $motifFile $myfile > fimo-output_${myfile}

