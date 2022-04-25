#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=6:30:00#$bd-hh:mm:ss, default: 00:15:00
#SBATCH --mem=2G

#SBATCH --array=382,1,11,13,16,19,21,381

#SBATCH --job-name=rerun_regPotential
#SBATCH --output=rerun_regPotential.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 

echo ----------------------------------
echo SBATCH_JOB_NAME:  $SBATCH_JOB_NAME
echo SLURM_JOBID:  $SLURM_JOBID

echo SLURM_ARRAY_TASK_ID:  $SLURM_ARRAY_TASK_ID
echo SLURM_ARRAY_JOB_ID:  $SLURM_ARRAY_JOB_ID
echo ----------------------------------

setpkgs -a anaconda3
source activate verBio


OUTPUT_DIR=/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/regulatory_potential/
INPUT_DIR=/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output/individualTE_fimo_output/ 

myfile=$(ls ${INPUT_DIR} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}') 
OUTPUT_FILE=${myfile/.txt/_regPotential.tsv}

cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/regulatory_potential
python calc_reg_potential.py -f $INPUT_DIR$myfile -o $OUTPUT_DIR$OUTPUT_FILE

