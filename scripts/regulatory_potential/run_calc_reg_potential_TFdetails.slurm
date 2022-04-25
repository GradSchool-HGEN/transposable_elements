#!/bin/bash


#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=5:00:00#$bd-hh:mm:ss, default: 00:15:00
#SBATCH --mem=5G

#SBATCH --array=1177,46,1006,1007,1013,1017,1042,1047,1048,1050,1052,1074,1081,1083,10,11,121,124,125,13,15,167,16,170,173,174,17,18,19,1,20,213,21,233,23,249,250,255,266,268,269,275,276,277,278,279,281,282,283,284,285,286,287,288,289,28,290,291,292,293,294,295,296,297,298,299,2,300,304,305,307,309,310,311,312,313,314,316,320,321,322,323,328,329,330,333,338,341,349,350,352,353,355,356,357,358,359,360,361,362,363,364,365,366,368,369,371,374,381,382,383,384,385,386,388,3,423,43,443,44,47,4,517,6,739,748,7,825,829,830,831,8,926,927,928,929,931,933,935,937,939,943,971,974,980,982,984,9%10

#SBATCH --job-name=rerun_byTF_regPotential
#SBATCH --output=rerun_byTF_regPotential.out
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


OUTPUT_DIR=/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/regulatory_potential/byTF/
INPUT_DIR=/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output/individualTE_fimo_output/ 
myfile=$(ls ${INPUT_DIR} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}') 
OUTPUT_FILE=${myfile/.txt/_regPotential_byTF.tsv}

cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/regulatory_potential
python calc_reg_potential_TFdetails.py -f $INPUT_DIR$myfile -o $OUTPUT_DIR$OUTPUT_FILE

#total 1177 jobs 