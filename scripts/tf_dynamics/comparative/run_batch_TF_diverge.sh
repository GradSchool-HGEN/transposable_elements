#!/bin/bash 

#SBATCH --nodes=1 # Nodes 
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20G
#SBATCH --array=1-6
#SBATCH --time=05:30:00

# Max job duration
#SBATCH --job-name=batch_TF_diverge
#SBATCH --output=%j_%a_TF_div.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL 

echo SBATCH_JOB_NAME:  $SBATCH_JOB_NAME 
echo SLURM_JOBID:  $SLURM_JOBID
echo SLURM_ARRAY_TASK_ID:  $SLURM_ARRAY_TASK_ID
echo SLURM_ARRAY_JOB_ID:  $SLURM_ARRAY_JOB_ID

TF_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run"

myfile=$(ls ${TF_DIR} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}' )

module load R-bundle-Bioconductor
module load GSL/2.1
module load R/3.4.3-X11-20160819


Rscript --no-save --no-restore parallel_analyze_TF_diverge.R /dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/hg19_MER20.fa /dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/mm9_MER20.fa ${myfile} /dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/motifDiverge MER20_${myfile%_HUMAN.H10MO.*.pspm}



#comments 

