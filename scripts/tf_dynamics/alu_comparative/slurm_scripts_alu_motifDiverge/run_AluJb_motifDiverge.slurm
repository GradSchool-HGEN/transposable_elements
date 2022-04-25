#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=20:00:00 #10:00:00
#SBATCH --mem=30G #50G

#SBATCH --array=244-274%10 #275

#SBATCH --job-name=AluJb_motifDiverge
#SBATCH --output=/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/alu_comparative/slurm.out/AluJb/AluJb_%A_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL


#module load GCC OpenMPI R
module load R-bundle-Bioconductor
module load GSL/2.1
module load R/3.4.3-X11-20160819

#  AluYa5, AluYb8, AluSp, AluY, AluSc, AluSg, AluSq, AluSx, AluJb, AluJo,
TE_NAME="AluJb"
homolog_header="rheMac2"

tf_dir="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/TFpspm_to_run"
myfile=$(ls ${tf_dir} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}' ) 
TF_NAME=${myfile%_HUMAN.H10MO.*.pspm}


hg19_elements=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/fasta/hg19_${TE_NAME}.fa
homolog_elements=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/fasta/${homolog_header}_${TE_NAME}.fa


output_dir=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/output_motifDiverge
output_name=hg19_${homolog_header}_${TE_NAME}_${TF_NAME}

echo  ${SLURM_ARRAY_TASK_DI}: $TF_NAME "about to start for AluJb at " $(date) >> /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/alu_comparative/slurm.out/log.txt

Rscript --no-save --no-restore /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/comparative/parallel_analyze_TF_diverge.R $hg19_elements $homolog_elements ${tf_dir}/${myfile} $output_dir $output_name 

echo ${SLURM_ARRAY_TASK_DI}: $TF_NAME "completed for AluJb at " $(date) >> /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/alu_comparative/slurm.out/log.txt
