#!/bin/bash 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=10:00:00 #10:00:00
#SBATCH --mem=50G #50G

#SBATCH --array=1-274 #275

#SBATCH --job-name=AluSg_motifDiverge
#SBATCH --output=/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/alu_comparative/slurm.out/AluSg_%A_%a.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL


#module load GCC OpenMPI R
module load R-bundle-Bioconductor
module load GSL/2.1
module load R/3.4.3-X11-20160819

#  AluYa5, AluYb8, AluSp, AluY, AluSc, AluSg, AluSq, AluSx, AluJb, AluJo,
TE_NAME="AluSg"
homolog_header="rheMac2"

tf_dir="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/TFpspm_to_run"
myfile=$(ls ${tf_dir} | awk -v line=${SLURM_ARRAY_TASK_ID} '{if (NR==line) print$0}' ) 
TF_NAME=${myfile%_HUMAN.H10MO.*.pspm}


hg19_elements=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/fasta/hg19_${TE_NAME}.fa
homolog_elements=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/fasta/${homolog_header}_${TE_NAME}.fa


output_dir=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/output_motifDiverge
output_name=hg19_${homolog_header}_${TE_NAME}_${TF_NAME}

Rscript --no-save --no-restore /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/comparative/parallel_analyze_TF_diverge.R $hg19_elements $homolog_elements ${tf_dir}/${myfile} $output_dir $output_name 
