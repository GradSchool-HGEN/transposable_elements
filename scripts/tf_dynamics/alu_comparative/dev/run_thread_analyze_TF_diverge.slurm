#!/bin/bash

#SBATCH --mail-user=abraham.abin13@gmail.com  
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=09:30:00
#SBATCH --mem=50G
#SBATCH --output=%j_rerun_CTCF_div.out

#module load GCC OpenMPI R
module load R-bundle-Bioconductor
module load GSL/2.1
module load R/3.4.3-X11-20160819

R --version


TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/CTCF_HUMAN.H10MO.A.pspm"
# TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/ETS1_HUMAN.H10MO.C.pspm"
# TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/FOXO1_HUMAN.H10MO.C.pspm"
# TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/CTCFL_HUMAN.H10MO.A.pspm"
# TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/ETS2_HUMAN.H10MO.C.pspm"
# TF_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/HOCOMOCOv10_HUMAN_mono_meme_tsv/TF_to_run/P53_HUMAN.H10MO.B.pspm"


Rscript --no-save --no-restore parallel_analyze_TF_diverge.R /dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/hg19_MER20.fa /dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/mm9_MER20.fa ${TF_PATH} /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/comparative MER20_CTCF_rerun
