#!/bin/bash
#SBATCH --mail-user=abraha1@vanderbilt.edu  
#SBATCH --mail-type=ALL
#SBATCH --nodes=6
#SBATCH --tasks-per-node=1
#SBATCH --time=012:30:00
#SBATCH --mem=10G
#SBATCH --output=full_benchmark.out

module load GCC OpenMPI R
R --version

srun --mpi=pmi2 Rscript analyze_TF_diverge.R
