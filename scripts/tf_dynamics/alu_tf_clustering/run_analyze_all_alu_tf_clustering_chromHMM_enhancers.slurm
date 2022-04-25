#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=./slurm_out/tf_enrich_test_%j.out

# Load the Anaconda distribution of Python, which comes
module load Anaconda3

# Pass your Python script to the Anaconda3 python intepreter for execution
python analyze_all_alu_tf_clustering_chromHMM_enhancers.py
#####SBATCH --ntasks=8
