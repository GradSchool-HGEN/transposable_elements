#!/bin/bash
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=heatmap_test_%j.out

# Load the Anaconda distribution of Python, which comes
module load Anaconda3

# Pass your Python script to the Anaconda3 python intepreter for execution
python create_full_heatmap.py
