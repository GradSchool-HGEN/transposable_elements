#!/bin/bash 

#SBATCH --nodes=1 # Nodes 
#SBATCH --ntasks=8
#SBATCH --mem=12G

# Max job duration
#SBATCH --time=-12:00:00 # $d-hh:mm:ss, default: 00:15:00
#SBATCH --job-name=mult_Align
#SBATCH --output=mult_align
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL

cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/alu_comparative

./extract_multiz46way.sh

