#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=0-04:00:00
#SBATCH --mem=10G


#SBATCH --job-name=MER20_fantEnhPro
#SBATCH --output=MER20_fantEnhProt.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=FAIL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


module load GCC/5.4.0-2.26 BEDTools/2.26.0

cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align

time python run_TE_annotate_hmm.py 

