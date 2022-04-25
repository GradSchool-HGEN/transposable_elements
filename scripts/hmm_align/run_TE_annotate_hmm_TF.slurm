#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=0-04:00:00
#SBATCH --mem=10G


#SBATCH --job-name=THE1B_anno_TF
#SBATCH --output=THE1B_anno_TF.out
#SBATCH --mail-user=abraham.abin13@gmail.com
#SBATCH --mail-type=ALL #BEGIN or END or FAIL or ALL(begin,end,fail) or TIME_LIMIT_X, where X is percent of allocated


module load GCC/5.4.0-2.26 BEDTools/2.26.0
# setpkgs -a anaconda3

# source activate verBio

# python TE_annotate.py -h

# bedtools --help
# python TE_annotate.py "MER20" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100.bedGraph" -o "/dors/capra_lab/abraha1/projects/transposable_elements/results/"

cd /dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align

time python run_TE_annotate_hmm.py TFmotif


# annotation = ["phastCon100", "phyloP100", "TFmotif", "fantEnh","fantEnh_Promoters"]


