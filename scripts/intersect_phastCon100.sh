#!/bin/bash

# This script will do bedtools intersect on dfam TE and phastCon100 scores. 
# Abin Abraham 


hg38_dfamfile=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy
phastConFilePATH=/dors/capra_lab/abraha1/data/hg38_phastCon100way/phastCon100.bedGraph
outputdir=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100-dfamTE/

cd /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon100-dfam

bedtools intersect -wb  -a $phastConFilePATH -b $hg38_dfamfile -sorted > phastCon100-dfam.out 



