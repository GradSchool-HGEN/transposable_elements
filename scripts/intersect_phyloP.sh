#!/bin/bash

# This script will do bedtools intersect on dfam TE and phylo P scores. 
# Abin Abraham 


hg38_dfamfile=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy
phyloPdirPATH=/dors/capra_lab/abraha1/data/hg38_phastCon100way/hg38.phyloP100way/hg38.phyloP100way_byChr/
outputdir=/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phyloP-dfamTE/new/

cd $phyloPdirPATH

for currChr in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}
#for currChr in {1,2}
do
    #grep -w "chr"$currChr $hg38_dfamfile
    #grep -w "chr"$currChr $hg38_dfamfile | bedtools intersect -loj -wb -a 'stdin' -b ${phyloPdirPATH}"chr"$currChr".bed" -sorted > "phyloP_Dfam_chr"$currChr"" 
    grep -w "chr"$currChr $hg38_dfamfile | bedtools intersect -wb  -a ${phyloPdirPATH}"chr"$currChr".bed" -b 'stdin' -sorted > $outputdir"phyloP_Dfam_chr"$currChr"" 

    
done



