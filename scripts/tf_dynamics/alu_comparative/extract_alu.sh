#!/bin/bash


REPEATMASKER_FILE=/dors/capra_lab/data/transposable_elements/repeatmasker/hg19.fa.out

TE_ARRAY=( AluYa5 AluYb8 AluSp AluY AluSc AluSg AluSq AluSx AluJb AluJo )

 


for Variable in ${TE_ARRAY[@]}
do  
    echo "$Variable"
    awk -v OFS="\t" -v TE=$Variable '{if($10 == TE){print $5,$6,$7,$10,$11,$9}}' $REPEATMASKER_FILE > ${Variable}_hg19.out
done