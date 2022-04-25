#!/bin/bash

TE_ARRAY=(AluJb AluJo AluSc AluSg AluSp AluSq AluSx AluYa5 AluYb8 AluY) 
FIMO_REPEATMASKER_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output/all_TE_fimo_out.txt"



for i in "${TE_ARRAY[@]}"; do
    echo "$i"

    # grep -w $i $FIMO_REPEATMASKER_FILE >> ${i}_store_output.txt 

done



awk -v TE=$i '{if($2=="AluJb" || $2=="AluJo" || $2=="AluSc" || $2=="AluSg" || $2=="AluSp" || $2=="AluSq" || $2=="AluSx" || $2=="AluYa5" || $2=="AluYb8" || $2=="AluY") print}' $FIMO_REPEATMASKER_FILE > matchedAlu_fimo_repeatmasker.txt