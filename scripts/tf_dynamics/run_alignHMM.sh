#!/bin/bash 


input_file=/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/common_LTRs.txt

while read LINE; 
    do echo $LINE

    python /dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align/alignHMM.py $LINE | tee $LINE_alignHMM_output.txt
done < $input_file
