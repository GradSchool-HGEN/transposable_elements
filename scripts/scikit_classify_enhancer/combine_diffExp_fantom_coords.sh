#!/bin/sh


cd /dors/capra_lab/data/enhancers/FANTOM/differentially_expressed/tissues

for i in UBERON_00*; do 

  echo $i
  newName=${i#*UBERON_00*_}
  finalName=${newName%_differentially_expressed_enhancers.bed}
  
  awk -v OFS="\t" -v awk_name=$finalName '{print $1,$2,$3,awk_name}' $i >> all_differentially_expressed_enhancers.tsv

done 

