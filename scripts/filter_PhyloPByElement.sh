#!/bin/sh
# This script will filter phyloP100 data that has already been intersected with dfam database based on a given TE. 
#
#
#
# Abin Abraham
# created on: 2017-12-28 22:27:02



interphyloP_dir=""
elementOfInterest="MER20"

cd $interphyloP_dir

awk -v OFS=‘\t’ '{ if($4==$elementOfInterest){print $1, $2, $3,$4, $10}}' phyloP_Dfam_chr* >> phyloP_MER20