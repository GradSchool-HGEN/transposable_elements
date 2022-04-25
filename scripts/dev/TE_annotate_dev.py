#!/bin/python
# This script will overlay a given annotation for a transposable element or a list of elements over its consensus sequence
#
# Abin Abraham      
# created on: 2017-12-07 20:15:16
# 
# Depends on: 
#       --> TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# 
# Inputs: 
#       1) transposable element of interest 
#       2) annotation file 
# 
# Output: 
#       1) 

import pandas as pd 
import subprocess
import argparse

# -----------
# Arguments 
# ----------- 

parser = argparse.ArgumentParser(description="a short description of what this does for help text ~")

parser.add_argument("-e", "--element_name",
                    dest = elem, action='store', type=str,
                    help="Select a transposable element") 
parser.add_argument("-a", "--annotation_file",
                    action='store', type=str,
                    dest = 'anno_file_path', 
                    help="Path or file name for annotation") 
parser.add_argument("-c", "--checkfile"
                    action='store_true', type=str,
                    dest='check_switch', default=False,  
                    help="checks input file format")
parser.add_argument("-s", "--skip_checkFile"
                    action='store_false', type=str,
                    dest='check_switch', default=False, 
                    help="skips checks input file format")
parser.add_argument("--sort", 
                    action='store_true', 
                    dest='sort_flag', default=False, 
                    help="sort files before analysis")                    

argsIN = parser.parse_args()

element_NAME = argsIN.elem
annotationFile_NAME = argsIN.anno_file_path
checkFile_FLAG = argsIN.check_switch
sort_FLAG = argsIN.sort_flag


# -----------
# Check Data Formats, and sort if necessary  
# ----------- 

if sort_FLAG: 
    subprocess.Popen(['sort', '-k1,1', '-k2,2n', annotationFile_NAME, '>', 'sorted_'+annotationFile_NAME ],  universal_newlines=True)

    

# Read in TE data, then filter based on TE of interest 
# TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# annotation_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

# -----------
# MAIN
# ----------- 

# idnefity an element 
# filter TEDfam_dataset on element 
# sort TEDfam_dataset sort -k1,1 -k2,2n in.bed > in.sorted.bed
# sort annotation file the same way 

# run intersection on TEDfam file and annotation file  

obs_intersect = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotation_file, '-b', TEbed_file],  universal_newlines=True)

# obs_intersect = subprocess.Popen(['bedtools', 'intersect', '-wb', '-a', TEbed_file, '-b', annotation_file], stdout=subprocess.PIPE, universal_newlines=True)
# obs_intersect2 = subprocess.check_output(['bedtools', 'intersect','-wa', '-wb', '-a', 'stdin', '-b', TEbed_file], stdin=obs_intersect.stdout, universal_newlines=True)

nl = obs_intersect.split('\n')
out_data =[x.split('\t') for x in nl]
col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "genomic_chr","genomic_start","genomic_end",  \
            "TE_model","model_start","model_end", "model_length" ] 
df = pd.DataFrame( data = out_data, columns=col_name)
df = df.apply(pd.to_numeric, errors = 'ignore')
df["mapToConsensus_start"] = df.intersect_start - df.genomic_start + df.model_start
df["mapToConsensus_end"] = df.intersect_end - 1- df.genomic_start + df.model_start


