#!/bin/python
# This script will filter out rows where there is a Alu element overlap with a FANTOM ENHANCER.  
#
#
#
# Abin Abraham
# created on: 2018-04-27 10:06:46


# ELEMENTS TO FILTER OUT 
    # - AluYa5
    # - AluYb8
    # - AluSp
    # - AluY
    # - AluSc
    # - AluSg
    # - AluSq
    # - AluSx
    # - AluJb
    # - AluJo



import pandas as pd 
import numpy as np 
import os, sys

#=========== FILE PATHS ===========
ENHANCER_TE_ROOT_PATH = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/fantom_enhancers"
ENHANCER_TE_FILE_NAME = "all_fimo_out_withTFCoord_diff_expressed_enhancers_TIDY_TE_intersect.tsv"

OUTPUT_DIR = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/tf_dynamics/alu_TF_cluster"

#=========== LOAD DATA & TIDY UP ===========
df = pd.read_csv(os.path.join(ENHANCER_TE_ROOT_PATH, ENHANCER_TE_FILE_NAME), sep="\t", header=None)
col_names = ['TF_genome_chr', 'TF_genome_start', 'TF_genome_end', 'TFname', 'sequence_name', 'tissue', 'score', 'p_value', 'q_value', 'matched_sequence', 'TFquality', 'chr_TE', 'start_TE', 'end_TE', 'TE', 'TE_family']
df.columns = col_names
df['enhancer_chr'], df['enhancer_start_end'] = df['sequence_name'].str.split(':', 1).str
df['enhancer_start'], df['enhancer_end'] = df['enhancer_start_end'].str.split('-', 1).str

new_cols = ['enhancer_chr','enhancer_start', 'enhancer_end', 'tissue', 'TF_genome_chr', 'TF_genome_start', 'TF_genome_end', 'TFname', 'score', 'p_value', 'q_value', 'matched_sequence', 'TFquality', 'chr_TE', 'start_TE', 'end_TE', 'TE', 'TE_family']
new_df = df.loc[:,new_cols]

#=========== FILTER & WRITE FILE ===========
elements_to_filer = ["AluYa5", "AluYb8", "AluSp", "AluY", "AluSc", "AluSg", "AluSq", "AluSx", "AluJb", "AluJo"]

for this_element in elements_to_filer: 

    df_to_write = new_df.loc[new_df.iloc[:,16] == this_element, :] 
    file_name = os.path.join(OUTPUT_DIR, "{}_FANTOM_enh_repeatmasker_TE_TF_intersection.tsv".format(this_element))
    df_to_write.to_csv(file_name, sep="\t", header=True, index=False)

print("completed filtering!")




