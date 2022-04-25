#!/bin/python
# This script will creates a tidy data set. 
#   final columns are: ['enhc_chr', 'enhc_start', 'enhc_end', 'tissue', 'TFname','TFstart', 'TFstop', 'strand', 'score', 'p_value', 'q_value', 'matched_sequence', 'TFquality' ]
#   see FIMO_FILE for input file 
#   output is saved to OUTPUT_DIR
# 
# Abin Abraham
# created on: 2018-02-14 14:21:50


### note: fimo output is 1 based start inclusive; this script will convert to genomic coordinates with BED12 like coordinates: 0 based start, end point not included. 

import os 
import pandas as pd 
import numpy as np 


FIMO_FILE = "/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/differentially_expressed/alltxt/all_fimo_out_diff_expressed.txt"
OUPUT_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fantom_enhancers"


#=========== main ===========
df = pd.read_csv(FIMO_FILE, sep='\t', header=None, skiprows=1)
df.columns =[ 'pattern_name', 'sequence_name', 'TFstart', 'TFstop', 'strand', 'score', 'p_value', 'q_value', 'matched_sequence', 'tissue']
df = df[df['TFstart'] !='start'] #multiple headers are spread throughout the file; those are removed

df['TFname'] = df['pattern_name'].str.split('.',0, expand=True).iloc[:, [0]]
df['TFquality'] = df['pattern_name'].str.split('.',0, expand=True).iloc[:, [2]]
df.drop(['pattern_name'], axis=1, inplace=True)

df['enhc_chr'] = df['sequence_name'].str.split(':',0, expand=True).iloc[:, [0]]
df['enhc_coord'] = df['sequence_name'].str.split(':',0, expand=True).iloc[:, [1]]
# df.drop(['sequence_name'], axis=1, inplace=True)

df['enhc_start'] = df['enhc_coord'].str.split('-',0, expand=True).iloc[:, [0]]
df['enhc_end'] = df['enhc_coord'].str.split('-',0, expand=True).iloc[:, [1]]
df.drop(['enhc_coord'], axis=1, inplace=True)


df[['enhc_start','enhc_end','TFstart', 'TFstop']] = df[['enhc_start','enhc_end','TFstart', 'TFstop']].apply(pd.to_numeric)
df['TF_genome_chr'] = df['enhc_chr']
df['TF_genome_start'] = df['enhc_start'] + df['TFstart'] - 1 
df['TF_genome_end'] = df['enhc_start'] + df['TFstop'] - 1 +1 # 0 based, non-inclusive end 


df = df.loc[:,['TF_genome_chr', 'TF_genome_start', 'TF_genome_end', 'TFname', 'sequence_name', 'tissue', 'score', 'p_value', 'q_value', 'matched_sequence', 'TFquality' ]]
# df = df.loc[:,['enhc_chr', 'enhc_start', 'enhc_end', 'tissue', 'TFname','TFstart', 'TFstop', 'strand', 'score', 'p_value', 'q_value', 'matched_sequence', 'TFquality' ]]

df['q_value'] = df['q_value'].apply(pd.to_numeric)
df.head()
df = df[((df['TFquality'] == "A") | (df['TFquality'] == "B") | (df['TFquality'] == "C")) & (df['q_value'] < 0.10)]

df.to_csv(os.path.join(OUPUT_DIR,'all_fimo_out_withTFCoord_diff_expressed_enhancers_TIDY.tsv'),sep="\t", header=None, index=False )