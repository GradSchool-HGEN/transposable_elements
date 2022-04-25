#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-03-12 21:59:42


import pandas as pd 
import subprocess
import os 
import sys 
import re

FILES_TO_SPLIT_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/motifDiverge/formatted_motifDiverge_output"
ENHANCER_FILE_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_enhancer/genome_wide_predicted_placental_enhancers.bed"

match_exp = re.compile('_[A-Z0-9p]*_')  

FILES_TO_SPLIT = ["outputMER20_FOXO1_deplete_motifDiverge_qvalue.tsv",
"outputMER20_CTCFL_rerundeplete_motifDiverge_qvalue.tsv",     
"outputMER20_FOXO1_enrichment_motifDiverge_qvalue.tsv",
"outputMER20_CTCFL_rerunenrichment_motifDiverge_qvalue.tsv",  
"outputMER20_p53_rerundeplete_motifDiverge_qvalue.tsv",
"outputMER20_ETS1_deplete_motifDiverge_qvalue.tsv",            
"outputMER20_p53_rerunenrichment_motifDiverge_qvalue.tsv",
"outputMER20_ETS1_enrichment_motifDiverge_qvalue.tsv",         
"outputMER20_ZSCA4_deplete_motifDiverge_qvalue.tsv",
"outputMER20_ETS2_deplete_motifDiverge_qvalue.tsv",            
"outputMER20_ZSCA4_enrichment_motifDiverge_qvalue.tsv",
"outputMER20_ETS2_enrichment_motifDiverge_qvalue.tsv"]

# =============  FUNCTION =============


# =============  MAIN =============
for one_file in FILES_TO_SPLIT: 
    this_file = os.path.join(FILES_TO_SPLIT_PATH, one_file)

    TF_name = match_exp.search(this_file).group()[1:-1]


    if re.search('enrichment',this_file) !=None: 
        enrich_deplete_name = 'enrichment'
    else: 
        enrich_deplete_name='deplete'
    df = pd.read_csv(this_file,sep="\t",header=0) 
    df.columns = ['p_value','hg19','mm9', 'q_value']
    df['hg19chr'], df['coords'] = df['hg19'].str.split(':').str
    df['hg19_start'], df['hg19_end'] = df['coords'].str.split('-').str
    df.drop(['hg19','coords'], axis=1, inplace=True) 
    new_col_order = ['hg19chr','hg19_start','hg19_end', 'p_value','q_value','mm9']
    df = df[new_col_order]

    output_file_name = os.path.join(FILES_TO_SPLIT_PATH, "split_{}".format(one_file)) 
    df.to_csv(output_file_name, sep="\t", header=False, index=False) 

    cmd_hg19 = "bedtools intersect -loj -a {} -b {}".format(output_file_name,ENHANCER_FILE_PATH).split()
    store_lines = subprocess.check_output(cmd_hg19,shell=False, universal_newlines=True)
    
    
    output_file_name=os.path.join(FILES_TO_SPLIT_PATH, "intersect_MER20_{}_{}_predicted_placenta_enhancers.tsv".format(TF_name, enrich_deplete_name))
    with open(output_file_name, 'w') as fh: 
        fh.writelines(store_lines)

    mod_df = pd.read_csv(output_file_name, sep="\t", header=None)
    mod_df['tf_name'] = TF_name
    mod_df['enrich_status'] = enrich_deplete_name

    mod_df.to_csv(output_file_name, sep="\t", header=None, index=False)