#!/bin/python
# This script will create a summary table for enhancer and TFBS intersection on MER20 for the list of files provided. 
# The enhancers come from predicted placental enhancers.  
# THE TFBS are the output of motifDiverge; MER20 hg19 and mm9 homologs were tested for enrichment/depletion for a given TFBS. 
# Then they were intersected with placental enhancers.
#
#
# Abin Abraham
# created on: 2018-03-13 08:23:17


import numpy as np
import pandas as pd 
import os 
import sys 

INTERSECT_FILE_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/motifDiverge/formatted_motifDiverge_output"
INTERSECT_FILE_LIST=["intersect_MER20_CTCFL_deplete_predicted_placenta_enhancers.tsv",     
                    "intersect_MER20_CTCFL_enrichment_predicted_placenta_enhancers.tsv",  
                    "intersect_MER20_ETS1_deplete_predicted_placenta_enhancers.tsv",      
                    "intersect_MER20_ETS1_enrichment_predicted_placenta_enhancers.tsv",   
                    "intersect_MER20_ETS2_deplete_predicted_placenta_enhancers.tsv",      
                    "intersect_MER20_ETS2_enrichment_predicted_placenta_enhancers.tsv",   
                    "intersect_MER20_FOXO1_deplete_predicted_placenta_enhancers.tsv",     
                    "intersect_MER20_FOXO1_enrichment_predicted_placenta_enhancers.tsv",  
                    "intersect_MER20_p53_deplete_predicted_placenta_enhancers.tsv",       
                    "intersect_MER20_p53_enrichment_predicted_placenta_enhancers.tsv",    
                    "intersect_MER20_ZSCA4_deplete_predicted_placenta_enhancers.tsv",     
                    "intersect_MER20_ZSCA4_enrichment_predicted_placenta_enhancers.tsv"]
OUTPUT_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/motifDiverge/summary_counts_MER20_placental_enhancer.tsv"


# =============  MAIN =============

store_rows = list()

### load data 

for this_file in INTERSECT_FILE_LIST: 
    this_file_full_path = os.path.join(INTERSECT_FILE_DIR, this_file)
    df = pd.read_csv(this_file_full_path, sep="\t", header=None)
    df.columns = ['hg19_chr', 'hg19_start','hg19_end','p_value','q_value','mm9', 'placental_enhancer_chr','placental_enhancer_start','placental_enhancer_end', 'TF_name','enrich_deplete']

    df['YES_TFenrich_YES_enhancer'] = np.where((df['q_value'] <=0.1) & (df['placental_enhancer_chr'] != "."), 1, 0)
    df['YES_TFenrich_NO_enhancer'] = np.where((df['q_value'] <=0.1) & (df['placental_enhancer_chr'] == "."), 1, 0)
    df['NO_TFenrich_YES_enhancer'] = np.where((df['q_value'] > 0.1) & (df['placental_enhancer_chr'] != "."), 1, 0)
    df['NO_TFenrich_NO_enhancer'] = np.where((df['q_value'] >  0.1) & (df['placental_enhancer_chr'] == "."), 1, 0)
    df['YES_TFenrich'] = np.where(df['q_value'] <=0.1, 1, 0)
    df['YES_enhancer'] = np.where(df['placental_enhancer_chr'] != ".", 1, 0)

    TF_name = df.TF_name[0]
    enrich_deplete_status = df.enrich_deplete[0]

    summary = (TF_name, enrich_deplete_status, len(df[df['YES_TFenrich_YES_enhancer'] == 1]) , len(df[df['YES_TFenrich_NO_enhancer'] == 1]), 
        len(df[df['NO_TFenrich_YES_enhancer'] == 1]), len(df[df['NO_TFenrich_NO_enhancer'] == 1]),
        len(df[df['YES_TFenrich'] == 1]), len(df[df['YES_enhancer'] == 1]), len(df))

    store_rows.append(summary)



### put data into a table 
col_labels_summary_df = ['TF', 'enrich_deplete', 'YES_TFchange_YES_enhancer', 'YES_TFchange_NO_enhancer', 'NO_TFchange_YES_enhancer', 'NO_TFchange_NO_enhancer',
   'YES_TFchange', 'YES_enhancer', 'total_count' ]

summary_df = pd.DataFrame(store_rows, columns=col_labels_summary_df)
summary_df.to_csv(OUTPUT_FILE, sep="\t", header=True, index=False)