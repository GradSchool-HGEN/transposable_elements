#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-03-07 08:54:39

import pandas as pd 
import numpy as np 
import os, sys 


INTERSECT_ENHANCER_FILE ="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_enhancer/hg19_MER20_intersect_enhancer_wOverlap.tsv"
CEBP_ENRICH_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/motifDiverge/enrich_CEBPb_MER20_hg19_mm9.out"

# load enrich/depelete 
# load enhancer intersect 
# calc odds ratio

# -----------
# MAIN
# ----------- 

### load data 
df = pd.read_csv(INTERSECT_ENHANCER_FILE,sep="\t",header=0)
df.columns = ['TE_chr', 'TE_start', 'TE_end', 'TE_strand', 'Villar_chr',
        'Villar_start', 'Villar_end', 'enhancer_activity', 'overlap_with_TE',
       'overlap_with_EnhancerActivity']

df_enrich = pd.read_csv(CEBP_ENRICH_FILE, sep="\t", header=0)
df_enrich.columns = ['q_value', 'p_value', 'hg19_TE_coords', 'mm9_TE_coords']

df['hg19_TE_coords'] = df["TE_chr"]+ ":" + df["TE_start"].map(str) + "-" + df["TE_end"].map(str)

#merge enrich w/ intersect enhancer 
merged_df = pd.merge(df, df_enrich, how='outer', on='hg19_TE_coords')

### ## TO DO: change p-values to q-values 
merged_df['noEnrich_noEnhc'] = np.where((df['q_value'] >=0.05) & (df['overlap_with_EnhancerActivity']<0.5), 1, 0)
merged_df['noEnrich_yesEnhc'] = np.where((df['q_value'] >=0.05) & (df['overlap_with_EnhancerActivity']>0.5), 1, 0)
merged_df['yesEnrich_noEnhc'] = np.where((df['q_value'] <=0.05) & (df['overlap_with_EnhancerActivity']<0.5), 1, 0)
merged_df['yesEnrich_yesEnhc'] = np.where((df['q_value'] <=0.05) & (df['overlap_with_EnhancerActivity']>0.5), 1, 0)
## need 4 counts 
## 