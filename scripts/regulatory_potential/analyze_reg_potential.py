#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-01-30 17:43:43

# look at how regpotential score varies by 
    # AGE, split apart by TE Main Family Type 


import pandas as pd 
import numpy as np 


# REG_POTENTIAL_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/regulatory_potential/all_fimo_TE_regPotential.tsv"
# TE_LINEAGE_FILE = "/dors/capra_lab/data/transposable_elements/repeatmasker/hg19_TE_counts_wlin.txt"

REG_POTENTIAL_FILE = "/Users/Abin/Google Drive/transfer/data/all_fimo_TE_regPotential.tsv"
TE_LINEAGE_FILE = "/Users/Abin/Google Drive/transfer/data/hg19_TE_counts_wlin.txt"

#-------
# functions
#-------



#-------
# main
#-------

df = pd.read_csv(REG_POTENTIAL_FILE, sep='\t',header=None,names=['TE','TEcoord','RegPotential_score',"numTF_perTE"])
tf = pd.read_csv(TE_LINEAGE_FILE, sep='\t',header=0)
mf= pd.merge(df,tf, how='inner', on='TE')
mf = mf.drop('Count', axis=1)
mf.to_csv("/Users/Abin/Google Drive/transfer/data/RegPotential_wLin.tsv", sep='\t', index=False, index_label=False)
