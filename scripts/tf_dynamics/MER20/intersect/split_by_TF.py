#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-26 22:00:14


import pandas as pd 

df_file = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/intersect/inter-TE-anno.tsv"

col_name = ["intersect_chr", "intersect_start", "intersect_end", "TF", "TE_genomic_chr","TE_genomic_start","TE_genomic_end", "TE_model","model_start","model_end", "model_length"]
df = pd.read_csv(df_file, sep="\t", header=None)
df.columns  = col_name

uniq_TF = df.iloc[:,3].unique() 

for this_TF in uniq_TF: 
    tf_df = df.loc[df['TF']==this_TF,:]
    tf_df['TF'] = 1 
    tf_df.to_csv("/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/intersect/inter-TE-anno_by{}.tsv".format(this_TF),
    sep="\t", index=False, header=False)
