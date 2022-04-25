#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-16 19:21:47


import os 
import time
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

# mainPath ="/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/facet_expressed_enhancers/alltxt/allfacets/"
ROOT_PATH  ="/Users/abin-personal/Google Drive/GD_transfer/data/"
OUTPUT_PATH ="/Users/abin-personal/Google Drive/GD_transfer/output/"
FILE_NAME = "summary_permuted_TF_byTissue_pValues.tsv"


# =============  main =============

df=pd.read_table(os.path.join(OUTPUT_PATH,FILE_NAME), sep="\t", header=0)

print(df)

depletedTF = df[(df['TE_overlap_ratio_pValue'] > 0.95) & (df['count'] > 10)]
len(depletedTF)
enrichedTF = df[(df['TE_overlap_ratio_pValue'] < 0.05) & (df['count'] > 10)]
len(enrichedTF)

depletedTF = df[(df['TE_overlap_ratio_pValue'] > 0.95) & (df['count'] > 5)]
len(depletedTF)
enrichedTF = df[(df['TE_overlap_ratio_pValue'] < 0.05) & (df['count'] > 5)]
len(enrichedTF)