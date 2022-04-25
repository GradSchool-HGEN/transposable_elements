
#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-11 10:11:18

import os
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import binom_test

ROOT_PATH = "/Users/Abin/Google Drive/GD_transfer/data"
INTERSECT_FILE = "intersect_loj_TFBS_facet_enhancers.tsv"


# =============  functions ========


# =============  main =============


df = pd.read_csv(os.path.join(ROOT_PATH, INTERSECT_FILE), header=None, sep="\t")
df.columns = ["TF_chr", "TF_start", "TF_end", "TF", "TE", "TE_locus", "motif_id", "strand", "pValue", "motifSeq", "enhc_chr", "enhc_start", "enhc_end", "enhc_overlap"]
df['enhc_overlap'] = np.where(df['enhc_chr']!=".", True, False)


# create array with rows = TFs, columns = TEs only if there is an enhancer overlap 
# create another array, same as above but without enhancer overlap 
enh_df = df.loc[df['enhc_overlap']==True, ['TF','TE']]
enh_df['count'] = 1
enh_df = enh_df.groupby(['TF','TE']).sum().unstack().fillna(0)['count']

no_enh_df = df.loc[df['enhc_overlap']==False, ['TF','TE']]
no_enh_df['count'] = 1
no_enh_df = no_enh_df.groupby(['TF','TE']).sum().unstack().fillna(0)['count']

# binomial test for each TE, for a given TF, calculate number of TF w/ and w/o enhancer overlap
bi_df = df.groupby(['TF','TE', 'enhc_overlap']).count()['TF_chr'].unstack().fillna(0)

# 
plt.scatter(bi_df.values[:,0],bi_df.values[:,1], s=1)    
plt.xlabel('Number of TF sites on TE without Enhancer Overlap')
plt.ylabel('Number of TF sites on TE with Enhancer Overlap')
plt.title('Count of TF sites on TEs by Enhancer Overlap\n (Each dot represents one TF on one TE)')
plt.show()


# calcualte proportion of TF from enhancers vs. not enhancers 
gb_df = df.groupby(['TF','enhc_overlap']).count()['TF_chr'].unstack()
gb_df['ratio'] = gb_df.fillna(0).iloc[:,1]/gb_df.fillna(0).iloc[:,0]
gb_df['total_n'] = gb_df.fillna(0).iloc[:,1]+gb_df.fillna(0).iloc[:,0]
gb_df= gb_df.sort_values('ratio', ascending=False)

byTF_Overlap_ratio = gb_df.iloc[:,2].values
TF_labels = [gb_df.index.values[x]+" ("+gb_df['total_n'].values.astype('int').astype('str')[x]+ ")" for x in np.arange(len(gb_df.index.values))]

width = 1
ind  = np.arange(len(byTF_Overlap_ratio))

fig, ax = plt.subplots(figsize=(14,5))
plt.bar(ind, byTF_Overlap_ratio )
plt.title('Ratio of TF found on Enhancers vs. Not Enhancers Collapsed Across TEs')
plt.ylabel('Count of TF found on Enhancers vs. Not Enhancers')
plt.xlabel('TFs from LTRs')
plt.xticks(ind, TF_labels, rotation='vertical')
plt.tight_layout()
plt.show()


# for each TF, claculate ratio of enhancer overlap with no enhancer overlap 
gb_tf = df.groupby(['TF','enhc_overlap']).count()['TF_chr']

ratio = bi_df.values[:,0]/bi_df.values[:,0]