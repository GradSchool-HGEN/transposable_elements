#!/bin/python
# # Analyzing TF motifs Across Different Alu Elements 
#
#
# Alu Elements Analyzed in Su, M. et al.
#     - AluYa5
#     - AluYb8
#     - AluSp
#     - AluY
#     - AluSc
#     - AluSg
#     - AluSq
#     - AluSx
#     - AluJb
#     - AluJo
# 
# Reference: Su, M. et al., 2014. Evolution of Alu elements toward enhancers. CellReports, 7(2), pp.376–385.
#
# Abin Abraham
# created on: 2018-04-30 07:15:44

import time
import os, sys 
import pandas as pd 
import numpy as np 
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
plt.ioff() 

#=========== FILE PATHS ===========


ROOT_PATH_chromHMMEnhancer = "/dors/capra_lab/users/abraha1/projects/transposable_elements/data/tf_dynamics/alu_TF_cluster/chromHMM_Alu_intersection"
CHROMHMM_ENHANCER_INTERACTION_FILE={
    "AluJb":"AluJb_chromHMMenhancer_intersect.out",
    "AluJo":"AluJo_chromHMMenhancer_intersect.out",
    "AluSc":"AluSc_chromHMMenhancer_intersect.out",
    "AluSg":"AluSg_chromHMMenhancer_intersect.out",
    "AluSp":"AluSp_chromHMMenhancer_intersect.out",
    "AluSq":"AluSq_chromHMMenhancer_intersect.out",
    "AluSx":"AluSx_chromHMMenhancer_intersect.out",
    "AluYa5":"AluYa5_chromHMMenhancer_intersect.out",
    "AluYb8":"AluYb8_chromHMMenhancer_intersect.out",
    "AluY":"AluY_chromHMMenhancer_intersect.out"}

OUTPUT_FIG_ROOT="/dors/capra_lab/users/abraha1/projects/transposable_elements/results/alu_tf_clustering"
OUTPUT_FIG_DIR = "{}_figs".format(time.strftime("%d_%m_%Y_%I-%M-%S"))
os.mkdir(os.path.join(OUTPUT_FIG_ROOT, OUTPUT_FIG_DIR))
OUTPUT_FIG_PATH = os.path.join(OUTPUT_FIG_ROOT,OUTPUT_FIG_DIR)


#=========== MAIN ===========
### LOAD FILE PATHS 
TE_ANALYZED = 'AluJo'
df = pd.read_csv(os.path.join(ROOT_PATH_chromHMMEnhancer, CHROMHMM_ENHANCER_INTERACTION_FILE[TE_ANALYZED]), sep="\t", header=None, usecols=np.arange(13))

### CLEAN UP DATA
df.columns = ["TF_chr", "TF_start","TF_end", "TF","TE_coordinate",
              "motif_score", "strand", "p_value","q_value","motif_seq","enhancer_chr","enhancer_start","enhancer_end"]

df = df.drop([ "motif_seq", "p_value", "strand"], axis=1)
nodups_df = df.drop_duplicates(subset=['TF','TE_coordinate'], keep=False)

df['enhancer_overlap'] = np.where(df['enhancer_chr'] ==".", False, True);

### PRELIMINARY ANALYSIS 
## histogram of number of TF counts on each instance of TE 
raw_tf_counts = df.groupby(by='TE_coordinate').count()['TF']
unique_tf_counts = nodups_df.groupby(by='TE_coordinate').count()['TF']
TE_carrying_tf = nodups_df.groupby(by='TF').count()['TE_coordinate']

# ### plot
# fig, ax = plt.subplots(1,3)
# _ = ax[0].hist(raw_tf_counts.values, alpha = 0.5, color = 'r')
# _ = ax[0].set_xlabel('# of Raw TF Motifs')
# _ = ax[0].set_ylabel('count')
# _ = ax[0].set_title('# of TF per TE')

# _ = ax[1].hist(unique_tf_counts.values, alpha = 0.5, color = 'g')
# _ = ax[1].set_xlabel('# of Unique TF Motifs');
# _ = ax[1].set_ylabel('count');
# _ = ax[1].set_title('# of unique TF per TE')

# _ = ax[2].hist(TE_carrying_tf.values, alpha = 0.5, color = 'b')
# _ = ax[2].set_xlabel('# of TEs with A Specific TF');
# _ = ax[2].set_ylabel('count');
# _ = ax[2].set_title('# of TEs with a Specific TF')

# plt.tight_layout()
# plt.savefig(os.path.join(OUTPUT_FIG_PATH,'{}_hist_TFcounts_on_TEs__{}.eps').format(time.strftime("%d_%m_%Y"),TE_ANALYZED))
# plt.close(fig)

### MARK TEs with at least one Enhancer Overlap 
TE_no_enhancer_once =  set(df.loc[df['enhancer_overlap']== False,'TE_coordinate'].unique())
TE_yes_enhancer =  set(df.loc[df['enhancer_overlap']== True,'TE_coordinate'].unique())
TE_with_no_enhancer = TE_no_enhancer_once.difference(TE_yes_enhancer)

len(TE_yes_enhancer)
len(TE_no_enhancer_once)
len(TE_with_no_enhancer)
len(TE_no_enhancer_once.intersection(TE_yes_enhancer))


### GROUPBY
agg_func = {'motif_score':['count','mean','median']}
gb_df = nodups_df.groupby(by=['TE_coordinate','TF']).agg(agg_func).unstack(fill_value=0)
count_gb_df = gb_df['motif_score']['count']
count_gb_df.shape

### Heatmap of TF count ### 
import seaborn as sns; sns.set(color_codes=True)

### enhancer_overlap_key
cmap = dict(zip(count_gb_df.index,'r'*len(count_gb_df.index)))

### change those TE with overlap into 'r' 
for this_te in TE_with_no_enhancer: 
    cmap[this_te] = 'b'
    
# 'r' will be enhancer overlap, 'b' will be no enhancer overlap at all 

### create column color labels 
row_colors = list()
for this_index in count_gb_df.index: 
    row_colors.append(cmap[this_index])

g = sns.clustermap(count_gb_df.iloc[:,:], row_colors=row_colors)
plt.yticks(rotation=0) 
plt.title(TE_ANALYZED) 
plt.savefig(os.path.join(OUTPUT_FIG_PATH,'{}_heatmap_TFvsTE_withEnhancer__{}.eps').format(time.strftime("%d_%m_%Y"),TE_ANALYZED))
plt.close(fig)


