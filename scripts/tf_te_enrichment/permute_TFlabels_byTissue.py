#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-13 07:28:44

import os 
import time
import pandas as pd 
import numpy as np
import pickle 
import matplotlib.pyplot as plt
from collections import defaultdict


# mainPath ="/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/facet_expressed_enhancers/alltxt/allfacets/"
ROOT_PATH ="/Users/abin-personal/Google Drive/GD_transfer/data/"
OUTPUT_PATH="/Users/abin-personal/Google Drive/GD_transfer/output/"
FILE_NAME ="all_fimo_out_withTFCoord_diff_expressed_enhancers_TIDY_TE_intersect.tsv"


# =============  main =============

### load data 
rdf=pd.read_table(os.path.join(ROOT_PATH,FILE_NAME), sep="\t", header=None)
rdf.columns = ["chr(TFmotif)", "start(TFmotif)"," end(TFmotif)","TF", "enhancer_sequence", "tissue", "motif_fimo_score", "p_value", "q_value","match_sequence",  "motif_quality", "chr(TE)","start(TE)","end(TE)" ,"TE", "TEfam"]

### subset 
df = rdf.loc[:,['TF','enhancer_sequence','TE','TEfam','tissue']]
df['TE_overlap'] = np.where(df['TE']!=".", True, False)

store_overlap_count_dict = defaultdict(dict)
store_ratio_dict = defaultdict(dict)
store_count = defaultdict(dict) 

### observed counts of each TF by TE overlap 
for this_tissue_actual in df['tissue'].unique(): 

    by_tissue_df_actual = df[df['tissue']==this_tissue_actual]
    by_tissue_df_actual_count = by_tissue_df_actual.groupby(['TF','TE_overlap']).size().unstack().fillna(0) 

    for this_row_tf_actual in np.arange(by_tissue_df_actual_count.shape[0]): 

        this_tf_label_actual = by_tissue_df_actual_count.index.values[this_row_tf_actual]
        print("this tissue: {}, for this tf: {}, the overlap count is {}".format(this_tissue_actual, this_tf_label_actual, by_tissue_df_actual_count.iloc[this_row_tf_actual,1] ))

        store_overlap_count_dict[this_tf_label_actual][this_tissue_actual] = by_tissue_df_actual_count.iloc[this_row_tf_actual,1]
        store_ratio_dict[this_tf_label_actual][this_tissue_actual] = (by_tissue_df_actual_count.iloc[this_row_tf_actual,1]+1)/(by_tissue_df_actual_count.iloc[this_row_tf_actual,0]+1) #overlap/no overlap 
        
        store_count[this_tf_label_actual][this_tissue_actual] = by_tissue_df_actual_count.iloc[this_row_tf_actual,1] + by_tissue_df_actual_count.iloc[this_row_tf_actual,0]  #count total number of TF per tissue 

### shuffle and calculate ratio of TE overlap to non-TE overlap
now = time.time()
for this_tissue in df['tissue'].unique(): 
    
    by_tissue_df = df[df['tissue']==this_tissue] ## remove TF that do not exist in this tissue context 

    shuffled_by_tissue_df = by_tissue_df.copy()
    NUM_ITERS_PERM = 10
    num_iters = NUM_ITERS_PERM

    while num_iters > 0:
        print(num_iters)
        permuted_index = np.random.permutation(np.arange(by_tissue_df.shape[0])) #permute 
        shuffled_by_tissue_df.iloc[:,0] = shuffled_by_tissue_df.iloc[permuted_index,0].values #shuffle TF labels found in this tissue context 
        shuffled_by_tissue_df_count = shuffled_by_tissue_df.groupby(['TF', 'TE_overlap']).size().unstack().fillna(0)

        for this_tf_row in np.arange(shuffled_by_tissue_df_count.shape[0]):
            #seletcion logic based on number... 
            this_tf_label = shuffled_by_tissue_df_count.index.values[this_tf_row]
            this_tf_count = shuffled_by_tissue_df_count.iloc[1,1] #count of True TE overlap 
            this_tf_ratio = (shuffled_by_tissue_df_count.iloc[1,1]+1)/(shuffled_by_tissue_df_count.iloc[1,0]+1) # ratio of true/false
            store_overlap_count_dict[this_tf_label][this_tissue] = np.append(store_overlap_count_dict[this_tf_label][this_tissue] ,this_tf_count )    
            store_ratio_dict[this_tf_label][this_tissue] = np.append(store_ratio_dict[this_tf_label][this_tissue] ,this_tf_ratio )    
        
        num_iters = num_iters -1 

print("finished! {}".format((time.time()-now)/60))

### (3576 total TF-Tissue pairs) 
### calculate empircal p_value 
store_emp_p_value_overlap_count = defaultdict(dict) 
for this_tf, this_tissue_dict in store_overlap_count_dict.items():
    for tissue_label, value in this_tissue_dict.items(): 
        
        emp_p_value = len(value[value>=value[0]])/len(value)
        store_emp_p_value_overlap_count[this_tf][tissue_label] = emp_p_value
        print("for {}, in {}, the emp_Pvalue is {:.4f}".format(this_tf, tissue_label,emp_p_value))

#repeat for ratio
store_emp_p_value_ratio = defaultdict(dict) 
for this_tf, this_tissue_dict in store_ratio_dict.items():
    for tissue_label, value in this_tissue_dict.items(): 
        
        emp_p_value = len(value[value>=value[0]])/len(value)
        print("for {}, in {}, the emp_Pvalue is {:.4f}".format(this_tf, tissue_label,emp_p_value))
        store_emp_p_value_ratio[this_tf][tissue_label] = emp_p_value


### create table summarizing results 
# create a list of dictionary, where each dictionary is a row with column:value as key:value pairs 

row_list = list() 
for this_TF, tissue_dict in store_count.items():
    for this_tissue, this_count in tissue_dict.items(): 

        this_dict_row = {'TF':this_TF, 'tissue':this_tissue, 'count':this_count, 'TE_overlap_count_pValue':format(store_emp_p_value_overlap_count[this_TF][this_tissue], '.5f'), 
        'TE_overlap_ratio_pValue':format(store_emp_p_value_ratio[this_TF][this_tissue],'.5f')}
        row_list.append(this_dict_row)
        

### save data table
summary_df = pd.DataFrame(row_list)
summary_df = summary_df.loc[:, ['TF','tissue','count','TE_overlap_count_pValue','TE_overlap_ratio_pValue']]
summary_df.to_csv(os.path.join(OUTPUT_PATH,"summary_df.tsv"), sep="\t", header=True, index=False)


# ### save dictionaries 
# save_dictionary = [store_count, store_overlap_count_dict, store_ratio_dict , store_emp_p_value_overlap_count, store_emp_p_value_ratio, store_normality, store_normality_ratio]
# pickle.dump(save_dictionary, open('permute_TFlabels_byTissue_dictionaries.pi', "wb"))

# ### load dictionaries 
# save_dictionary = pickle.load( open( os.path.join(OUTPUT_PATH,"permute_TFlabels_byTissue_dictionaries.pi"), "rb" ) )
# store_count, store_overlap_count_dict, store_ratio_dict , store_emp_p_value_overlap_count, store_emp_p_value_ratio, store_normality, store_normality_ratio = save_dictionary





### plot histograms

# ### raw data characteristics 
# # Number of TF hits by tissue 
# ind = np.arange(len(df.tissue.value_counts()))
# plt.bar(ind, df.tissue.value_counts().values)
# plt.title("Raw Number of TF Hits by Tissue")
# plt.xticks( ind, df.tissue.value_counts().index.values, rotation='vertical')
# plt.tight_layout()
# plt.grid()
# plt.show()
