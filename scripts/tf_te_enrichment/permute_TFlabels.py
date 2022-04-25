#!/bin/python
# This script will permute TF labels across all enhancers and create a ratio of TE overlap to non-TE overlap for each TF. 
#
#
#
# Abin Abraham
# created on: 2018-02-13 07:28:44

import os 
import time 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import normaltest


# mainPath ="/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/facet_expressed_enhancers/alltxt/allfacets/"
ROOT_PATH ="/Users/Abin/Google Drive/GD_transfer/data/"
OUTPUT_PATH="/Users/Abin/Google Drive/GD_transfer/output/"
FILE_NAME ="all-fimo-out-withTFCoord-filtered-TE-intersect.tsv"

# =============  functions =============

def plot_hist_per_TF(perm_values):
    plt.hist(perm_values)
    plt.show()


# =============  main =============
rdf = pd.read_table(os.path.join(ROOT_PATH,FILE_NAME), sep="\t", header=None)
rdf.columns = ["chr(TFmotif)", "start(TFmotif)"," end(TFmotif)","TF", "enhancer_sequence", "q-value", "strand", "motif_quality", "motif-fimo-score", "tissue", "chr(TE)","start(TE)","end(TE)" ,"TE", "TEfam"]

### subset 
df = rdf.loc[:,['TF','enhancer_sequence','TE',]]
df['TE_overlap'] = np.where(df['TE']!=".", True, False)

### counts of each TF by TE overlap 
tf_te_counts = df.groupby(['TF','TE_overlap']).size().unstack()
tf_te_n = dict(zip(tf_te_counts.index, np.sum(tf_te_counts.fillna(0).values,axis =1).astype('int')))

### dictionaries to store permutation values; index 0 is the original sample value 
tf_te_ratios = dict(zip(tf_te_counts.index.values, tf_te_counts.loc[:,True]/tf_te_counts.loc[:,False]))
tf_te_count = dict(zip(tf_te_counts.index.values, tf_te_counts.loc[:,True]))

### shuffle and calculate ratio of TE overlap to non-TE overlap
shuffle_df = df.copy()
NUM_ITERS_PERM = 10000
num_iters = NUM_ITERS_PERM
now = time.time()
while num_iters > 0: 
    print(num_iters)
    permuted_index = np.random.permutation(np.arange(df.shape[0]))
    shuffle_df.iloc[:,0] =  shuffle_df.iloc[permuted_index, 0].values # shuffle TF labels 
    shuffle_df_counts = shuffle_df.groupby(['TF','TE_overlap']).size().unstack()
    shuffle_df_counts.fillna(0)
    for this_row in np.arange(shuffle_df_counts.shape[0]):
        this_tf = shuffle_df_counts.index.values[this_row]
        ratio = (shuffle_df_counts.iloc[this_row,1]+1)/(shuffle_df_counts.iloc[this_row, 0]+1) 

        
        perm_true_hits = shuffle_df_counts.iloc[this_row,1] # ratio of true/false w/ TE overlap 
        tf_te_ratios[this_tf] = np.append(tf_te_ratios[this_tf], ratio)
        tf_te_count[this_tf] = np.append(tf_te_count[this_tf], perm_true_hits)
        # print(tf_te_ratios[this_tf])

    num_iters = num_iters -1 
        
print("finished! {}".format((time.time()-now)/60))

### open dictionary and calculate percentiles
store_labels = np.empty(0)
store_percentiles = np.empty(0)
store_emp_p_value = np.empty(0)
store_n = np.empty(0)
test_counter = 1
for k,v in tf_te_ratios.items(): 
    store_labels = np.append(store_labels,k)

    emp_p_value = len(v[v>=v[0]])/len(v)
    store_percentiles = np.append(store_percentiles, len(v[v<v[0]])/len(v))
    store_n = np.append(store_n, tf_te_n[k])
    store_emp_p_value = np.append(store_emp_p_value, emp_p_value)
    
    # print("for {}, ratios are {}".format(k,v[~np.isnan(v)]))
    print("for {}, ratios are {}".format(k,v))

    plt.clf()

    plt.hist(v[~np.isnan(v)])
    plt.axvline(x=v[0], color='r' )
    plt.title("Hsitogram of {} n=({}), perm = ".format(k, tf_te_n[k], NUM_ITERS_PERM))
    plt.xlabel('Ratio of number of TFs with and without enhancer overlap'.format(k))
    # plt.show()
    # plt.savefig(os.path.join(OUTPUT_PATH,"TF_Figs/ratio/{}.png".format(k)))
    plt.clf()

    # counter + 1 
    # if test_counter == 13 :break


### open dictionary for hits of enhancers overlap and calc percentiles
store_percentiles_counts = np.empty(0)
store_n_counts = np.empty(0)
for k,v in tf_te_count.items(): 
    store_labels_count = np.append(store_labels,k)

    store_percentiles_counts = np.append(store_percentiles_counts,len(v[v<v[0]])/len(v))
    store_n_counts = np.append(store_n, tf_te_n[k])
    print("for {}, counts are {}".format(k,v[~np.isnan(v)]))

    plt.clf()
    plt.hist(v[~np.isnan(v)])
    plt.axvline(x=v[0], color='r' )
    plt.title("Histogram of {} n=({}), perm={}".format(k, tf_te_n[k], NUM_ITERS_PERM))
    plt.xlabel('Number of TFs from {} that Overlap Enhancers'.format(k))
    # plt.show()
    # plt.savefig(os.path.join(OUTPUT_PATH,"TF_Figs/enhc_overlap_counts/{}_counts.png".format(k)))
    plt.clf()

### create a df with percentile and total_num_of_TFs
rounded_store_percentiles = [ round(elem, 4) for elem in store_percentiles.tolist() ]
rounded_store_percentiles_counts = [ round(elem, 4) for elem in store_percentiles_counts.tolist() ]
tf_percentile_df = pd.DataFrame( {"emp_P_value_ratio":store_emp_p_value,"percentile_from_ratio" : rounded_store_percentiles,"percentile_from_counts" : rounded_store_percentiles_counts, "total_num_of_TFs" :store_n},
            index = store_labels)

tf_percentile_df.sort_values('total_num_of_TFs').to_csv(os.path.join(OUTPUT_PATH, 'tf_percentile_from_perm_test_wPvalue.tsv'),sep="\t")





# plt.scatter(tf_percentile_df.total_num_of_TFs, tf_percentile_df.percentile,s=1 )
# plt.show()

### plot histogram/boxplots of count of TFs 
# counts_tf = np.sum(tf_te_counts.fillna(0).values,axis =1).astype('int')
# q75, q25 = np.percentile(counts_tf, [75 ,25])
# iqr = q75 - q25
 
# minIQR = q25 - (iqr*1.5)
# maxIQR = q75 + (iqr*1.5)

# plt.boxplot(counts_tf)
# plt.title('distribution of count of each TFs in Enhancer Datasets, with outliers')
# plt.ylim(-1,maxIQR)
# plt.show()

# #plot 
# width = 1
# ind = np.arange(len(store_labels))
# plt.bar(ind,store_percentiles, width)
# plt.title('Title')
# plt.ylabel('Percentile Score')
# plt.xlabel('TFs')
# plt.xticks(ind, store_labels, rotation='vertical')
# plt.tight_layout()
# plt.show()




# plt.hist(store_values, 50, density=True, facecolor='g', alpha=0.75)
# plt.xlabel('Smarts')
# plt.ylabel('Probability')
# plt.title('Histogram of IQ')
# plt.show()