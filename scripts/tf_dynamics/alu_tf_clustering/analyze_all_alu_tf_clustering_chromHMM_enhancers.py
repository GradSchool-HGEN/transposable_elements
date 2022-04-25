#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-05-06 15:01:21



### heatmap with enhancer overlap (subset of data)
### heatmap with TE that have enhancer overlap (subset of data)
### calculate OR for TF across all TEs 
### calculate OR for TF across TEs that overlap enhancer only. 


import os, sys
import numpy as np 
import pandas as pd

import time

import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
# =============  FILE PATHS =============
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


OUTPUT_DIR = "/dors/capra_lab/users/abraha1/projects/transposable_elements/results/alu_tf_clustering"

# =============  FUNCTIONS =============
def calc_OR(dataframe, this_tf):
    true_TF_true_enhancer = dataframe[(dataframe.TF == this_tf) & (dataframe.enhancer_overlap == True)].shape[0]
    true_TF_false_enhancer = dataframe[(dataframe.TF == this_tf) & (dataframe.enhancer_overlap == False)].shape[0]
    false_TF_true_enhancer = dataframe[(dataframe.TF != this_tf) & (dataframe.enhancer_overlap == True)].shape[0]
    false_TF_false_enhancer = dataframe[(dataframe.TF != this_tf) & (dataframe.enhancer_overlap == False)].shape[0]
    if (true_TF_false_enhancer/false_TF_false_enhancer) != 0: 
        odds_ratio = (true_TF_true_enhancer/false_TF_true_enhancer)/(true_TF_false_enhancer/false_TF_false_enhancer)
    else: 
        odds_ratio = np.nan
   
    return odds_ratio

# =============  MAIN =============
for this_elem in CHROMHMM_ENHANCER_INTERACTION_FILE: 
    print(this_elem)
    ### LOAD AND TIDY DATA
    TE_ANALZED = this_elem
    df = pd.read_csv(os.path.join(ROOT_PATH_chromHMMEnhancer, CHROMHMM_ENHANCER_INTERACTION_FILE[TE_ANALZED]), sep="\t", header=None, usecols=np.arange(13))

    df.columns = ["TF_chr", "TF_start","TF_end", "TF","TE_coordinate",
                "motif_score", "strand", "p_value","q_value","motif_seq","enhancer_chr","enhancer_start","enhancer_end"]

    df = df.drop([ "motif_seq", "p_value", "strand"], axis=1)
    df['enhancer_overlap'] = np.where(df['enhancer_chr'] ==".", False, True)

    ### CALCULATE TE-ENHANCER OVERLAP 
    num_TE = len(df.TE_coordinate.unique())
    TE_no_enhancer_once =  set(df.loc[df['enhancer_overlap']== False,'TE_coordinate'].unique())
    TE_yes_enhancer =  set(df.loc[df['enhancer_overlap']== True,'TE_coordinate'].unique())
    TE_with_no_enhancer = TE_no_enhancer_once.difference(TE_yes_enhancer)

    # =============  WRITE OUTPUT SUMMARY =============
    try:
        output_dir_name = "output_{}_{}".format(TE_ANALZED, time.strftime("%d_%m_%Y_%S"))
        os.mkdir(os.path.join(OUTPUT_DIR,output_dir_name))
    except  FileExistsError:
        print("output_dir_name exists, contents will be overwritten")

    summary_file_name = os.path.join(OUTPUT_DIR,output_dir_name,"output_{}.txt".format(TE_ANALZED))
    with open(summary_file_name, 'w') as fh:
        fh.write("SUMMARY OUTPUT FROM TF_CLUSTERING ANALYSIS OF {}\n".format(TE_ANALZED))
        fh.write("-------------------------------------------------------------------\n")
        fh.write("\n")
        fh.write("Analysis generated on {} using script {}.\n".format(time.strftime("%d_%m_%Y"),os.path.realpath(__file__)))
        fh.write("Total Number of TE:{}\n".format(num_TE))
        fh.write("Total Number of TE with Enhancer Overlap:{} (%allTE: {})\n".format(len(TE_yes_enhancer), len(TE_yes_enhancer)/num_TE)*100)
        
    # =============  CALCULATE ENRICHMENT =============
    unique_tf = df.TF.unique()

    # subset TE to those that have any overlap with an enhancer
    unique_TE_with_enhancers = set(df[df['enhancer_overlap']==True].TE_coordinate)
    enhancer_df = df.loc[df['TE_coordinate'].isin(unique_TE_with_enhancers)]

    # calc OR 
    store_OR_enhancers = np.empty([len(unique_tf)])
    store_OR = np.empty([len(unique_tf)])

    for this_index,this_tf in enumerate(unique_tf): 
        print(this_index)
        odds_ratio = calc_OR(df, this_tf)
        store_OR[this_index] = odds_ratio
        
        odds_ratio_enh = calc_OR(enhancer_df, this_tf)
        store_OR_enhancers[this_index] = odds_ratio_enh


    # initialize dataframe to hold OR 
    series_or = pd.Series(store_OR, index=unique_tf).sort_values()
    series_or_enh = pd.Series(store_OR_enhancers, index=unique_tf).sort_values()

    # plot 
    fig, ax = plt.subplots(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    bar_ax = ax.bar(np.arange(len(series_or)), series_or,0.5)
    _ = plt.title("TF Enrichment for Enhancer Overlap for all {}".format(TE_ANALZED))
    _ = plt.xlabel("Transcription Factors")
    _ = plt.ylabel("Odds Ratio")
    plt.savefig(os.path.join(OUTPUT_DIR, output_dir_name, "{}_TF_TE_Enhancer_Enhrichment_{}.eps".format(time.strftime("%d_%m_%Y"),TE_ANALZED)))
    plt.close()

    fig, ax = plt.subplots(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
    bar_ax = ax.bar(np.arange(len(series_or_enh)), series_or_enh,0.5)
    _ = plt.title("TF Enrichment for Enhancer Overlap for {} that Overlap Enhancers".format(TE_ANALZED))
    _ = plt.xlabel("Transcription Factors")
    _ = plt.ylabel("Odds Ratio")
    plt.savefig(os.path.join(OUTPUT_DIR,output_dir_name, "{}_TF_TEwEnhancer_Enhancer_Enhrichment_{}.eps".format(time.strftime("%d_%m_%Y"),TE_ANALZED)))
    plt.close()

    data_to_plot = series_or_enh[series_or_enh > 2]
    if len(data_to_plot) > 0: 
        fig, ax = plt.subplots(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
        bar_ax = ax.bar(np.arange(len(data_to_plot)), data_to_plot,0.5)
        _ = plt.xticks(np.arange(len(data_to_plot)), data_to_plot.index)
        _ = plt.title("TF with OR > 2.0 for Enhancer Overlap for {} that Overlap Enhancers".format(TE_ANALZED))
        _ = plt.xlabel("Transcription Factors")
        _ = plt.ylabel("Odds Ratio")
        plt.savefig(os.path.join(OUTPUT_DIR, output_dir_name, "{}_enrichedTF_TEwEnhancer_Enhancer_Enhrichment_{}.eps".format(time.strftime("%d_%m_%Y"),TE_ANALZED)))
        plt.close()
