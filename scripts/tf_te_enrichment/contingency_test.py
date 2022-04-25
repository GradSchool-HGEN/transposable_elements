

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

mainPath ="/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/facet_expressed_enhancers/alltxt/allfacets/"
#mainPath1 ="/dors/capra_lab/abraha1/results/TFanalysis-human-fantom5/facet_expressed_enhancers/alltxt/allfacets"
file="all-fimo-out-withTFCoord-filtered-TE-intersect.tsv"

df=pd.read_table(file, sep="\t", header=None)
df.columns = ["chr(TFmotif)", "start(TFmotif)"," end(TFmotif)","TF", "enhancer_sequence", "q-value", "strand", "motif_quality", "motif-fimo-score", "tissue", "chr(TE)","start(TE)","end(TE)" ,"TE", "TEfam"]

sdf = df[["TF","enhancer_sequence","tissue","TE","TEfam"]]

TEfamilies_toRemove = ["Low_complexity", "Low_complexity", "RNA", "rRNA", "Satellite", "Satellite/acro", "Satellite/centr", "Satellite/telo", "scRNA", "Simple_repeat","snRNA","srpRNA", "tRNA", "Unknown"]

to_drop = sdf[sdf['TEfam'].isin(TEfamilies_toRemove)].index

nsdf = sdf.drop(to_drop)

# ---- TIDY UP DATA ----
# split data table into two: one with TE overlap; one without 

TE_Absent = nsdf[nsdf.TE=="."]
TE_Positive = nsdf[nsdf.TE!="."]

TE_Absent.drop(["TE", "TEfam", "enhancer_sequence"], axis=1, inplace=True)
TE_Positive.drop(["TE", "TEfam",  "enhancer_sequence"], axis=1, inplace=True) #dropping these columns, because already know there is TE overlap 

# ---- Create counts of TFs by tissue type 

grouped_TE_Absent = TE_Absent.groupby(["TF","tissue"]).count()
grouped_TE_Absent = grouped_TE_Absent.unstack(-1)
grouped_TE_Absent = grouped_TE_Absent.fillna(0)


grouped_TE_Positive = TE_Positive.groupby(["TF","tissue"]).count()
grouped_TE_Positive = grouped_TE_Positive.unstack(-1)
grouped_TE_Positive = grouped_TE_Positive.fillna(0)


grouped_TE_Positive.to_csv("grouped_TE_Positive", sep="\t")
grouped_TE_Absent.to_csv("grouped_TE_Absent", sep="\t")
# -- Chi Test 
# easier to do this in R to accomodate a table with dimensions larger than 2x2


TEposCount = grouped_TE_Positive.iloc[1].values
TEnegCount = grouped_TE_Absent.iloc[1].values

#TEposCount = [10, 10]
#TEnegCount = [10, 10]

p_val = chi2_contingency([TEposCount, TEnegCount])[1]
