#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-04 00:20:23

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import train_test_split
from sklearn import metrics
from sklearn.cross_validation import cross_val_score
from sklearn.decomposition import PCA
from sklearn.metrics import precision_recall_curve

ROOT_PATH = "/Users/abin-personal/Google Drive/GD_transfer/data"
RAW_FILE = "intersect_HERV_TFBS_facet_enhancers.tsv"

# =============  functions =============
def calc_overlap(df):
    numOverlap = len(set(np.arange(df['TE_start'], df['TE_end'])).intersection(set(np.arange(df["enhc_start"],df["enhc_end"]))))
    totalOverlap  = numOverlap/(df['TE_end'] - df['TE_start'])
    print(totalOverlap)
    return totalOverlap

# =============  main =============
df = pd.read_csv(os.path.join(ROOT_PATH,RAW_FILE), sep ='\t',header=None,index_col=None)

### TIDY DATA
df.columns

df.columns = ["TE_chr", "TE_start", "TE_end", "TE", "TE_family",  "TF_chr", "TF_start", "TF_end", "TF", "TE_info", "enhc_chr", "enhc_start", "enhc_end","enhc_overlap" ]
df['enhancer_id'] = df['enhc_chr'].map(str) + ":" +  df['enhc_start'].map(str)  + "-" + df['enhc_end'].map(str) 
df['TE_overlap'] = np.where(df['TE_chr']!=".", True, False)
print(df.head())

### create input features 
X = df.groupby(['enhancer_id', 'TF']).count()[['TE_info']].unstack().fillna(0).values
Y = df.drop_duplicates(['enhancer_id'])['TE_overlap']

# pca = PCA(n_components=2, svd_solver='full')
# pca.fit(x)                 
# print(pca.explained_variance_ratio_)  
# print(pca.singular_values_)  



## TRAIN DATA 
clf = RandomForestClassifier(n_estimators=1, random_state=0)

Xtrain, Xtest, ytrain, ytest = train_test_split(X, Y, random_state=0)
clf = RandomForestClassifier(max_depth=11)
clf.fit(Xtrain, ytrain)
ypred = clf.predict(Xtest)

cv = cross_val_score(clf, X, Y, cv=10)
cv.mean()

precision, recall, _ = precision_recall_curve(ytest, ypred)

