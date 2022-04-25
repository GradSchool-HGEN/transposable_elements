#!/bin/python
# train, test, evaluate a random forest to predict enhancer activity in a specific tissue. 
# given a TE with TF binding sites on TE as features.
#
#
#
# Abin Abraham
# created on: 2018-02-04 00:20:23


import os
import datetime
import pickle 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.cross_validation import train_test_split
from sklearn.cross_validation import cross_val_score
from sklearn.externals import joblib

from imblearn.over_sampling import RandomOverSampler

from treeinterpreter import treeinterpreter as ti 
from scipy.cluster.hierarchy import dendrogram, linkage

#=========== PATHS ===========
ROOT_PATH = "/Users/abin-personal/Google Drive/GD_transfer/data/"
SAVE_PATH = "/Users/abin-personal/Google Drive/GD_transfer/output/"
RAW_FILE = "intersect_HERV_TFBS_facet_enhancers_withTissue.tsv"


# =============  functions =============
def calc_overlap(df):
    num_overlap = len(set(np.arange(df['TE_start'], df['TE_end'])).intersection(set(np.arange(df["enhc_start"],df["enhc_end"]))))
    totalOverlap  = num_overlap/(df['TE_end'] - df['TE_start'])
    print(totalOverlap)
    return totalOverlap

def plot_pca_rawdata(pca_df): 
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA\n (0=No Enhancer Overlap, 1= Enhancer Overlap)', fontsize = 20)
    targets = [0,1]
    colors = ['k', 'r']

    for target, color in zip(targets, colors): 
        indicies_to_keep = pca_df['enhc_overlap'] == target
        ax.scatter(pca_df.loc[indicies_to_keep, 'pca1'], pca_df.loc[indicies_to_keep, 'pca2'], c = color, s = 10, alpha =0.5)

    ax.legend(targets)
    ax.grid()
    plt.show()

def plot_feature_summary(features_summary):    
    plt.figure(figsize=(18,6)) 
    plt.title("Feature Importances for Random Forest Classification of Enhancer Overlap Given LTRs w/ TFBS (Resampled Data)")
    plt.bar(range(len(features_summary)), features_summary[:,1], color='k')
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xticks(range(len(features_summary)),features_summary[:,0], rotation='vertical', fontsize=9)
    plt.savefig(os.path.join(ROOT_PATH,"{}_FeatureImportance_RandomForest_ResampledData.eps".format(datetime.datetime.today().strftime('%Y-%m-%d'))))
    plt.show()

# =============  main =============
rdf = pd.read_csv(os.path.join(ROOT_PATH,RAW_FILE), sep ='\t',header=None,index_col=None, dtype='str')

### TIDY UP DATA
rdf.columns = ["TE_chr", "TE_start", "TE_end", "TE", "TE_family", "TF_chr", "TF_start", "TF_end", "TF", "TE_info", "enhc_chr", "enhc_start", "enhc_end", "enhc_tissue", "enhc_overlap" ]
rdf['TE_id'] = rdf['TE_chr'].map(str) + ":" +  rdf['TE_start'].map(str)  + "-" + rdf['TE_end'].map(str) 
rdf[['TE_start', 'TE_end','enhc_start', 'enhc_end']].apply(pd.to_numeric)
rdf['enhc_overlap'] = np.where(rdf['enhc_chr']!=".", True, False)

### CREATE TRAIN AND TEST SET
# grouby TE and TF binding site motifs 
df = rdf.groupby(['TE_id', 'TF']).count()[['TE_info']].unstack().fillna(0)
mdf = rdf.groupby(['TE_id']).count()[['enhc_tissue']]
X = df.values
Y = mdf.values
Y = Y.astype(int).ravel()

df.to_csv(os.path.join(SAVE_PATH,'TFdf'),sep="\t", header=False)
mdf.to_csv(os.path.join(SAVE_PATH,'overlap_df'),sep="\t", header=True)

### VISUALIZE DATA 
pca = PCA(n_components=2)
X_r = pca.fit_transform(X)
pca_df = pd.DataFrame(data =np.concatenate((X_r, Y.reshape(-1,1)), axis=1), columns = ['pca1', 'pca2', 'enhc_overlap'])
print(pca.explained_variance_)
# plot_pca_rawdata(pca_df)

### TRAIN DATA 
Xtrain, Xtest, ytrain, ytest = train_test_split(X, Y, random_state=0)
clf = RandomForestClassifier(max_depth=11)
clf.fit(Xtrain, ytrain)
ypred = clf.predict(Xtest)

### EVALUATE MODEL
# metrics.confusion_matrix(ytest,ypred)
# cv = cross_val_score(clf, X, Y, cv=10)
# print("10 fold cross validation mean: {}.".format(cv.mean()))

#plot roc  
# fpr, tpr, thresholds = metrics.roc_curve(ytest, ypred, pos_label=1)
# plt.plot(fpr,tpr, lw =1)

### OVERSAMPLE, correct for imbalanced data set 
ros = RandomOverSampler(random_state=0)
X_resampled, Y_resampled = ros.fit_sample(X, Y)

### TRAIN RESAMPLED DATA 
Xtrain_r, Xtest_r, ytrain_r, ytest_r = train_test_split(X_resampled, Y_resampled, random_state=0)
clf_re = RandomForestClassifier(max_depth=11,random_state=0)
clf_re.fit(Xtrain_r, ytrain_r)
ypred_r = clf_re.predict(Xtest_r)

### EVALUATE MODEL
# cv_r = cross_val_score(clf_re, X_resampled, Y_resampled, cv=10)
# print("10 fold cross validation mean of oversampled data : {}.".format(cv_r.mean()))

print(metrics.confusion_matrix(ytest_r,ypred_r))
print(metrics.classification_report(ytest_r,ypred_r))
print(metrics.precision_score(ytest_r,ypred_r))

### SAVE MODEL 
# joblib.dump(clf_re, os.path.join(SAVE_PATH,"{}_classify_enhancerOverlap_OverSampled.pkl".format(datetime.datetime.today().strftime('%Y-%m-%d'))))

### PCA on resampled data  
pca_r = PCA(n_components=2)
X_r = pca_r.fit_transform(X_resampled)
pca_df_r = pd.DataFrame(data =np.concatenate((X_r, Y_resampled.reshape(-1,1)), axis=1), columns = ['pca1', 'pca2', 'enhc_overlap'])
print(pca_r.explained_variance_ratio_)
# plot_pca_rawdata(pca_df_r)

### FEATURE IMPORTANCE
feature_labels = df.columns.levels[1].values
feature_importance = clf_re.feature_importances_
features_summary = np.sort(np.concatenate((feature_labels.reshape(-1,1), feature_importance.reshape(-1,1)), axis =1), axis=0)
# plot_feature_summary(features_summary)

### TREE INTERPRETATION
# store_path = np.empty(Xtest_r.T.shape)
# for ind, instance in enumerate(Xtest_r): 
#     prediction, _, contributions = ti.predict(clf_re, instance.reshape((-1,106)))

#     if prediction[0][0] > prediction[0][1]:
#         store_path[:, ind] = contributions[0][:,0]
#     elif prediction[0][0] < prediction[0][1]:
#         store_path[:, ind] = contributions[0][:,1]
    
#     print("{} of {}".format(ind, len(Xtest_r)))
#     # if ind == 1000: break 

# pickle.dump(store_path, open('node_weights_resampled_test.pickle','wb'))

### plot node weights as a function of classicfication 
node_weight = pickle.load( open( 'node_weights_resampled_test.pickle', "rb" ) )
# ytest_r
# ypred_r
true_positive_mask  = np.logical_and((ytest_r==1), (ypred_r==1))
true_negative_mask  = np.logical_and((ytest_r==0), (ypred_r==0))
false_positive_mask = np.logical_and((ytest_r==0), (ypred_r==1))
false_negative_mask = np.logical_and((ytest_r==1), (ypred_r==0))

true_positive_node_means = np.mean(node_weight[:,true_positive_mask ], axis=1)
true_postive_node_std = np.std(node_weight[:,true_positive_mask ], axis=1)
true_negative_node_means = np.mean(node_weight[:,true_negative_mask ], axis=1)
true_negative_node_std = np.std(node_weight[:,true_negative_mask ], axis=1)

#plot
fig, ax = plt.subplots(figsize=(14,5))
width = 1
index = np.arange(len(true_positive_node_means))
tp_bars = ax.bar(index, true_positive_node_means, width, alpha=1, color='b', yerr=true_postive_node_std, label='TP')
# tn_bars = ax.bar(index+width, true_negative_node_means, width, alpha=1, color='r', yerr=true_negative_node_std, label='TN')
ax.set_xlabel('TFs')
ax.set_ylabel('Node Weights')
ax.set_title('Node Weights by Classification')
ax.set_xticks(index+width)
ax.set_xticklabels(feature_labels,rotation='vertical', fontsize=7)
ax.legend()
fig.tight_layout()
plt.savefig(os.path.join(SAVE_PATH,"Feb_11_2018_NodeWeights_TruePos.eps"))
plt.show()

# filter store_path to select True Positive compared to true negatives 
# store_path[:, ytest_r>0]
# dendrogram(linkage(store_path), labels=feature_labels)


# # one example: 
# print("Evaluating Node Feature Importance....\n")
# one_instance = X[1:2]
# clf_re.predict_proba(one_instance)

# prediction, bias, contributions = ti.predict(clf_re, one_instance)
# print("Prediction ",  prediction)
# print("Bias (trainset prior) " , bias)
# print("Feature contributions:")
# for c, feature in zip(contributions[0], 
#                              feature_labels):
#     print(feature," ", c)



