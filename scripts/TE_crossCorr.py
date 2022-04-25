#!/bin/python
# This script will calculate the cross-correlation between two annotations for a given transposable element. 
#   Uses the dfam database. 
#
#
# Abin Abraham          
# created on: 2017-12-18 12:03:45

import numpy as np
import pickle 

#-------
# functions
#-------

def plotRawData(annotationArray, annoName, elemName):
    import matplotlib.pyplot as plt


    consBases = np.arange(0.0, annotationArray.shape[1], 1, dtype=int)
    plt.figure(figsize=(8, 5), dpi=160)

    for i, oneLocus in enumerate(annotationArray):
        plt.plot(consBases, oneLocus, linewidth = 1, alpha=0.4, label="an instance of  "+elemName)
    
    plt.grid()
    plt.title(annoName +" annotation for all genomic " + elemName)
    plt.ylabel(annoName)
    plt.xlabel("aligned to consensus")
    plt.show()
    


#-------
# main
#-------
resultsPATH = "/dors/capra_lab/abraha1/projects/transposable_elements/results/"


#laod data -using pickle 
with open(resultsPATH+'MER20_phyloP100_2017-12-18.p', 'rb') as fp:
    phyloP100_data = pickle.load(fp)

with open(resultsPATH+'MER20_TF_on_Enhancer_2017-12-18.p', 'rb') as fs:
    enhTF_data = pickle.load(fs)


phastCon_Array = [k for k in phyloP100_data.values()]
phastCon_Array = np.asarray(phastCon_Array)

enhTF_Array = [k for k in enhTF_data.values()]
enhTF_Array = np.asarray(enhTF_Array)

plotRawData(phastCon_Array, 'phastCon', "MER20")
plotRawData(enhTF_Array, 'enhancer-TF-overlap', "MER20")

# find both signals for a given TE element
sharedSeq = phyloP100_data.keys() & enhTF_data.keys()

# thisSeq = 'chr1:14054921-14055081'
import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5), dpi=160)

for thisSeq in sharedSeq: 
    
    phastConValues = phyloP100_data[thisSeq]/max(phyloP100_data[thisSeq])
    enhTFValues = enhTF_data[thisSeq]/max(enhTF_data[thisSeq])
    np.correlate(phastConValues,enhTFValues)

    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(8, 5), dpi=160)
    plt.plot( np.arange(0.0, len(phastConValues), 1, dtype=int), phastConValues, linewidth = 1, alpha=0.8, label="Conservation")
    plt.plot( np.arange(0.0, len(enhTFValues), 1, dtype=int), enhTFValues, linewidth = 1, alpha=0.8, label="enhTFValues")

    plt.grid()
    plt.title("Comparing Conservation and enhancer TF for one instance of MER20")
    plt.ylabel('annotation Score')
    plt.legend()
    plt.xlabel("aligned to consensus")
    plt.show()
        