#!/bin/python
# This script will # of bases aligned over the multalig consensus
#
#
#
# Abin Abraham
# created on: 2017-12-20 11:09:41

#consider 2-bit logograms..

import getMappingDict as md 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

rootPATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/multAlign/"
align_FILE = "muscle-output_MER20_maxiter2"
align_FILE2 = "muscle-output_MER20_fastest"
align_FILE3 = "muscle-output_subset_MER20"
align_FILE4 = "mafft-v1_MER20_subset"
align_FILE5 = "mafft-default_MER20"

# create dicitonary that has mapping to alignmnet
mappedDict = md.createMapDict(rootPATH+align_FILE5)
consensusLenght = list(mappedDict.values())[0][1]

sumAlign = np.zeros(consensusLenght)
for v in mappedDict.values():
        sumAlign[v[0]] += 1

#=========== tag sequences with > 20 bases found in regions with low coverage ===========
#plot histogram of coverage 
axS = sns.distplot(sumAlign, hist=False, rug=True)
axS.set_xlabel = "Distribution of Number of Sequences Hits to A Base in Alignment"
axS.set_ylabel = "Density"
plt.show()


# define threshold for bases in alignment with low hits 
highCutoff = 1000
lowCutoff = 200 

belowInd_high = set(np.where(sumAlign<highCutoff)[0])
belowInd_low = set(np.where(sumAlign<lowCutoff)[0])


numBases_inDessert_threshold = 20
storeSeq = set()
for i,v in enumerate(mappedDict):
                thisSeq = set(list(mappedDict.values())[i][0])
                numDesertBases = len(belowInd_high.intersection(thisSeq))
                if numDesertBases > numBases_inDessert_threshold: 
                        storeSeq.add(list(mappedDict.keys())[i])
 
# filepath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/multAlign/MER20_more20desertBases.txt"
# with open(filepath, 'w') as file1:
#     for item in storeSeq:
#         file1.write("{}\n".format(item))


#=========== Analyzing Low Coverage Regions in Alignment ===========

### plot number of bases vs. num of bases aligned to desert regions 
storeNumBases = np.empty([len(mappedDict),2])
for i,v in enumerate(mappedDict):
                thisSeq = set(list(mappedDict.values())[i][0])
                numDesertBases = len(belowInd_high.intersection(thisSeq))
                numBasesTE = len(thisSeq)
                storeNumBases[i][0] = numBasesTE
                storeNumBases[i][1] = numDesertBases

figB, axB = plt.subplots()
axB.scatter(storeNumBases[...,0],storeNumBases[...,1], s=5, alpha=0.4)
axB.set_xlabel("num bases in TE", fontsize=15)
axB.set_ylabel("num bases in desert regions in Alignment", fontsize=15)
axB.set_title('Num Bases in Desert Regions in Alignment Per TE for MER20')
axB.grid(True)
plt.show()

#same plot but with histograms
sns.set_style("whitegrid")
sns.set(style="ticks")
sns.jointplot(storeNumBases[...,0], storeNumBases[...,1], kind="scatter", color="#4CB391", s=5).set_axis_labels("num bases in a TE","num bases found in desert regions in alignment for a TE")
plt.show()

### Calculate Number of TEs with X% of bases in desert region
matchCutoff_array = np.linspace(0,0.5,25)
storeValue = np.empty([len(matchCutoff_array),2])
# storeValue = np.empty([len(100),2])
for ind,m in enumerate(matchCutoff_array): 
        matchCutoff = m
        desertCount = 0
        for i,v in enumerate(mappedDict):
                thisSeq = set(list(mappedDict.values())[i][0])
                perMatch = len(belowInd_high.intersection(thisSeq))/len(thisSeq)
                if perMatch > matchCutoff: 
                        desertCount += 1

        storeValue[ind][0] = m
        storeValue[ind][1] = desertCount

###plot number of sequence as function of % aligned bases found in desert regions 
fig, ax = plt.subplots()
ax.plot(storeValue[...,0], storeValue[...,1])
ax.set(xlabel='Cutoff: Num of Aligned Bases in Consensus Desert/Total Num bases in TE ', ylabel='Count of TEs Above Cutoff',
       title='Num of Bases aligned to Consensus with low coverage for MER20')
ax.grid()
plt.show()        


### Plot Coverage of  Alignment
fig, ax = plt.subplots()
ax.plot(np.arange(consensusLenght), sumAlign)
ax.set(xlabel='bases in consensus', ylabel='number of hits',
       title='MER20 Alignment to Consensus (using MAFFT-default)')
ax.grid()
#fig.savefig("MER20_alignment.png")
plt.show()