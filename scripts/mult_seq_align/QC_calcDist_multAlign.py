#!/bin/python
# This script will calculate the distance between a subset of all aligned sequences. 
#
#
#
# Abin Abraham
# created on: 2017-12-20 12:46:44


from scipy import spatial as sp
from Bio import AlignIO 
import itertools
import numpy as np
import matplotlib.pyplot as plt



def calcHamming(seq1, seq2):
    val = sp.distance.hamming(seq1,seq2)
    return val*len(seq1)


def benchMark(seqRec,storeHAM,seqIND):
    for i,k in enumerate(seqIND):
        hamDist = calcHamming(list(seqRec[k[0]].seq), list(seqRec[k[1]].seq)) 
        storeHam[i] = hamDist
        print(str(i) + ":" + str(k))

def calcDist(alignedFile, numSubSample):
    ### LOAD 
    # alignedFile = "/Volumes/capra_lab/abraha1/projects/transposable_elements/data/dfam/temp.muscle.out"
    # alignedFile  ="/Users/abin-personal/Desktop/transfer/data/temp.muscle.out"
    # alignedFile = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/temp.muscle.out"

    multAlign = AlignIO.read(alignedFile, 'fasta')

    #numSubSample = 10
    randomIndex1 = np.random.choice(len(multAlign), numSubSample, replace=False)
    randomIndex2 = np.random.choice(len(multAlign), numSubSample, replace=False)
    sameIndex_toRemove = [i for i, x in enumerate(randomIndex1==randomIndex2) if x]
    randomIndex1 = np.delete(randomIndex1, sameIndex_toRemove)
    randomIndex2 = np.delete(randomIndex2, sameIndex_toRemove)
    seqInd = [randomIndex1, randomIndex2]

    storeHam = np.empty(len(seqInd[0]))

    for i in range(len(storeHam)):
        firstIndex = int(seqInd[0][i])
        secondIndex = int(seqInd[1][i])
        hamDist = calcHamming(list(multAlign[firstIndex].seq), list(multAlign[secondIndex].seq)) 
        storeHam[i] = hamDist

    return HammingDistances



'''
#Violing Plot of Store Ham Values 
fontsz =10 
data = storeHam

fig, ax = plt.subplots()
ax.violinplot(data, points=100, widths=0.1, showmeans=True, showextrema=True, showmedians=True)
ax.set_title("Distribution of Hamming Scores for \nRandom Subsample (" + str(numSubSample)+ " no replacement) of Aligned Sequences")

plt.show()

'''