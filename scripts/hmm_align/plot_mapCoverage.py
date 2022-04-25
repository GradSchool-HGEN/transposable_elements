#!/bin/python
# This script will plot the coverage of the agliment for TE to its consensus
#
#
#
# Abin Abraham
# created on: 2017-12-25 18:15:19

import matplotlib.pyplot as plt
import parse_nhmmerOutput as pnh
import createMapDict as crmd
import numpy as np
from Bio import SeqIO


def plot_mappingCoverage(element, nhmmerparsedOutput, consensusFile, allElemMappedDict):
    allSeq = list(nhmmerparsedOutput.keys())

    confa = SeqIO.read(consensusFile,'fasta')
    consLength = len(confa.seq)
    sumCoverage = np.zeros(consLength)
    sumDeletion = np.zeros(consLength)
    sumInsertions = np.zeros(consLength)

    for oneseq in allSeq: 
        one_elemTOcons = allElemMappedDict[oneseq]
        # bases of shared homology
        bases = crmd.getHomologBases_con(one_elemTOcons)
        bases = [int(x) for x in bases]
        delbases =  crmd.getDeletedBases_con(one_elemTOcons)
        delbases = [int(x) for x in delbases]
        ins_bases = crmd.getInsertBases_con(one_elemTOcons)
        ins_bases = [int(x) for x in ins_bases]
        sumCoverage[bases] += 1
        sumDeletion[delbases]  += 1 
        sumInsertions[ins_bases] += 1 

    # plot alignment coverage over consensus bases 
    fig, ax = plt.subplots()
    ax.plot(np.arange(consLength), sumCoverage,label ='Coverage')
    ax.plot(np.arange(consLength), sumDeletion,label ='Deletion')
    ax.plot(np.arange(consLength), sumInsertions,label ='Insertion')
    ax.set(xlabel='consensus bases for '+element, ylabel='Number of Genomics Hits Per Base',
        title=element+' Check of Mapping to Consensus')
    ax.grid()
    plt.legend()
    plt.show()



'''
parsedOutput = pnh.parseTermOut("nhmmer-output-terminal_allMER20")
allSeq = list(parsedOutput.keys())

conMER20fa = SeqIO.read('hmmemit-ouput-MER20.fa','fasta')
consLength = len(conMER20fa.seq)
sumCoverage = np.zeros(consLength)
sumDeletion = np.zeros(consLength)
sumInsertions = np.zeros(consLength)
for oneseq in allSeq: 
    #load consensus sequence 

    bases = crmd.getHomologBases_con(elemTOcons)
    bases = [int(x) for x in bases]
    delbases =  crmd.getDeletedBases_con(elemTOcons)
    delbases = [int(x) for x in delbases]
    ins_bases = crmd.getInsertBases_con(elemTOcons)
    ins_bases = [int(x) for x in ins_bases]
    sumCoverage[bases] += 1
    sumDeletion[delbases]  += 1 
    sumInsertions[ins_bases] += 1 

# plot alignment coverage over consensus bases 
fig, ax = plt.subplots()
ax.plot(np.arange(consLength), sumCoverage,label ='Coverage')
ax.plot(np.arange(consLength), sumDeletion,label ='Deletion')
ax.plot(np.arange(consLength), sumInsertions,label ='Insertion')
ax.set(xlabel='consensus bases for MER20', ylabel='Number of Genomics Hits',
       title='MER20 Alignment Coverage Across Consenus Bases')
ax.grid()
plt.legend()
plt.show()
'''

