#!/bin/python
# This script will plot TF motifs onto a consensus sequnce for a TE. Consenus derived from pHMM models from DFAM using HMMMER tools. 
#
#
#
# Abin Abraham
# created on: 2018-01-12 19:15:14

import pandas as pd
import numpy as np 
from Bio import SeqIO
from cycler import cycler
import datetime


element = "MER20"
TFannotationfile  = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/fimo-dfam/fimo-out-consensusSeq_{}/formatted_fimo_out_{}.tsv".format(element,element)
outputDir = '/dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align'
def getConsensusLength(element):
    #get the length of the consensus sequence 
    TEconsensusFasta_File = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-input_{}/consensus_{}.fa".format(element,element)
    consensusSeq = SeqIO.read(TEconsensusFasta_File,'fasta')
    
    return len(consensusSeq.seq)


def mapTF_toConsensus(element):
    # chr, start_relative, end_relative, TF, score SeqName, start, stop, strand, pvalue, qvalue,matched, TF quality
    #create one row per TF
   # TFannotationfile = 
    df = pd.read_table(TFannotationfile, sep='\t')
    df.columns = [ "start_relative", "end_relative", "score", "TF",  "SeqName"," start", "stop", "strand", "pvalue", "qvalue","matched", "TFquality"]
    gbTF = df.groupby('TF')
    consensusLength = getConsensusLength(element)
    consArray = np.zeros((1, consensusLength))
    tmp_consArray = np.zeros((1,consensusLength))
    mask_consArray = np.zeros((1,consensusLength))
    
    TFlist = [ list(gbTF)[i][0] for i in range(len(gbTF))]

    #load array with score of each TF, one row per TF
    for i in range(len(gbTF)):
        spec = list(gbTF)[i][1][["start_relative","end_relative", "score"]].values
        print(i)
        for s,e,v in spec: 
            #check if these indices/bases already had a value, if so take the average 
            tmp_consArray[0, int(s):int(e)]= (tmp_consArray[0, int(s):int(e)] + v)
            mask_consArray[0, int(s):int(e)]+=1

        tmp_consArray[mask_consArray != 0] = tmp_consArray[mask_consArray != 0]/mask_consArray[mask_consArray != 0]
        consArray = np.vstack((consArray, tmp_consArray))
        tmp_consArray = np.zeros((1,consensusLength))
        mask_consArray = np.zeros((1,consensusLength))

    consArray = consArray[1:,:] # reomve initiliziling row

    return(consArray)


def plot_TFtoConsensus(map_consArray, element_NAME, outputDir):

    import matplotlib.pyplot as plt
    from cycler import cycler 
    import os 

    consensusLength = len(map_consArray[0])
    consBases = np.arange(0.0, consensusLength, 1, dtype=int)
    meanConsArray = np.mean(map_consArray, axis=0)
    stdConsArray = np.std(map_consArray, axis=0)

    plt.subplot(2, 1, 1) #RAW DATA 
    plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
    for i in range(len(map_consArray)): #each row in map_consArray gets a individual plot
        plt.plot(consBases, map_consArray[i], linewidth = 1, alpha=0.3)
    plt.grid()
    plt.title('TF Motif Indentified on Consensus'+ element_NAME + '\n top (one line per TF), bottom (average across all TF)')
    plt.ylabel('TF Motif Score')
    [x1,x2,y1,y2] = plt.axis()

    plt.subplot(2, 1, 2) #MEAN and STD 
    plt.plot(consBases, meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(consBases, meanConsArray + 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.plot(consBases, meanConsArray - 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.xlabel('consensus bases (1start)')
    plt.ylabel('TF Motif Score')

    plt.subplots_adjust(hspace=0.3)
    saveFigName = os.path.join(outputDir, element_NAME+"_TFmotif_on_Consensus_" + str(datetime.date.today())+ ".png")
    plt.savefig(saveFigName)
    # plt.show()
    
    return()