#!/bin/python
# This script contains a functions for plotting TF motifs on a TE mapped to its consensus. 
#
# Dependencies: 
#       1) /dors/capra_lab/abraha1/results/TE_TF-heatmaps/formated_fimo.txt (1start, inclusive)
# 
# Abin Abraham
# created on: 2017-12-11 16:40:01


# element = "MER20-consensus"

import pandas as pd 
import numpy as np 
import subprocess
import datetime

def listElementNames():
    TFannotationfile = "/dors/capra_lab/abraha1/results/TE_TF-heatmaps/fimo_out/formated_fimo.txt"
    getTE = subprocess.check_output(['cut','-f', '1' , TFannotationfile ], universal_newlines=True)
    elementList = set(getTE.split('\n'))
    return(elementList)
    
def getConsensusLength(element):
    #get the length of the consensus sequence 
    TEconsensusFasta_File = "/dors/capra_lab/data/transposable_elements/dfam/length_Dfam.cons.fa"
    getConsLength = subprocess.check_output(['grep', element, TEconsensusFasta_File ], universal_newlines=True)
    consLength = int(getConsLength.split('\n')[0].split()[1])
    return(consLength)

def mapTF_toConsensus(element):
    TFannotationfile = "/dors/capra_lab/abraha1/results/TE_TF-heatmaps/fimo_out/formated_fimo.txt"
    filteredFile = subprocess.check_output(['grep', '-w', element, TFannotationfile ], universal_newlines=True)
    out_data = [x.split('\t') for x in filteredFile.split('\n')]
    out_data = out_data[:-1]
    out_data = [ [x[0], int(x[1]), int(x[2]), x[3], float(x[4]) ]for x in out_data]
    
    #create one row per TF
    df = pd.DataFrame(out_data, columns=["TE", "start", "end", "TF", "score"])
    gbTF = df.groupby('TF')
    consensusLength = getConsensusLength(element)
    consArray = np.zeros((1, consensusLength))
    tmp_consArray = np.zeros((1,consensusLength))
    mask_consArray = np.zeros((1,consensusLength))
    
    TFlist = [ list(gbTF)[i][0] for i in range(len(gbTF))]

    #load array with score of each TF, one row per TF
    for i in range(len(gbTF)):
        spec = list(gbTF)[i][1][["start","end", "score"]].values
        # check for overlap --> NEED TO MODIFY, this is not actually checking for overlap
        
        for s,e,v in spec: 
            #check if these indices/bases already had a value, if so take the average 
            tmp_consArray[0, int(s):int(e)+1]= (tmp_consArray[0, int(s):int(e)+1] + v)
            mask_consArray[0, int(s):int(e)+1]+=1

        tmp_consArray[mask_consArray != 0] = tmp_consArray[mask_consArray != 0]/mask_consArray[mask_consArray != 0]
        consArray = np.vstack((consArray, tmp_consArray))
        tmp_consArray = np.zeros((1,consensusLength))
        mask_consArray = np.zeros((1,consensusLength))

    consArray = consArray[1:,:]

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


