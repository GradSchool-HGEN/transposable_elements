#!/bin/python
# This script will analyze the output of TE_annotate_hmm.py for various annotations for MER20. 
#
#
#
# Abin Abraham
# created on: 2018-01-04 08:21:52




import os
import numpy as np 
import matplotlib.pyplot as plt
from datetime import date 
from scipy import signal

rootpath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/TE_annotate_hmm_output"     
saveFigDir = "/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/MER20_hmmAnno"

def getAnnoPath(annotation):
    annoFiles = {"enhan":"THE1B_fantEnh_mappedData_mult_2018-01-09.npy",
                "phyloP":"THE1B_phyloP100_mappedData_mult_2018-01-09.npy",
                "TF":"THE1B_TFmotif_mappedData_mult_2018-01-09.npy" }
                # "phastCon":"THE1B_phastCon100_mappedData_mult_2018-01-09.npy",
    
    return os.path.join(rootpath, annoFiles[annotation])

def loadDict(filepath):
    unpack = np.load(filepath)
    annoDict = unpack.item()
    return annoDict


def getTwoSignals(seqtosearch,dict1, dict2):
    twoSigs = dict() 
    for seq in seqtosearch: 
        sig1 = dict1.get(seq, None)
        sig2 = dict2.get(seq, None)
    
        if (sig1 != None) and (sig2 != None): 
            twoSigs[seq] = [sig1, sig2]

    return twoSigs        

def plot_corr(sig1, sig2,num): 
    plt.subplot(311)
    plt.plot(sig1, 'r', label="sig1")
    plt.plot(sig2, 'b', label="sig2")
    plt.legend()

    plt.subplot(312)
    plt.plot(signal.correlate(sig1,sig2,mode='same'),label='correlate')
    plt.legend()
    
    plt.subplot(313)
    plt.plot(signal.convolve(sig1,sig2,mode='same'),label='convolve')
    plt.legend()
    plt.savefig("/dors/capra_lab/abraha1/projects/transposable_elements/results/convolve/"+str(num)+".png")

def plot_raw(annotDict, element_NAME, annotation_NAME):

    # annotDict = phyloPArray
    # element_NAME = "MER20"
    # annotation_NAME = "phast"
    consLen  = len(list(annotDict.values())[0].ravel())
    allArray = np.zeros((1,consLen))

    counter = 0  # count number of sequences that are not included in plot
    a1 = plt.figure(figsize=(10,10))

    #plot raw data 
    plt.subplot(3,1,1)
    ktrck =0 

    for k,v in annotDict.items(): 
        if ktrck == 100: 
            break
        print("plotting seq {} out of {} for {} with annotation {}".format(ktrck, len(annotDict),element_NAME,annotation_NAME))
        ktrck += 1
        v = v.ravel()
        allArray = np.vstack((allArray,v))
        plt.plot(np.arange(consLen), v.ravel(), 'b-', linewidth = 0.7, alpha=0.1)
        if ~all(np.isnan(v)): counter +=1
    plt.grid()
    plt.title('Genomic Instances of {} annotated with  {} (n={})'.format(element_NAME, annotation_NAME, str(counter)), weight='bold')
    plt.ylabel(annotation_NAME)
    [x1,x2,y1,y2] = plt.axis()
    # frame1 = plt.gca()
    # frame1.axes.xaxis.set_ticklabels([])
    
    #plot mean and stdv
    allArray = allArray[1:] #remove initializing zeros
    meanConsArray = np.nanmean(allArray, axis=0)
    stdConsArray = np.nanstd(allArray, axis=0)
    plt.subplot(3, 1, 2) #MEAN and STD 
    plt.plot(np.arange(consLen), meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(np.arange(consLen), meanConsArray + 2*stdConsArray, 'k:', linewidth = 0.7, alpha=1)
    plt.plot(np.arange(consLen), meanConsArray - 2*stdConsArray, 'k:', linewidth = 0.7, alpha=1)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.title("Mean (Red) and 2*SD (black)", weight='bold')
    plt.ylabel(annotation_NAME)
    # frame1 = plt.gca()
    # frame1.axes.xaxis.set_ticklabels([])

    #plot coverage at each base
    plt.subplot(3, 1, 3) 
    count_allArray = allArray
    count_allArray[~np.isnan(count_allArray)] = 1 
    sumConsArray = np.nansum(count_allArray, axis=0)
    plt.plot(np.arange(consLen), sumConsArray, 'k-', linewidth = 1.0, alpha=1)
    plt.grid()
    plt.title('Data Coverage at Each Base', weight='bold')
    plt.xlabel('consensus bases')
    plt.ylabel("Num Values at Each Base")
    a1.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.3)
    pltname = "{}_{}_rawHMManno_{}.eps".format(element_NAME, annotation_NAME, date.today())
    # plt.savefig(os.path.join(saveFigDir,pltname))
    plt.show()

def run_plot_raw(annoList):
    for oneAnno in annoList: 
        thisDict = loadDict(getAnnoPath(oneAnno))
        plot_raw(thisDict, elem, oneAnno)

def plot_perHit( sig1, sig2, sig1Name, sig2Name,title, *savePath):
    plt.plot(sig1,'r',label=sig1Name)
    plt.plot(sig2,'b',label=sig2Name)
    plt.title(title)
    plt.legend()
    if savePath: 
        plt.savefig(savePath[0])

#-------
# main
#-------

### create raw data plots: 
elem = "THE1B"
annoList =  ["phyloP"]
# annoList =  ["phastCon","TF", "enhan","phyloP"]
run_plot_raw(annoList)

'''
### load annotations 
# phastArray = loadDict(getAnnoPath("phastCon"))
TFArray = loadDict(getAnnoPath("TF"))
enhArray = loadDict(getAnnoPath("enhan"))
phyloPArray = loadDict(getAnnoPath("phyloP"))

## compare those with enhancer overlap to those without overlap
seqToSearch = set(enhArray.keys())
twoSigDict = getTwoSignals(seqToSearch, enhArray, TFArray)


#plot enhancer activity and TF activity 
for k,v in twoSigDict.items():
    print()
    print(k) 
    print(v)

    plt.figure()
    sig1 = v[0]
    sig2 = v[1]
    sig1[np.isnan(sig1)] = 0
    sig2[np.isnan(sig2)] = 0

    # plt.plot(sig1,'r',label='Enhancer Activity')
    # plt.plot(sig2,'b',label='TF Activity')
    # plt.title(k)
    # plt.legend()
    # plt.savefig("/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/MER20_enh_TF_singleSeq/"+k+"_2.eps")


# sig1, sig2 = twoSigDict[ 'chr6:110361036-110361252']
# sig1[np.isnan(sig1)] = 0 
# sig2[np.isnan(sig2)] = 0 


### compare two signals 

# 

# to do 
#   same plots as above, but split into w/ enhc overlap and without 
'''