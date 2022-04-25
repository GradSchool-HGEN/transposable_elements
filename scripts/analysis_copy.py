#!/bin/python
# This script will ...  
#
#
#
# Abin Abraham
# created on: ...⇧+⌘+I

#!/bin/python
# This script will analyze the output of TE_annotate_hmm.py for various annotations for MER20. 
#
#
#
# Abin Abraham
# created on: 2018-01-04 08:21:52


# STILLL IN DEV

import os
import numpy as np 
import matplotlib.pyplot as plt
from datetime import date 
from scipy import signal
from scipy import stats
from cycler import cycler



rootpath = "/Users/abin-personal/Desktop/transfer/data/"     
saveFigDir = "/Users/abin-personal/Desktop/transfer/figs/"
# rootpath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/TE_annotate_hmm_output"     
# saveFigDir = "/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/MER20_hmmAnno"

def getAnnoPath(annotation):
    annoFiles = {"enhan":"MER20_fantEnh_mappedData_mult_2018-01-03.npy",
                "phastCon":"MER20_phastCon100_mappedData_mult_2018-01-04.npy",
                "phyloP":"MER20_phyloP100_mappedData_mult_2018-01-03.npy",
                "TF":"MER20_TFmotif_mappedData_mult_2018-01-03.npy" }
    
    return os.path.join(rootpath, annoFiles[annotation])

def loadDict(filepath):
    unpack = np.load(filepath)
    annoDict = unpack.item()
    return annoDict

def getSignals(enhSeq,*dictlist):
    twoSigs = dict() 

    for seq in enhSeq:
        # [thisdict[seq] for thisdict in *dictlist]
        reqSignals = [(thisdict.get(seq, None)) for thisdict in dictlist]
        twoSigs[seq] = reqSignals
    return twoSigs        

def plot_corr(sig1, sig2, num): 
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
    plt.savefig(saveFigDir+str(num)+".png")

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
        # if ktrck == 100: 
        #     break
        # ("plotting seq {} out of {} for {} with annotation {}".format(ktrck, len(annotDict),element_NAME,annotation_NAME))
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

def NaNToZero(nparray):
    nparray[np.isnan(nparray)] = 0 
    return nparray

def plot_bxplt_corr(noEnh_data, Enh_data, funcName, plttitle, *saveFlag):
    plt.figure()
    plt.boxplot([noEnh_data, Enh_data], sym="",labels=["No Enhancer Overlap", "Enhancer Overlap"],widths=0.5)
    for i in np.arange(len([noEnh_data, Enh_data])): 
        y = [noEnh_data, Enh_data][i]
        x = np.random.normal(1+i, 0.04, size=len(y))
        plt.plot(x,y,'r+', alpha=0.4, markersize=2)
    plt.title(plttitle)
    plt.ylabel(funcName)
    if saveFlag: 
        plt.savefig(saveFigDir+"MER20_bxplt_"+funcName+".eps")
    else: 
        plt.show()

def calc_phyloP_TF_assoc(sequences, func):
    store_dict = dict()
    store_dict_norm = dict()
    for i, thisSeq in enumerate(sequences):
        sig1 = phyloPArray[thisSeq]
        sig2 = NaNToZero(TFArray[thisSeq])
        sig1m = sig1[~np.isnan(sig1)]
        sig2m = sig2[~np.isnan(sig1)]

        if (sig1m.size != 0) and (sig2m.size != 0): 
            store_dict[thisSeq] = func(sig1m,sig2m)[0] 
            if max(sig2m) == 0: 
                store_dict_norm[thisSeq] = func(sig1m,sig2m)[0] 
            else:     
                store_dict_norm[thisSeq] = func(sig1m,sig2m/max(sig2m))[0]

    return [ np.asarray(list(store_dict.values())).ravel(), np.asarray(list(store_dict_norm.values())).ravel()]


#-------
# main
#-------

### load annotations 
phastArray = loadDict(getAnnoPath("phastCon"))
TFArray = loadDict(getAnnoPath("TF"))
enhArray = loadDict(getAnnoPath("enhan"))
phyloPArray = loadDict(getAnnoPath("phyloP"))

#set up a set of seq w/ and w/o enhc overlap
enhSeq = set(enhArray.keys())
phyloP_TF_Seq = set(phyloPArray.keys()).intersection(TFArray.keys())
noEnhc_phyloP_TF_Seq = phyloP_TF_Seq - enhSeq
enhc_phyloP_TF_Seq = phyloP_TF_Seq.intersection(enhSeq)  

# =============  RAW DATA PLOTS =============
### create raw data plots: 
elem = "MER20"
annoList =  ["phastCon","TF", "enhan","phyloP"]
run_plot_raw(annoList)

# =============  compare TE w/ and w/o enhc overlap =============
#initialize 

[corr_noEnhc, norm_corr_noEnhc] = calc_phyloP_TF_assoc(noEnhc_phyloP_TF_Seq, np.correlate)
[corr_Enhc, norm_corr_Enhc] = calc_phyloP_TF_assoc(enhc_phyloP_TF_Seq, np.correlate)

[corrcoef_noEnhc, norm_corrcoef_noEnhc] = calc_phyloP_TF_assoc(noEnhc_phyloP_TF_Seq, stats.pearsonr)
[corrcoef_Enhc, norm_corrcoef_Enhc] = calc_phyloP_TF_assoc(enhc_phyloP_TF_Seq, stats.pearsonr)


# =============  BOXPLOT of Cross-Corr of PhyloP and TF Motif Score =============
plot_bxplt_corr(corr_noEnhc, corr_Enhc,"CrossCorr","Distribution of CrossCorrelation Values \nBetween (NOT Normalized) TFmotif Score and PhyloP for MER20")
plot_bxplt_corr(norm_corr_noEnhc, norm_corr_Enhc,"CrossCorr","Distribution of CrossCorrelation Values \nBetween normalized TFmotif Score and PhyloP for MER20")

plot_bxplt_corr(corrcoef_noEnhc, corrcoef_Enhc,"PearsonCorrCoef","Distribution of Pearson Correlation Coeff Values \nBetween (NOT Normalized) TFmotif Score and PhyloP for MER20",True)
plot_bxplt_corr(norm_corrcoef_noEnhc, norm_corrcoef_Enhc,"PearsonCorrCoef_(norm)","Distribution of Pearson Correlation Coeff Values \nBetween normalized TFmotif Score and PhyloP for MER20", True)

# plt.figure()
# plt.boxplot([corr_noEnhc, corr_Enhc], sym="",labels=["No Enhancer Overlap", "Enhancer Overlap"],widths=0.5)
# for i in np.arange(len([corr_noEnhc, corr_Enhc])): 
#     y = [corr_noEnhc, corr_Enhc][i]
#     x = np.random.normal(1+i, 0.04, size=len(y))
#     plt.plot(x,y,'r+', alpha=0.4, markersize=2)
# plt.title("Distribution of CrossCorrelation Values \nBetween (NOT Normalized) TFmotif Score and PhyloP for MER20")
# plt.ylabel("CrossCorrelation(no lag)")
# plt.show()

# =============  Violin Plots of Cross-Corr of PhyloP and TF Motif Score =============
plt.figure()
plt.violinplot([corr_noEnhc_norm, corr_Enhc_norm], [1,2], points = 100, widths=0.5, showextrema=True, showmeans=True, showmedians=True)
plt.show()

# =============  SCATTER PLOTS OF PhyloP AND TF Motif Score =============
#scatter plot of phyloP and TFmotif score per Base for all Enhancer Overlapping TEs
plt.figure()
for n,seq in enumerate(enhc_phyloP_TF_Seq): 
    phdata = phyloPArray[seq]
    TFdata = NaNToZero(TFArray[seq])
    mask_enhanc = enhArray[seq]
    #remove NaN indices from phData, and remove same indiceis from TFData
    mod_phdata = phdata[~np.isnan(phdata)]
    mod_TFdata = TFdata[~np.isnan(phdata)]
    mod_mask_enhanc = mask_enhanc[~np.isnan(phdata)]

    if n ==1: 
        labelOV= "Bases Overlapping Enhancer"
        labelV= "Bases NOT overlapping Enhancer"
    else:
        labelOV= ""
        labelV= ""

    plt.plot(mod_TFdata[mod_mask_enhanc==1], mod_phdata[mod_mask_enhanc==1],"r+",markersize=4, alpha=0.4, label=labelOV)
    plt.plot(mod_TFdata[np.isnan(mod_mask_enhanc)], mod_phdata[np.isnan(mod_mask_enhanc)],"k.",markersize=4, alpha=0.4, label=labelV)
    
plt.legend(loc='best',handletextpad=0.3,framealpha=0.4)
plt.title('Conservation vs. TFmotif Score PER BASE for MER20s \nWith Any Overlap With Enhancer')
plt.xlabel('TFmotif Score')
plt.ylabel('PhyloP Score')
plt.show()

#scatter plot of phyloP and TFmotif score per Base for MER20 Without Enhancer Overlap.
plt.figure()
for n,seq in enumerate(noEnhc_phyloP_TF_Seq): 
    phdata = phyloPArray[seq]
    TFdata = NaNToZero(TFArray[seq])
    #remove NaN indices from phData, and remove same indiceis from TFData
    mod_phdata = phdata[~np.isnan(phdata)]
    mod_TFdata = TFdata[~np.isnan(phdata)]
        
    if n==1:
        plt.plot(mod_TFdata, mod_phdata,"k.",markersize=4, alpha=0.4, label="one base in a given MER20")
    else: 
        plt.plot(mod_TFdata, mod_phdata,"k.",markersize=4, alpha=0.4, label="")
    
plt.legend(loc='best',handletextpad=0.3,framealpha=0.4)
plt.title('Conservation vs. TFmotif Score PER BASE for MER20s \nWITHOUT Any Overlap With Enhancer')
plt.xlabel('TFmotif Score')
plt.ylabel('PhyloP Score')
plt.show()


# # =============  Boxplot of Distribution of PhyloP by TF motif(Boolean) and Enhancer Overlap  =============
# store_phdata_NoTF_YesEnhc = np.array(0)
# store_phdata_YesTF_YesEnhc = np.array(0)
# store_phdata_YesTF_NoEnhc = np.array(0)
# store_phdata_NoTF_NoEnhc = np.array(0)

# for n,seq in enumerate(enhc_phyloP_TF_Seq): 
#     phdata = phyloPArray[seq]
#     TFdata = NaNToZero(TFArray[seq])
#     mask_enhanc = enhArray[seq]
#     #remove NaN indices from phData, and remove same indiceis from TFData
#     mod_phdata = phdata[~np.isnan(phdata)]
#     mod_TFdata = TFdata[~np.isnan(phdata)]
#     mod_mask_enhanc = mask_enhanc[~np.isnan(phdata)]

#     mod_TFdata[mod_mask_enhanc==1]
#     phdata_NoTF_YesEnhc = mod_phdata[mod_mask_enhanc==1 and mod_TFdata ==0 ] #bases with NO TF motif hit & Enhancer Overlap
#     phdata_YesTF_YesEnhc = mod_phdata[mod_mask_enhanc==1 and mod_TFdata !=0 ] 
#     phdata_YesTF_NoEnhc = mod_phdata[np.isnan(mod_mask_enhanc) and mod_TFdata !=0 ] 
#     phdata_NoTF_NoEnhc = mod_phdata[np.isnan(mod_mask_enhanc) and mod_TFdata ==0 ] 

#     store_phdata_NoTF_YesEnhc = np.append(store_phdata_NoTF_YesEnhc, phdata_NoTF_YesEnhc)
#     store_phdata_YesTF_YesEnhc = np.append(store_phdata_YesTF_YesEnhc, phdata_YesTF_YesEnhc)
#     store_phdata_YesTF_NoEnhc = np.append(store_phdata_YesTF_NoEnhc, phdata_YesTF_NoEnhc)
#     store_phdata_NoTF_NoEnhc = np.append(store_phdata_NoTF_NoEnhc, phdata_NoTF_NoEnhc)








# plt.title('MER20 with any Enhancer Overlap: Distribution of PhyloP Values per Base (Collapsed Across MER20)')
# plt.ylabel('PhyloP Score')
# plt.xlabel('...')
# plt.legend()



'''
plt.figure()
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
for i in forscatter_Enhc: 
    plt.plot(i[1], i[0], ".", label="MER20")

plt.title("PhyloP and TF motif Score per Base for MER20s Overlapping a Enhancer")
plt.xlabel("TFmotif Score")
plt.ylabel("PhyloP Score")

plt.show()
'''



#plot enhancer activity and TF activity 
# twoSigDict = getSignals(enhSeq, phyloPArray, TFArray)
# for k,v in twoSigDict.items():
#     print()
#     print(k) 
#     print(v)

#     plt.figure()
#     sig1 = v[0]
#     sig2 = v[1]
#     sig1[np.isnan(sig1)] = 0
#     sig2[np.isnan(sig2)] = 0

    # plt.plot(sig1,'r',label='Enhancer Activity')
    # plt.plot(sig2,'b',label='TF Activity')
    # plt.title(k)
    # plt.legend()
    # plt.savefig("/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/MER20_enh_TF_singleSeq/"+k+"_2.eps")


# sig1, sig2 = twoSigDict[ 'chr6:110361036-110361252']
# sig1[np.isnan(sig1)] = 0 
# sig2[np.isnan(sig2)] = 0 


### compare two signals 

# plot_raw(phyloPArray, "MER20", "phyloP")

# to do 
#   same plots as above, but split into w/ enhc overlap and without 