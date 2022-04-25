#!/bin/python
# This script will 
#
#
#
# Abin Abraham
# created on: ...⇧+⌘+I

import pandas as pd 
import numpy as np 
import os 
import glob 
import re
import matplotlib.pyplot as plt
from datetime import date 



saveFigDir = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/figs/"

# =============  functions =============
def plot_TFraw(annotDict, element_NAME, annotation_NAME, saveFigFlag=False):
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
        ("plotting seq {} out of {} for {} with annotation {}".format(ktrck, len(annotDict),element_NAME,annotation_NAME))
        ktrck += 1
        v = v.ravel()
        v[np.isnan(v)] = 0
        allArray = np.vstack((allArray,v))
        plt.plot(np.arange(consLen), v.ravel(), 'b-', linewidth = 0.7, alpha=0.1)
        if ~all(np.isnan(v)): counter +=1
    plt.grid()
    plt.title('Genomic Instances of {} annotated with  {} (n={})'.format(element_NAME, annotation_NAME, str(counter)), weight='bold')
    plt.ylabel(annotation_NAME)
    [x1,x2,y1,y2] = plt.axis()
    # frame1 = plt.gca()
    # frame1.axes.xaxis.set_ticklabels([])
    
  
    #plot sum of TF score at each base
    allArray = allArray[1:] #remove initializing zeros
    meanConsArray = np.nansum(allArray, axis=0)
    # stdConsArray = np.nanstd(allArray, axis=0)
    plt.subplot(3, 1, 2) #MEAN and STD 
    plt.plot(np.arange(consLen), meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    # plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.title("Sum of TF Scores Per Base", weight='bold')
    plt.ylabel(annotation_NAME)
    # frame1 = plt.gca()
    # frame1.axes.xaxis.set_ticklabels([])

    #plot average TF score at each base
    plt.subplot(3,1,3)
    meanTFscore = np.mean(allArray, axis=0)
    stdConsArray = np.nanstd(allArray, axis=0)
    plt.plot(np.arange(consLen), meanTFscore, 'r-', linewidth = 1.2, alpha=1)
    
    plt.plot(np.arange(consLen), meanTFscore + 2*stdConsArray, 'k:', linewidth = 0.7, alpha=1)
    plt.plot(np.arange(consLen), meanTFscore - 2*stdConsArray, 'k:', linewidth = 0.7, alpha=1)


    plt.grid()
    plt.title("Mean of TF Hits Per Base", weight='bold')
    plt.ylabel(annotation_NAME)


    pltname = "{}_{}_rawHMManno_{}.eps".format(element_NAME, annotation_NAME, date.today())
    if saveFigFlag: 
        plt.savefig(os.path.join(saveFigDir,pltname))
    else: 
        plt.show()



# =============  main =============



files_to_load = glob.glob("/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/mapped_TF/MER20_*_mappedData_mult_2018-02-26.npy")


store_mapping_dict = {}
for this_file in files_to_load: 
    this_tf = re.search('_[A-Z0-9]*_', this_file).group()[1:-1]
    print(this_tf)
    this_tf_mapping_dict = np.load(this_file)
    store_mapping_dict[this_tf] = this_tf_mapping_dict
    
    plot_TFraw(store_mapping_dict[this_tf].item(), 'MER20', this_tf, saveFigFlag=True)