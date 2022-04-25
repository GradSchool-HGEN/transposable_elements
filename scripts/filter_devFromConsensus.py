#!/bin/python
# This script will compare # of bases deviation between aligned consensus and genomic hit for a TE; then filter based on percent difference or number of different bases. 
#
# Inputs: 
#   python filter_devFromConsensus.py <element_name, default=MER20> -n <#basesThreshold> -p <PercentBasesThreshold> 
# Outputs: 
#   filter_elem funciton returns, a dictionary object with sequence:filterBoolean, where 1 = under threshold
#    
# Dependencies: 
#   > TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# 
# Abin Abraham
# created on: 2017-12-18 17:18:00   

import numpy as np 
import pandas as pd 
import subprocess 
import matplotlib.pyplot as plt 
import argparse
import datetime
import sys
import os
#-------
# functions
#-------    
def filter_elem(element, perBase_thres, numBase_thres):
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    filterForElement = subprocess.check_output(['grep', '-w', element, TEbedfile_PATH], universal_newlines=True)
    out_data =[x.split('\t') for x in filterForElement.split('\n')]
    col_name = ["chr", "start", "end", "TE", "model_start","model_end","model_length"]
    df = pd.DataFrame(data=out_data[:-1], columns=col_name)
    df = df.apply(pd.to_numeric, errors='ignore')

    df['elemSeq_name'] = df['chr'] +":"+ df['start'].astype(str) +"-"+ df['end'].astype(str)
    df['genomicNumBases'] = df.end - df.start 
    df['alignConsensNumBases'] = df.model_end - df.model_start
    df['deviationFromConsensus'] = abs(df.genomicNumBases-df.alignConsensNumBases)
    df['per_deviationFromConsensus'] = abs(df.genomicNumBases-df.alignConsensNumBases)/df.model_length

    #=========== filter ===========    
    # by number of bases 
    df['filter_byNumber'] = 0
    df.loc[df.deviationFromConsensus < numBase_thres, 'filter_byNumber'] = 1
    ndf = df[df.deviationFromConsensus < numBase_thres].copy() 
    ndf["genomicNumBases"] = ndf.end - ndf.start 
    ndf['alignConsensNumBases'] = ndf.model_end - ndf.model_start
    ndf['deviationFromConsensus'] = abs(ndf.genomicNumBases-ndf.alignConsensNumBases)
    ndf['per_deviationFromConsensus'] = abs(ndf.genomicNumBases-ndf.alignConsensNumBases)/ndf.model_length
    
    # by bases deviation compared to total length (percent)
    df['filter_byPercent'] = 0
    df.loc[df.per_deviationFromConsensus < perBase_thres, 'filter_byPercent'] = 1
    pdf = df[df.per_deviationFromConsensus < perBase_thres].copy() 
    pdf["genomicNumBases"] = pdf.end - pdf.start 
    pdf['alignConsensNumBases'] = pdf.model_end - pdf.model_start
    pdf['deviationFromConsensus'] = abs(pdf.genomicNumBases-pdf.alignConsensNumBases)
    pdf['per_deviationFromConsensus'] = abs(pdf.genomicNumBases-pdf.alignConsensNumBases)/pdf.model_length

    #=========== create dictionary to return ===========    
    filteredByNum_dict = dict(zip(list(df.elemSeq_name),list(df.filter_byNumber)))
    filteredByPer_dict = dict(zip(list(df.elemSeq_name),list(df.filter_byPercent)))


    return filteredByNum_dict, filteredByPer_dict, df

def QCplots(filteredDataFrame, perThreshold, numThreshold, element): 
    fdf = filteredDataFrame
    #=========== scatter plots ===========
    plt.figure(0)
    plt.scatter(fdf.loc[fdf.filter_byPercent == 1, 'genomicNumBases'], fdf.loc[fdf.filter_byPercent == 1, 'alignConsensNumBases'], color='r',  alpha=0.6,  marker=".", s =2, label=element+" <"+str(perThreshold))
    plt.scatter(fdf.loc[fdf.filter_byPercent == 0, 'genomicNumBases'], fdf.loc[fdf.filter_byPercent == 0, 'alignConsensNumBases'], color='b',  alpha=0.6,  marker=".", s =2, label=element+" >="+str(perThreshold))
    plt.title("Filtering "+element+ " by Percent: \n  #Bases Diff Between Genomic Hit & Consensus/Length of TE")
    plt.xlabel("Number of Bases in Genomic "+element) 
    plt.ylabel("Number of Bases in Aligned Consensus "+element) 
    plt.legend(loc=2)
    plt.plot([0, max(max(fdf.genomicNumBases), max(fdf.alignConsensNumBases))], [0, max(max(fdf.genomicNumBases), max(fdf.alignConsensNumBases))], '--k') #diagonal line 
    plt.text(0.4, 0.2, "Percent Threshold: " + str(perThreshold))
    plt.show()

    plt.figure(1)
    plt.scatter(fdf.loc[fdf.filter_byNumber == 1, 'genomicNumBases'], fdf.loc[fdf.filter_byNumber == 1, 'alignConsensNumBases'], color='r',  alpha=0.6,  marker=".", s =2, label=element+" <"+str(numThreshold))
    plt.scatter(fdf.loc[fdf.filter_byNumber == 0, 'genomicNumBases'], fdf.loc[fdf.filter_byNumber == 0, 'alignConsensNumBases'], color='b',  alpha=0.6,  marker=".", s =2, label=element+" >="+str(numThreshold))
    plt.title("Filtering "+element+ " by Number of Bases: \n  #Bases Diff Between Genomic Hit & Consensus")
    plt.xlabel("Number of Bases in Genomic "+element) 
    plt.ylabel("Number of Bases in Aligned Consensus "+element) 
    plt.legend(loc=2)
    plt.plot([0, max(max(fdf.genomicNumBases), max(fdf.alignConsensNumBases))], [0, max(max(fdf.genomicNumBases), max(fdf.alignConsensNumBases))], '--k') #diagonal line 
    plt.text(0.4, 0.2, "Number of Bases Threshold: " + str(numThreshold))
    plt.show()

    #=========== cumulative distribution plots ===========
    # percent 
    perThresholds = np.linspace(0,1,200)
    totalNum = fdf.shape[0]
    yCounts = [fdf[fdf.per_deviationFromConsensus < perThresh].shape[0] for perThresh in perThresholds]
    yPercent = [x/totalNum for x in yCounts]

    plt0, ax = plt.subplots()
    ax.plot(perThresholds, yPercent)
    title0 = 'QC: Number of ' +element+ ' as function of Filter Threshold (Percent)'
    xlab0 = '#Bases Diff Between Genomic Hit & Consensus/Length of TE'
    ylab0 = '# of TE below threshold/total num of TE \n(percent)  '
    ax.set(xlabel=xlab0, ylabel=ylab0, title=title0)
    ax.grid()
    ax.text(0.4, 0.2, 'Total Number of '+element+" before \nfiltering: " +str(totalNum))
    plt0.show()
    
    #count 
    countThresholds = np.linspace(0,50,500)
    cCounts = [fdf[fdf.deviationFromConsensus < countThresh].shape[0] for countThresh in countThresholds]
    cPercent = [x/totalNum for x in cCounts]

    plt2, ax2 = plt.subplots()
    ax2.plot(countThresholds, cPercent)
    title2 = 'QC: Number of ' +element+ ' as function of Filter Threshold (#Bases)'
    xlab2 = '#Bases Diff Between Genomic Hit & Consensus'
    ylab2 = '# of TE below threshold/total num of TE \n(percent)  '
    ax2.set(xlabel=xlab2, ylabel=ylab2, title=title2)
    ax2.grid()
    ax2.text(10, 0.2, 'Total Number of '+element+" before \nfiltering: " +str(totalNum))
    plt2.show()
    


#-------
# main
#-------
def main(argv):
    # print header
    print('{:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))

    # -----------
    # Arguments 
    # ----------- 
    parser = argparse.ArgumentParser(description="a short description of what this does for help text ~")
    parser.add_argument("element", default='MER20',
                    action='store', type=str,
                    help="Select a transposable element") 

    parser.add_argument("-n", type=int, action='store',
                    default=10, dest='nThres', 
                    help='Number of Bases Different From Consensus')  
    parser.add_argument("-p", type=float, action='store',
                    default=0.1, dest='pThres', 
                    help='Percent of Bases Different From Consensus as funciton of TE length')

      
    argsIN = parser.parse_args()
    element= argsIN.element
    numBase_thres= argsIN.nThres
    perBase_thres= argsIN.pThres
    

    [filteredByNum_dict, filteredByPer_dict, fdf] = filter_elem(element, perBase_thres, numBase_thres)
    QCplots(fdf, perBase_thres, numBase_thres, element)
    

#-------
# Main function run 
#-------

if __name__ == "__main__":
    main(sys.argv[1:])
    #create output folder and put contents there 

