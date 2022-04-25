#!/bin/python
# This script will take an annotation for a given TE, and map it onto the consensus sequence. 
#   > inputs: 
#       - element : TE Element of interest
#       - annotation : annotation of interest 
#       - intersectBedFile_PATH : path for output file from running bedtools intersect on TE of interest with annotation 
#       - MappingDictFile_PATH : path for dictionary file used to map a TE to its consensus 
#       - consensusFastaFile_PATH : path for consensus sequence for TE of interest 
#           
#   > depends: 
#       - /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy
#
#   > outputs: 
#       - outputDir : output directory to store annotation 
#
#   > Expected Names for Input Files 
#       - consensusSeqfile - consensus sequence: "consensus_<element>.fa"
#       - dictname - mapping dictionary: "all<element>_mappedDict.pi"
# Abin Abraham
# created on: 2017-12-26 00:12:31
# modified  : 2018-01-08 12:45:29


from Bio import SeqIO
import pandas as pd
import numpy as np
import subprocess
import datetime
import pickle
import sys
import os
import time


# -----------
# Functions
# -----------
def getConsLength(element, consfile_PATH):
    # consensusSeqfile = consensusFastaPath+"consensus_"+element+".fa"

    if os.path.isfile(consfile_PATH) == False: 
        # print("Consensus sequence fasta not found. File must be named consensus_"+element+".fa")
        print("Consensus sequences not found. Looked for: ".format(consfile_PATH))
        return
    else: 
        consensusSeq = SeqIO.read(consfile_PATH,'fasta')
    
    return len(consensusSeq.seq) 

def getDictMapping(element, dict_PATH):
    # dictname = "all"+element+"_mappedDict.pi"

    if os.path.isfile(dict_PATH) == False: 
        print("Mapping to consensus dictionary file not found.\n")
        print("File does not exist: {}".format(dict_PATH))
    else: 
        allElemMapToCons_dict = pickle.load(open(dict_PATH, "rb" ))
    
    return allElemMapToCons_dict

def getBedtoolsIntersect(element, annotation, intersectBedFile_PATH ):
    
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    
    #check if intersect-file already exists, if not exit
    if os.path.isfile(intersectBedFile_PATH) == False:
        print("\n ***Bedtools intersect file not found in "+intersectBedFile_PATH+"***\n")  
        raise SystemExit("Exited Script because intersect file was not found.")
    else: 
        with open(intersectBedFile_PATH,'r') as f:
            intersectOutput = f.read()
        print("\nBedtools intersect file loaded from "+intersectBedFile_PATH+"\n") 

    return intersectOutput
    
def annotateConsensus(element, annotation, intersectBedFile_PATH, MappingDictFile_PATH, consensusFastaFile_PATH):
    # =============  summary of run =============
    start = time.time()
    print("---------------------------------------------------------------")
    print("Running annotateConsensus with the following parameters: \n" \
           "element = " +element+ "\n" \
           "annotation = " +annotation+ "\n" \
           "intersectBedFile_PATH = " +intersectBedFile_PATH+ "\n" \
           "MappingDictFile_PATH = " +MappingDictFile_PATH+ "\n" \
           "consensusFastaFile_PATH = " +consensusFastaFile_PATH+ "\n" ) 
    print("---------------------------------------------------------------")

    #=========== load intersect file =============
    intersectOutput = getBedtoolsIntersect(element, annotation, intersectBedFile_PATH)
   
    #=========== load data into data frame ===========
    out_data =[x.split('\t') for x in intersectOutput.split('\n')]
    col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "TE_genomic_chr","TE_genomic_start","TE_genomic_end", "TE_model","model_start","model_end", "model_length"] 
    df = pd.DataFrame(data=out_data[:-1], columns=col_name)
    print("Dataframe is loaded.")

    #=========== norm index of intesection to genomic instance of element ===========
    df["TE_genomicHit_name"] = df["TE_genomic_chr"] + ":" + df["TE_genomic_start"].map(str) + "-" + df["TE_genomic_end"].map(str)
    # df = df.apply(pd.to_numeric, errors='ignore')
    df[["intersect_start", "intersect_end", "annotation_score", "TE_genomic_start", "TE_genomic_end", "model_start", "model_end", "model_length"]]= df[["intersect_start", "intersect_end", "annotation_score", "TE_genomic_start", "TE_genomic_end", "model_start", "model_end", "model_length"]].apply(pd.to_numeric, errors='ignore')
    df["normToGenomicTE_start"] = df.intersect_start - df.TE_genomic_start 
    df["normToGenomicTE_end"] = df.intersect_end - df.TE_genomic_start 
    uniq_genomicHits = df.TE_genomicHit_name.unique()
    print("Completed mapping annotation coordinates relative to TE start coordinate.")
    
    # =============  Organize Data for Output =============
    consensusLength = getConsLength(element, consensusFastaFile_PATH) # might have to link this to the profile HMM concensus
    consArray = np.zeros(consensusLength)
    tmp_consArray = np.zeros(consensusLength)
    mask_consArray = np.zeros(consensusLength)
    print("Done with organizing data.")

    # =============  Set Up Mapping Dictionary  =============
    # reformat dictionary key 
    allMapDict = getDictMapping(element,MappingDictFile_PATH)
    renamedMapDict = dict()
    for k in allMapDict.keys():
        newk = k.split(":")[2] +":"+ k.split(":")[3][:-3]
        renamedMapDict[newk] = allMapDict[k]
    print("Mapping Dictionary Loaded")

    # =============  Annotate Consensus For Each TE instance =============
    annotatedDict = dict()
    for ind,thisElem in enumerate(uniq_genomicHits):
        #filter for one TE
        oneElem_mapped = df.loc[df['TE_genomicHit_name'] == thisElem,["normToGenomicTE_start", "normToGenomicTE_end","annotation_score"]]
        print("currently on: {}, total n = {}".format(str(ind+1), str(len(uniq_genomicHits))))
        
        #get mapping to consensus 
        if renamedMapDict.get(thisElem):
            thisMapDict = renamedMapDict.get(thisElem)
            
            #for this TE, loop through each annotation value (one row per dataframe)
            for ind in range(len(oneElem_mapped)): 
                s,e,v = oneElem_mapped.iloc[ind]
                indexarray = np.arange(s,e)
                #get the indices that match to conensus seq 
                thisconsarrayind = [int(thisMapDict[x]) for x in indexarray if (thisMapDict.get((x)) and int(thisMapDict[x]) >=0)]
                #add annotation to those indices 
                tmp_consArray[thisconsarrayind] = (tmp_consArray[thisconsarrayind] + v)
                #update mask to keep count 
                mask_consArray[thisconsarrayind] += 1

            #divide by mask to get average
            tmp_consArray[mask_consArray>0] = tmp_consArray[mask_consArray>0]/mask_consArray[mask_consArray>0]
            #those indicies without any annotations are converted to nan
            tmp_consArray[mask_consArray==0] = np.nan          
            annotatedDict[thisElem] = tmp_consArray
            
            #reset temporary arrys 
            tmp_consArray = np.zeros(consensusLength)
            mask_consArray = np.zeros(consensusLength)
        else:
            tmp_consArray[mask_consArray==0] = np.nan
            annotatedDict[thisElem] = tmp_consArray
    end = time.time()
    print("Finished annotating consensus sequence in {} minutes.\n".format(round((end-start)/60),3))
    return annotatedDict

def plot_MapToConsensus(annotDict, element_NAME, annotation_NAME):

    import matplotlib.pyplot as plt
    consLen  = len(list(annotDict.values())[0].ravel())
    allArray = np.zeros((1,consLen))
    counter = 0 
    plt.figure()

    #plot raw data 
    plt.subplot(2,1,1)
    for k,v in annotDict.items(): 
        v = v.ravel()
        allArray = np.vstack((allArray,v))
        plt.plot(np.arange(consLen), v.ravel(), 'b-', linewidth = 0.7, alpha=0.3)
        if ~all(np.isnan(v)): counter +=1
    
    plt.grid()
    #plt.title('Genomic Instances of ' +'element_NAME'+ ' annotated \n with ' + "annotation_NAME"+"(n="+str(counter)+")" )
    plt.title('Genomic Instances of {} annotated \n with  {} (n={})'.format(element_NAME, annotation_NAME, str(counter)))
    plt.ylabel(annotation_NAME)
    [x1,x2,y1,y2] = plt.axis()
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])
    
    #plot mean and stdv
    allArray = allArray[1:] #remove initializing zeros
    meanConsArray = np.nanmean(allArray, axis=0)
    stdConsArray = np.nanstd(allArray, axis=0)
    plt.subplot(2, 1, 2) #MEAN and STD 
    plt.plot(np.arange(consLen), meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(np.arange(consLen), meanConsArray + 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.plot(np.arange(consLen), meanConsArray - 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.title('Mean (Red) and 2*SD (blue) of ' + annotation_NAME + ' annotation of ' + element_NAME)
    plt.xlabel('consensus bases (1start)')
    plt.ylabel(annotation_NAME)
    
    plt.show()

def test(argv):
    # print header
    # print('{:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))

    # -----------
    # Arguments 
    # ----------- 

    element = "MER20"
    annotation = "phastCon100test"

    # Require Inputs 
    intersectBedFile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/intersectBed-output/intersect-output_"+element+"_"+annotation
    MappingDictFile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-output_12_27_2017/all"+element+"_mappedDict.pi"
    consensusFastaFile_PATH =  "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-input/"+"consensus_"+element+".fa"

    #output
    outputDir = "/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/"

    annotDict = annotateConsensus(element, annotation, intersectBedFile_PATH, MappingDictFile_PATH, consensusFastaFile_PATH)
    #plot_MapToConsensus(annotDict, element, annotation)
    return annotDict

    

#-------
# Main function run 
#-------


if __name__ == "__main__":
    annotDict = test(sys.argv[1:])
