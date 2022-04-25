#!/bin/python
# This script will take an annotation for a given TE, and map it onto the consensus sequence. 
#
#  Path Dependencies
#       1) TEbedfile_PATH
#  File Paths Passed to functions: 
#       1) consensusFastaPath (only path, no file name)
#       2) MapDictPath (only path, no file name)
#       3) annotationfile_PATH (file name incl)
#       4) intersectBedFile_PATH (file name incl)
#  Expected Names for Input Files 
#       1) consensusSeqfile - consensus sequence: "consensus_<element>.fa"
#       2) dictname - mapping dictionary: "all<element>_mappedDict.pi"
# 
# 
# 
#       TO DO: 
#           - remove dependence on unecessary datafram columns; convert script to accept a list of sequences of TE
# 
# Abin Abraham
# created on: 2017-12-27 14:58:11
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
def getConsLength(element, consensusFastaPath):
    consensusSeqfile = consensusFastaPath+"consensus_"+element+".fa"

    if os.path.isfile(consensusSeqfile) == False: 
        print("Consensus sequence fasta not found. File must be named consensus_"+element+".fa")
        print("File should be located at:"+consensusSeqfile )
        return
    else: 
        consensusSeq = SeqIO.read(consensusSeqfile,'fasta')
    
    return len(consensusSeq.seq) 

def getDictMapping(element, MapDictPath):
    dictname = "all"+element+"_mappedDict.pi"

    if os.path.isfile(MapDictPath+dictname) == False: 
        print("Mapping to consensus dictionary file not found.\n")
        print("File does not exist: "+MapDictPath+dictname )
    else: 
        allElemMapToCons_dict = pickle.load(open(MapDictPath+dictname, "rb" ))
    
    return allElemMapToCons_dict

def checkInputs(annotationfile_PATHc):
    numOfColumns = subprocess.check_output(['awk', '{print NF; exit}', annotationfile_PATHc], universal_newlines=True).splitlines()
    cmd = str("sort -k1,1 -k2,2n " + annotationfile_PATHc + " -o " + annotationfile_PATHc).split()
    subprocess.call(cmd)

    fileModifiedFLAG = False
    if int(numOfColumns[0]) > 4:
        fileModifiedFLAG = True
        print("\nAnnotation file has >4 columns. a new annotation file will be created and extra columns will be deleted.\n")
        cmd1 = 'cut -f1,2,3,4 '+annotationfile_PATHc
        cmd1 = cmd1.split()
        f = open("formated_annotationfile.tsv", "w")
        subprocess.call(cmd1, stdout=f)
        f.close()

    return fileModifiedFLAG

def runBedtoolsIntersect(element, annotation, intersectBedFile_PATH, annotationfile_PATH):
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    
    #check if intersect-file already exists, if not create it 
    if os.path.isfile(intersectBedFile_PATH) == False:
        print("\nBedtools intersect file not found in "+intersectBedFile_PATH+"\n")  
        raise SystemExit("Exited Script because intersect file was not found.")
        # print("File must be named 'intersect-output_<element>_<annotation>' \n")
        # print("Running bedtools intersect on annotation file:{}".format(annotationfile_PATH))
        # print('\n')
        # #check columns, modify to have only 4 columns 
        # fileFLAG = checkInputs(annotationfile_PATH) 
        # if fileFLAG:
        #     filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)
        #     intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', "formated_annotationfile.tsv", '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)  
        # else:
        #     filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)
        #     intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotationfile_PATH, '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, universal_newlines=True)
        
        # #write output to file
        # with open(intersectBedFile_PATH,'w') as f:
        #     f.write(intersectOutput)

        # print("File Saved in: "+intersectBedFile_PATH+"\n")
    else: 
        with open(intersectBedFile_PATH,'r') as f:
            intersectOutput = f.read()
        print("\nBedtools intersect file loaded from "+intersectBedFile_PATH+"\n") 

    return intersectOutput
    
def annotateConsensus(element, annotation, intersectBedFile_PATH, MappingDict_PATH, consensusFasta_PATH, annotationfile_PATH):
    # =============  summary of run =============
    print("---------------------------------------------------------------")
    print("Running annotateConsensus with the following parameters: \n" \
           "element = " +element+ "\n" \
           "annotation = " +annotation+ "\n" \
           "intersectBedFile_PATH = " +intersectBedFile_PATH+ "\n" \
           "MappingDict_PATH = " +MappingDict_PATH+ "\n" \
           "consensusFasta_PATH = " +consensusFasta_PATH+ "\n" \
           "annotationfile_PATH = " +annotationfile_PATH+ "\n" ) 
            
    print("---------------------------------------------------------------")

    # =============  Bedtools Intersect on DFam TE locations & Annotation Locations =============
    intersectOutput = runBedtoolsIntersect(element, annotation, intersectBedFile_PATH, annotationfile_PATH)
   
    #=========== load data into data frame ===========
    out_data =[x.split('\t') for x in intersectOutput.split('\n')]
    col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "TE_genomic_chr","TE_genomic_start","TE_genomic_end", "TE_model","model_start","model_end", "model_length" ] 
    df = pd.DataFrame(data=out_data[:-1], columns=col_name)
    print("Dataframe is loaded.")

    #=========== norm index of intesection to genomic instance of element ===========
    df["TE_genomicHit_name"] = df["TE_genomic_chr"] + ":" + df["TE_genomic_start"].map(str) + "-" + df["TE_genomic_end"].map(str)
    df = df.apply(pd.to_numeric, errors='ignore')
    df["normToGenomicTE_start"] = df.intersect_start - df.TE_genomic_start 
    df["normToGenomicTE_end"] = df.intersect_end - df.TE_genomic_start 
    uniq_genomicHits = df.TE_genomicHit_name.unique()
    print("Completed grouping of data frame.")
    
    # =============  Organize Data for Output =============
    consensusLength = getConsLength(element, consensusFasta_PATH) # might have to link this to the profile HMM concensus
    consArray = np.zeros(consensusLength)
    tmp_consArray = np.zeros(consensusLength)
    mask_consArray = np.zeros(consensusLength)
    print("Done with organizing data.")

    # =============  Set Up Mapping Dictionary  =============
    # reformat dictionary key 
    allMapDict = getDictMapping(element,MappingDict_PATH)
    renamedMapDict = dict()
    for k in allMapDict.keys():
        newk = k.split(":")[2] +":"+ k.split(":")[3][:-3]
        renamedMapDict[newk] = allMapDict[k]
    print("Dictionary Loaded")

    # =============  Annotate Consensus For Each TE instance =============
    annotatedDict = dict()
    ##tempCounter = 0
    for ind,thisElem in enumerate(uniq_genomicHits):
        ##t0 = time.time()
        #filter for one TE
        oneElem_mapped = df.loc[df['TE_genomicHit_name'] == thisElem,["normToGenomicTE_start", "normToGenomicTE_end","annotation_score"]]
        print("currently on: {}, total n = {}".format(str(ind), str(len(uniq_genomicHits))))
        if renamedMapDict.get(thisElem):
            #get mapping to consensus 
            thisMapDict = renamedMapDict.get(thisElem)
            
            #for this TE, loop throgh indices with annotation data
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
            # print(thisElem)
        else:
            tmp_consArray[mask_consArray==0] = np.nan
            annotatedDict[thisElem] = tmp_consArray
        
        ##print("one outer loop took:"+str(time.time()-t0))
        ##tempCounter += 1 
        ##if tempCounter == 3: break 
    
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
    MappingDict_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-output_12_27_2017/"
    consensusFasta_PATH =  "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-input/"
    annotationfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

    #output
    outputDir = "/dors/capra_lab/abraha1/projects/transposable_elements/results/hmm_align/"

    annotDict= annotateConsensus(element, annotation, intersectBedFile_PATH, MappingDict_PATH, consensusFasta_PATH, annotationfile_PATH)
    plot_MapToConsensus(annotDict, element, annotation)
    return annotDict

    

#-------
# Main function run 
#-------


if __name__ == "__main__":
    annotDict = test(sys.argv[1:])
