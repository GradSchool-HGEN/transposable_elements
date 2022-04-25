#!/bin/python
# This script will take an annotation for a given TE, and map it onto the consensus sequence. 
#
#
#
# Abin Abraham
# created on: 2017-12-27 14:58:11
from Bio import SeqIO
import pandas as pd
import numpy as np
import subprocess
import argparse
import datetime
import pickle
import sys
import os

# annotationfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

# -----------
# Functions
# -----------
def getConsLength(element):
    rootpath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/"
    input_dir = "align-input/"
    consensusSeqfile = rootpath+input_dir+"consensus_"+element+".fa"
    if os.path.isfile(consensusSeqfile) == False: 
        print("Consensus sequence fasta not found. File must be named consensus_"+element+".fa")
        print("File should be located at:"+consensusSeqfile )
        return
    else: 
        consensusSeq = SeqIO.read(consensusSeqfile,'fasta')
    
    return len(consensusSeq.seq) 

def getDictMapping(element):
    mappingDictFile = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-output_12_27_2017/"
    dictname = "all"+element+"_mappedDict.pi"
    allElemMapToCons_dict = pickle.load(open(mappingDictFile+dictname, "rb" ))
    return allElemMapToCons_dict

def checkInputs(annotationfile_PATHc):
    numOfColumns = subprocess.check_output(['awk', '{print NF; exit}', annotationfile_PATHc], universal_newlines=True).splitlines()
    cmd = str("sort -k1,1 -k2,2n " + annotationfile_PATHc + " -o " + annotationfile_PATHc).split()
    subprocess.call(cmd)

    fileModifiedFLAG = False
    if int(numOfColumns[0]) > 4:
        fileModifiedFLAG = True
        print("\nAnnotation file has >4 columns, extra columns will not be analyzed.\n")
        cmd1 = 'cut -f1,2,3,4 '+annotationfile_PATHc
        cmd1 = cmd1.split()
        f = open("formated_annotationfile.tsv", "w")
        subprocess.call(cmd1, stdout=f)
        f.close()

    return fileModifiedFLAG

def runBedtoolsIntersect(element, annotation, inputPath, annotationfile_PATH):
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    # inputPath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/"
    intersectBed_dir = "intersectBed-output/"
    intersectBed_filename = inputPath+intersectBed_dir+"intersect-output_"+element+"_"+annotation

    if os.path.isdir(inputPath+intersectBed_dir) == False:
        os.mkdir(inputPath+intersectBed_dir)
    
    if os.path.isfile(intersectBed_filename) == False:
        print("\nBedtools intersect file not found in "+intersectBed_filename+"\n")  
        print("File must be named 'intersect-output_<element>_<annotation>' \n")
        print("Running bedtools intersect....\n")

        fileFLAG = checkInputs(annotationfile_PATH) #check columns, modify to have only 4 columns 
        if fileFLAG:
            filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)
            intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', "formated_annotationfile.tsv", '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)  
        else:
            filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)
            intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotationfile_PATH, '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, universal_newlines=True)
    
        #write output to file
        with open(intersectBed_filename,'w') as f:
            f.write(intersectOutput)
        print("File Saved in: "+intersectBed_filename+"\n")
    else: 
        with open(intersectBed_filename,'r') as f:
            intersectOutput = f.read()

    return intersectOutput
    

def annotateConsensus(element, annotation, inputPath, annotationfile_PATH):
   
    intersectOutput = runBedtoolsIntersect(element, annotation, inputPath, annotationfile_PATH)
   
    #=========== load data into data frame ===========
    out_data =[x.split('\t') for x in intersectOutput.split('\n')]
    col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "TE_genomic_chr","TE_genomic_start","TE_genomic_end", "TE_model","model_start","model_end", "model_length" ] 
    df = pd.DataFrame(data=out_data[:-1], columns=col_name)

    #=========== norm index of intesection to genomic instance of element ===========
    df["TE_genomicHit_name"] = df["TE_genomic_chr"] + ":" + df["TE_genomic_start"].map(str) + "-" + df["TE_genomic_end"].map(str)
    df = df.apply(pd.to_numeric, errors='ignore')
    df["normToGenomicTE_start"] = df.intersect_start - df.TE_genomic_start 
    df["normToGenomicTE_end"] = df.intersect_end - df.TE_genomic_start 
    df = df.sort_values(by='TE_genomicHit_name')
    gb = df.groupby('TE_genomicHit_name')

    # Organize Data for Output
    consensusLength = getConsLength(element) # might have to link this to the profile HMM concensus
    consArray = np.zeros((1, consensusLength))
    tmp_consArray = np.zeros((1, consensusLength))
    mask_consArray = np.ones((1, consensusLength))

    allMapDict = getDictMapping(element)
    renamedMapDict = dict()
    for k in allMapDict.keys():
        newk = k.split(":")[2] +":"+ k.split(":")[3][:-3]
        renamedMapDict[newk] = allMapDict[k]

    nameList = [list(gb)[i][0] for i in range(len(gb))]

    for i in range(len(gb)):
        oneElem_mapped = list(gb)[i][1][["normToGenomicTE_start", "normToGenomicTE_end","annotation_score"]].values.astype(float)
        seqName = list(gb)[i][0]
        thisMapDict = renamedMapDict[seqName]
        for s, e, v in oneElem_mapped: # loop through each intersection for a given TE 
            indexarray = np.arange(s,e)
            thisconsarray = [int(thisMapDict[x]) for x in indexarray if (thisMapDict.get((x)) and int(thisMapDict[x]) >=0)]

            tmp_consArray[0, thisconsarray] = (tmp_consArray[0, thisconsarray] + v)
            mask_consArray[0, thisconsarray] += 1

        tmp_consArray = tmp_consArray/mask_consArray        
        consArray = np.vstack((consArray, tmp_consArray))
        tmp_consArray = np.zeros((1, consensusLength))
        mask_consArray = np.ones((1, consensusLength))

    consArray = consArray[1:, :] #remove the initializing row 

    return [consArray, nameList, consensusLength]

def plot_MapToConsensus(consensusLength, consArray, element_NAME, annotation_NAME, outputDir):

    import matplotlib.pyplot as plt
    #from cycler import cycler

    consBases = np.arange(0.0, consensusLength, 1, dtype=int)
    meanConsArray = np.mean(consArray, axis=0)
    stdConsArray = np.std(consArray, axis=0)

    plt.figure(figsize=(8, 5), dpi=160)
    
    plt.subplot(2, 1, 1) #RAW DATA l
    #plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
    for i in range(len(consArray)):
        plt.plot(consBases, consArray[i], 'b-', linewidth = 1, alpha=0.1)
    plt.grid()
    plt.title('Genomic Instances of ' + element_NAME + ' Annotated \n with ' + annotation_NAME )
    plt.ylabel(annotation_NAME)
    [x1,x2,y1,y2] = plt.axis()
    frame1 = plt.gca()
    frame1.axes.xaxis.set_ticklabels([])

    plt.subplot(2, 1, 2) #MEAN and STD 
    plt.plot(consBases, meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(consBases, meanConsArray + 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.plot(consBases, meanConsArray - 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.title('Mean (Red) and 2*SD (blue) of ' + annotation_NAME + ' annotation of ' + element_NAME)
    plt.xlabel('consensus bases (1start)')
    plt.ylabel(annotation_NAME)
    saveFigName = os.path.join(outputDir, element_NAME+"_"+annotation_NAME+"_MappedToConsensus_" + str(datetime.date.today())+ ".png")
    # plt.savefig(saveFigName)
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    plt.close()

    

def main(argv):
    # print header
    print('{:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))

    # -----------
    # Arguments 
    # ----------- 
    parser = argparse.ArgumentParser(description="a short description of what this does for help text ~")
    parser.add_argument("element", default='None',
                    action='store', type=str,
                    help="Select a transposable element") 
    parser.add_argument("annotationFilePath", 
                    action='store', type=str, 
                    help="Path or file name for annotation")
    
    # parser.add_argument("-s", type=bool, action='store_true',
    #                 default=False, dest='boolean_switch', 
    #                 help='save outputs to TE_annotate_output') 
    parser.add_argument("-o", type=str, action='store',
                    default="TE_annotate_output/", dest='outPath', 
                    help='path to store outputs; remember to append / to dirs')  
    parser.add_argument("-a", type=str, action='store',
                    default="annotations", dest='annoName', 
                    help='name of the annotation')  

    
    argsIN = parser.parse_args()
    element_NAME= argsIN.element
    annotationFile_NAME= argsIN.annotationFilePath
    annotationName= argsIN.annoName
    # saveOutputFLAG= argsIN.boolean_switch
    outputPath= argsIN.outPath


    #Set Up Output Directory 
    dir = os.path.dirname(outputPath)
    if not os.path.exists(dir): 
            print("\n NOTE: output directory does not exist; will default to TE_annotate_output in current directory")
            os.makedirs("TE_annotate_output")
            outputPath= os.path.join(os.getcwd(), 'TE_annotate_output') 


    [consArray, nameList, consensusLength] = annotateConsensus(element_NAME, annotationFile_NAME)
    plot_MapToConsensus(consensusLength, consArray, element_NAME, annotationName, outputPath)
    

#-------
# Main function run 
#-------


if __name__ == "__main__":
    main(sys.argv[1:])
    #create output folder and put contents there 

