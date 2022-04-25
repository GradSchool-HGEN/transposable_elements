#!/bin/python
# This script will overlay a given annotation for a transposable element or a list of elements over its consensus sequence
#
# Abin Abraham      
# created on: 2017-12-07 20:15:16
# 
# Depends on: 
#       --> TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# 
# Required Inputs: 
#       1) transposable element of interest 
#       2) annotation file e.g: annotationFile_NAME = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"
# 
# Output:  
# 
# To Do: 
#       > check if overlap occurs when mapping the annotation score for a given sequence 
#
# run TE_annotate.py "SVA_E" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

import pandas as pd
import numpy as np
import subprocess
import argparse
import datetime
import sys 
import os 

# -----------
# Functions 
# ----------- 
def checkInputs(annotationfile_PATH):
    numOfColumns = subprocess.check_output(['awk', '{print NF; exit}', annotationfile_PATH  ],  universal_newlines=True).splitlines()

    fileModifiedFLAG = False; 
    if int(numOfColumns[0]) > 4: 
        fileModifiedFLAG = True; 
        cmd = 'cut -f1,2,3,4'+annotationfile_PATH.split()
        f = open("formated_annotationfile.tsv", "w")
        numOfColumns = subprocess.call(cmd, stdout=f)
        f.close()

    return(fileModifiedFLAG)    

def annotateConsensus(element, annotationfile_PATH):
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    filterForElement = subprocess.Popen(['grep', element, TEbedfile_PATH ], stdout=subprocess.PIPE, universal_newlines=True)

    fileFLAG = checkInputs(annotationfile_PATH)
    if fileFLAG: 
        intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', "formated_annotationfile.tsv", '-b', 'stdin'], stdin = filterForElement.stdout,  universal_newlines=True)
    else:   
        intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotationfile_PATH, '-b', 'stdin'], stdin = filterForElement.stdout,  universal_newlines=True)

 
    out_data =[x.split('\t') for x in intersectOutput.split('\n')]
    col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "genomic_chr","genomic_start","genomic_end", "TE_model","model_start","model_end", "model_length" ] 

    df = pd.DataFrame( data = out_data[:-1], columns=col_name)
    df["genomicHit_name"] = df["genomic_chr"] + ":" + df["genomic_start"].map(str) + "-" + df["genomic_end"].map(str)
    df = df.apply(pd.to_numeric, errors = 'ignore')
    df["mapToConsensus_start"] = df.intersect_start - df.genomic_start + df.model_start
    df["mapToConsensus_end"] = df.intersect_end - 1 - df.genomic_start + df.model_start
    df = df.sort_values(by='genomicHit_name')
    gb = df.groupby('genomicHit_name')

    # Organize Data for Output 
    consensusLength = df.model_length.iloc[1]
    consArray = tmp_consArray = np.zeros((1, consensusLength))
    nameList = [list(gb)[i][0] for i in range(len(gb))]

    for i in range(len(gb)):
        spec = list(gb)[i][1][["mapToConsensus_start", "mapToConsensus_end","annotation_score"]].values.astype(float)
        # check for overlap 
        
        for s,e,v in spec: tmp_consArray[0, int(s):int(e)+1]=v
        consArray = np.vstack((consArray, tmp_consArray))

    consArray = consArray[1:, :]

    return([consArray, nameList, consensusLength])

def plot_MapToConsensus(consensusLength, consArray, element_NAME, annotation_NAME, outputDir):

    import matplotlib.pyplot as plt
    

    consBases = np.arange(0.0, consensusLength, 1, dtype=int)
    meanConsArray = np.mean(consArray, axis=0)
    stdConsArray = np.std(consArray, axis=0)

    plt.figure(figsize=(8, 5), dpi=160)
    
    plt.subplot(2, 1, 1) #RAW DATA l
    for i in range(len(consArray)):
        plt.plot(consBases, consArray[i], 'b-', linewidth = 1, alpha=0.1)
    plt.grid()
    plt.title('All genomic instances of ' + element_NAME + ' mapped to its consensus \n against ' + annotation_NAME + ' annotation')
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
    saveFigName = os.path.join(outputDir, "MappedToConsensus" + element_NAME + "_" + str(datetime.date.today())+ ".png")
    plt.savefig(saveFigName)
    # plt.show()

    

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
    
    argsIN = parser.parse_args()
    element_NAME= argsIN.element
    annotationFile_NAME= argsIN.annotationFilePath
    # saveOutputFLAG= argsIN.boolean_switch
    outputPath= argsIN.outPath


    #Set Up Output Directory 
    dir = os.path.dirname(outputPath)
    if not os.path.exists(dir): 
            print("\n NOTE: output directory does not exist; will default to TE_annotate_output in current directory")
            os.makedirs("TE_annotate_output")
            outputPath= os.path.join(os.getcwd(), 'TE_annotate_output') 


    [consArray, nameList, consensusLength] = annotateConsensus(element_NAME, annotationFile_NAME)
    plot_MapToConsensus(consensusLength, consArray, element_NAME, 'phastCon', outputPath)
    

#-------
# Main function run 
#-------


if __name__ == "__main__":
    main(sys.argv[1:])
    #create output folder and put contents there 

