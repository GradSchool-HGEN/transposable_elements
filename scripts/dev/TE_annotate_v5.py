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
#       > reformat dfam file input to be 0 based exclusive (bedtools requires bed format)
#       > figure out a way to sort files before passing it onto bedtools intersect (or check first then sort)
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
def checkInputs(annotationfile_PATHc):
    numOfColumns = subprocess.check_output(['awk', '{print NF; exit}', annotationfile_PATHc], universal_newlines=True).splitlines()

    cmd = str("sort -k1,1 -k2,2n " + annotationfile_PATHc + " -o " + annotationfile_PATHc).split()
    subprocess.call(cmd)

    fileModifiedFLAG = False
    if int(numOfColumns[0]) > 4:
        fileModifiedFLAG = True

        cmd1 = 'cut -f1,2,3,4 '+annotationfile_PATHc
        cmd1 = cmd1.split()
        f = open("formated_annotationfile.tsv", "w")
        subprocess.call(cmd1, stdout=f)
        f.close()

    return fileModifiedFLAG

def annotateConsensus(element, annotationfile_PATH):
    TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
    filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)

    fileFLAG = checkInputs(annotationfile_PATH)
    if fileFLAG:
        intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', "formated_annotationfile.tsv", '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)  
    else:
        intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotationfile_PATH, '-b', 'stdin', '-sorted'], stdin = filterForElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)
   
    #=========== load data into data frame ===========
    out_data =[x.split('\t') for x in intersectOutput.split('\n')]
    col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "TE_genomic_chr","TE_genomic_start","TE_genomic_end", "TE_model","model_start","model_end", "model_length" ] 
    df = pd.DataFrame(data=out_data[:-1], columns=col_name)

    #=========== map TE to consensus ===========
    df["TE_genomicHit_name"] = df["TE_genomic_chr"] + ":" + df["TE_genomic_start"].map(str) + "-" + df["TE_genomic_end"].map(str)
    df = df.apply(pd.to_numeric, errors='ignore')
    df["mapToConsensus_start"] = df.intersect_start - df.TE_genomic_start + df.model_start
    df["mapToConsensus_end"] = df.intersect_end - df.TE_genomic_start + df.model_start
    df = df.sort_values(by='TE_genomicHit_name')
    gb = df.groupby('TE_genomicHit_name')

    # Organize Data for Output
    consensusLength = df.model_length.iloc[1]
    consArray = np.zeros((1, consensusLength))
    tmp_consArray = np.zeros((1, consensusLength))
    mask_consArray = np.ones((1, consensusLength))

    nameList = [list(gb)[i][0] for i in range(len(gb))]

    for i in range(len(gb)):
        oneElem_mapped = list(gb)[i][1][["mapToConsensus_start", "mapToConsensus_end","annotation_score"]].values.astype(float)

        for s, e, v in oneElem_mapped:
            #check if these indices/bases already had a value, if so take the average 
            tmp_consArray[0, int(s):int(e)] = (tmp_consArray[0, int(s):int(e)] + v)/mask_consArray[0, int(s):int(e)]
            #check if these indices/bases already had a value, if so take the sum 
            # tmp_consArray[0, int(s):int(e)] = (tmp_consArray[0, int(s):int(e)] + v)
            
            mask_consArray[0, int(s):int(e)] += 1
        
        #append sequence to conArray, THEN reset for next TE sequence
        consArray = np.vstack((consArray, tmp_consArray))
        tmp_consArray = np.zeros((1, consensusLength))
        mask_consArray = np.ones((1, consensusLength))

    consArray = consArray[1:, :]

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
    plt.savefig(saveFigName)
    plt.subplots_adjust(hspace=0.5)
    # plt.show()
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

