#!/bin/python
# """For a given TE and an annotation file, this script will map the intersecting annotation to the consensus sequence of TE. """"
#
# Abin Abraham
# created on: 2017-12-20
#
# Depends on:
#        > TEbed_file: "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# 
# Required Inputs:
#        > element: name of transposable element
#        > annotation file path: this file should be tab seperated with following columns: chr, start, end, annotation value
#        > multilpleAlignment file: "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/testMSA_muscle_maxiter2.out"
#
# Output:
#        > consArray: numpy array with annotation values; each row = one locus of TE
#        > nameList: list of TE genomic locations (chr:start-end) in order of consArray rows
# Note:
#       ** all input files are assumed be 0-start exclusive

import pandas as pd
import numpy as np
import subprocess
import getMappingDict as getmap



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

def annotateConsensus(element, annotationfile_PATH, multSeqAlig_PATH):
    #***** dependencies *******
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
    df.drop( ["model_start", "model_end", "model_length"], axis=1, inplace=True)

    #=========== map TE to consensus ===========
    df["TE_genomicHit_name"] = df["TE_genomic_chr"] + ":" + df["TE_genomic_start"].map(str) + "-" + df["TE_genomic_end"].map(str)
    df = df.apply(pd.to_numeric, errors='ignore')
    df['offset_start'] = df['intersect_start'] - df['TE_genomic_start']
    df['offset_end'] = df['intersect_end'] - df['TE_genomic_start']
    gb = df.groupby('TE_genomicHit_name')
    nameList = [list(gb)[i][0] for i in range(len(gb))]

    mappingDictList = getmap.createMapDict(multSeqAlig_PATH)

    # Organize Data for Output
    consensusLength = list(mappingDictList.values())[0][1]
    consArray = np.zeros((1, consensusLength))
    tmp_consArray = np.zeros((1, consensusLength))
    mask_consArray = np.ones((1, consensusLength))


    for i in range(len(gb)): #loop over every loci for selected TE
        oneSeq = element+"::"+list(gb)[i][0]
        oneElem = list(gb)[i][1][["offset_start", "offset_end", "annotation_score"]].values.astype(float)
        oneDict = mappingDictList[oneSeq][0]

        for s, e, v in oneElem: #loop over every annotation at this specific locus
            cons_indices = oneDict[int(s):int(e)]
            #add to previous value to keep a total sum
            tmp_consArray[0, cons_indices] = (tmp_consArray[0, cons_indices] + v)
            mask_consArray[0, cons_indices] += 1
        
        # print(i)
        tmp_consArray = tmp_consArray/mask_consArray
        #tmp_consArray = np.divide(tmp_consArray/mask_consArray)

        consArray = np.vstack((consArray, tmp_consArray)) #append sequence to conArray
        tmp_consArray = np.zeros((1, consensusLength)) #reset for next locus of element


    consArray = consArray[1:, :]

    return [consArray, nameList]



def main(element, annotationFilePATH, multipleAlignment_PATH):

    [consensusArray, nameList] = annotateConsensus(element, annotationFilePATH, multipleAlignment_PATH)

    return consensusArray
    #multSeqAlig_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/temp.muscle.out"

'''
  #Set Up Output Directory
    dir = os.path.dirname(outputPATH)
    if not os.path.exists(dir):
            print("\n NOTE: output directory does not exist; will default to TE_annotate_output in current directory")
            os.makedirs("TE_annotate_output")
            fullOutputPath= os.path.join(os.getcwd(), 'TE_annotate_output')
'''