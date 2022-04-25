#!/bin/python
# This script will overlay a given annotation for a transposable element or a list of elements over its consensus sequence
#
# Abin Abraham      
# created on: 2017-12-07 20:15:16
# 
# Depends on: 
#       --> TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
# 
# Inputs: 
#       1) transposable element of interest 
#       2) annotation file e.g: annotationFile_NAME = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"
# 
# Output: 
#       1) numpy array where each row is one genomic instance of the transposable element 
#       2) list of the names of  genomic sequences that correspond to each row in the above row

# run TE_annotate.py "SVA_E" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

import pandas as pd
import numpy as np
import subprocess
import argparse
import datetime

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

argsIN = parser.parse_args()


element_NAME = argsIN.element
annotationFile_NAME = argsIN.annotationFilePath

# -----------
# MAIN
# ----------- 
# File Paths 
TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"

# Filter for specific TE element & run intersect; e.g: element_NAME = "SVA_E"

filteredFile = subprocess.Popen(['grep', element_NAME, TEbed_file ], stdout=subprocess.PIPE, universal_newlines=True)
obs_intersect = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', annotationFile_NAME, '-b', 'stdin'], stdin = filteredFile.stdout,  universal_newlines=True)

out_data =[x.split('\t') for x in obs_intersect.split('\n')]
col_name = ["intersect_chr", "intersect_start", "intersect_end", "annotation_score", "genomic_chr","genomic_start","genomic_end", 
            "TE_model","model_start","model_end", "model_length" ] 

df = pd.DataFrame( data = out_data[:-1], columns=col_name)
df["genomicHit_name"] = df["genomic_chr"] + ":" + df["genomic_start"].map(str) + "-" + df["genomic_end"].map(str)
df = df.apply(pd.to_numeric, errors = 'ignore')
df["mapToConsensus_start"] = df.intersect_start - df.genomic_start + df.model_start
df["mapToConsensus_end"] = df.intersect_end - 1 - df.genomic_start + df.model_start


# Organize Data for Output 
consensusLength = df.model_length.iloc[1]
df = df.sort_values(by='genomicHit_name')
tmp_consArray = np.zeros((1, consensusLength))
consArray = tmp_consArray
gb = df.groupby('genomicHit_name')
nameList = [list(gb)[i][0] for i in range(len(gb))]

for i in range(len(gb)):
    spec = list(gb)[i][1][["mapToConsensus_start", "mapToConsensus_end","annotation_score"]].values.astype(float)
    # check for overlap 
    st = np.sort(spec[:,0].astype(int))
    ed = np.sort(spec[:,1].astype(int))
    if (st[:-1][st[1:] == st[:-1]].size != 0) | (ed[:-1][ed[1:] == ed[:-1]].size != 0): 
        print("there are duplicate annotations for a given genomic instance of the TE")

    for s,e,v in spec: tmp_consArray[0, int(s):int(e)+1]=v
    consArray = np.vstack((consArray, tmp_consArray))

consArray = consArray[1:, :]


if False:

    #-------
    # PLOT RAW DATA 
    #-------
    import matplotlib.pyplot as plt

    consBases = np.arange(0.0, consensusLength, 1, dtype=int)
    meanConsArray = np.mean(consArray, axis=0)
    stdConsArray = np.std(consArray, axis=0)

    plt.subplot(2, 1, 1) #RAW DATA 
    for i in range(len(consArray)):plt.plot(consBases, consArray[i], 'b-', linewidth = 1, alpha=0.1)
    plt.grid()
    plt.title('Raw Data (Top) and Mean +/- 2SD (bottom) for \n Mapping to Consensus for '+ element_NAME)
    plt.ylabel('Annotation')
    [x1,x2,y1,y2] = plt.axis()

    plt.subplot(2, 1, 2) #MEAN and STD 
    plt.plot(consBases, meanConsArray, 'r-', linewidth = 1.2, alpha=1)
    plt.plot(consBases, meanConsArray + 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.plot(consBases, meanConsArray - 2*stdConsArray, 'b:', linewidth = 0.7, alpha=0.8)
    plt.axis((x1,x2,-1*y2,y2))
    plt.grid()
    plt.xlabel('consensus bases (1start)')
    plt.ylabel('Annotation')
    plt.savefig("RawMapping"+element_NAME+"_"+ datetime.date.today()+".png")
    plt.show()



    # TO DO 
    # demand specification of at least one TE 
    # what to do if there is overlap when filling consensus numpy array 
    # assumes that dfam coordinates are 1 start inclusive ends 
    # assumes that annotation data is 0 start exclusive ends 

    # Input file requirements 
        # annotation file: first three columns must be chr, start, end (in hg38 genomic coordinates)
        # annotation file can then have 

    # parser.add_argument("-a", 
    #                 action='store', type=str,
    #                 dest='anno_file_path', 
    #                 help="Path or file name for annotation")    