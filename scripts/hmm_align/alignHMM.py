#!/bin/python
# This script will create a mapping of a given TE to its consensus and output it in dictionary format.
#   > inputs: 
#       element      : user specified element
#       align-input/ : dir for input files with the following files: 
#           - "nhmmer-output_<element>"
#           - "pHMM_<element>"
#           - "<element>_torun.fa" **user must provide this file**
#           - "consensus_"<element>".fa"
# 
#   > depends: /dors/capra_lab/data/transposable_elements/dfam/Dfam.hmm
#   > depends: /dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta
#   > depends on the following python scripts: 
#           - getFasta_dfamTE
#           - parse_nhmmerOutput
#           - createMapDict
#   
#   > outputs: 
#           - "align-input_<element>/"         : required input files will be created in this dir
#           - "align-output_<element>_<date>/" : outputs will be written here
#       
# 
# Abin Abraham
# created on: 2017-12-26 00:12:31
# modified  : 2018-02-24 07:19:19


#pass element to analyze via command line call 
#create fasta for input TE sequence 
#create pHMM if it doesn't exist 
#run nhmmer 
#parse output 
#map indices to consensus 
#pickle mapping dictionary 
#save QC plot  

import os 
import subprocess
import time
import pickle 
import sys
import argparse
import warnings
import numpy as np 
import matplotlib.pyplot as plt
from Bio import SeqIO 

# these are scripts that Abin made
import getFasta_dfamTE as getFasta
import parse_nhmmerOutput as pnh
import createMapDict as crmd

# -----------
# dependencies 
# ----------- 

### root and dfam 
ROOT_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/LTRs"
DFAM_HMM_FILE = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.hmm"
TE_FASTA_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"

# ROOT_PATH = "/Users/abin-personal/Desktop/desk_transfer/data/LTRs"
# DFAM_HMM_FILE = "/Users/abin-personal/Desktop/desk_transfer/data/Dfam.hmm"
# TE_FASTA_FILE = "/Users/abin-personal/Desktop/desk_transfer/data/hg38_dfam.nrph.hits_tidy.fasta"

# -----------
# functions
# ----------- 

def checkIODir(dirpath):
    if os.path.isdir(dirpath) == False: 
        os.mkdir(dirpath)
        print("created... {}".format(dirpath))

def checkFileExist(filepath, cmd):
    if os.path.isfile(filepath) == False:
        subprocess.run(cmd, shell=True)
        print("\ninput required not found: {}".format(""))
        print("ran cmd: {}".format(cmd)) 

def checkOutputFileExist(filepath):

    if os.path.isfile(filepath) and (os.path.getsize(filepath) > 0): 
        print("\ngood, file exists! {}".format(filepath))
    else: 
        print("\nfile missing or empty: {}".format(filepath))
        warnings.warn('missing file!')
        print("exiting...alignHMM.py")
        os.sys.exit(-1) 
     
def plot_mappingCoverage(element, allSeq, consensusFile, allElemMappedDict,*savePath):
    #allSeq = parsed_nhmmer keys

    confa = SeqIO.read(consensusFile,'fasta')
    consLength = len(confa.seq)
    sumCoverage = np.zeros(consLength)
    sumDeletion = np.zeros(consLength)
    sumInsertions = np.zeros(consLength)

    for oneseq in allSeq: 
        one_elemTOcons = allElemMappedDict[oneseq]
        # bases of shared homology
        bases = crmd.getHomologBases_con(one_elemTOcons)
        bases = [int(x) for x in bases]
        delbases =  crmd.getDeletedBases_con(one_elemTOcons)
        delbases = [int(x) for x in delbases]
        ins_bases = crmd.getInsertBases_con(one_elemTOcons)
        ins_bases = [int(x) for x in ins_bases]
        sumCoverage[bases] += 1
        sumDeletion[delbases]  += 1 
        sumInsertions[ins_bases] += 1 
    print("I am about to call plt")
    # plot alignment coverage over consensus bases 
    fig, ax = plt.subplots()
    ax.plot(np.arange(consLength), sumCoverage,label ='Coverage')
    ax.plot(np.arange(consLength), sumDeletion,label ='Deletion')
    ax.plot(np.arange(consLength), sumInsertions,label ='Insertion')
    ax.set(xlabel='consensus bases for '+element, ylabel='Number of Genomics Hits Per Base',
        title='{} Coverage of Mapping to Consensus Using hmmalign (n={})'.format(element, len(allSeq)))
    ax.grid()
    plt.legend()
    print("I made it to almost saving..")
    if savePath: 
        plt.savefig(savePath[0])
    else:
        # plt.show()
        print("Iam in else")

def main(element):
    start = time.time()
    print("------------------------------------------------------------------------------------------")
    print("creating mapping to consensus for: {}".format(element))


    # input and output directory names 
    OUTPUT_DIR = "align-output_{}__{}/".format(element,time.strftime("%m_%d_%Y"))
    INPUT_DIR = "align-input_{}/".format(element)
    
    #check if dictionary file already exists 
    if os.path.isfile(os.path.join(ROOT_PATH,OUTPUT_DIR,"all{}_mappedDict.pi".format(element))) == True: 
        print("\n{} mapping dictionary already exists, therefore exiting program!".format(element))
        os.sys.exit(-1) 
        
    # input files names 
    NHMMER_OUTPUT_FILE = "nhmmer-output_{}".format(element)
    pHMM_FILE = "pHMM_{}".format(element)
    INPUT_ELEMENT_FASTA_FILE = "{}_torun.fa".format(element)
    CONSENSUS_FASTA_FILE = "consensus_{}.fa".format(element)

    # create input and output directory 
    checkIODir(os.path.join(ROOT_PATH, INPUT_DIR))

    # full paths for input files 
    elemFasta_filepath = os.path.join(ROOT_PATH,INPUT_DIR, INPUT_ELEMENT_FASTA_FILE)
    pHMM_filepath = os.path.join(ROOT_PATH,INPUT_DIR,pHMM_FILE)
    nhmm_filepath = os.path.join(ROOT_PATH,INPUT_DIR,NHMMER_OUTPUT_FILE)
    consSeq_filepath = os.path.join(ROOT_PATH,INPUT_DIR,CONSENSUS_FASTA_FILE)

    # create input TE fasta file if it does not exists 
    if os.path.isfile(elemFasta_filepath) == False:
        getFasta.extract_sequences(element,TE_FASTA_FILE, elemFasta_filepath)

    checkOutputFileExist(elemFasta_filepath)

    # create requried input files and check if they were made properly  
    checkFileExist(pHMM_filepath, "hmmfetch -o {} {} {}".format(pHMM_filepath, DFAM_HMM_FILE, element))
    checkOutputFileExist(pHMM_filepath)
    
    checkFileExist(nhmm_filepath, "nhmmer -o {} --notextw {} {}".format(nhmm_filepath, pHMM_filepath, elemFasta_filepath))
    checkOutputFileExist(nhmm_filepath)
    
    checkFileExist(consSeq_filepath, "hmmemit -c -o {} {}".format(consSeq_filepath, pHMM_filepath))
    checkOutputFileExist(consSeq_filepath)

    # parse the nhmmer output
    parsed_nhmmer = pnh.parseTermOut(nhmm_filepath)
    allSeq = list(parsed_nhmmer.keys())
    
    # for each sequence, create a dictionary with {sequence:mapping}
    allElemMapped_Dict = dict()
    for oneseq in allSeq: 
        elemMap_dict = crmd.createMapDict(oneseq, parsed_nhmmer)
        allElemMapped_Dict[oneseq] = elemMap_dict


    # create output directory and pickle dictionary
    checkIODir(os.path.join(ROOT_PATH,OUTPUT_DIR))
    pickle.dump(allElemMapped_Dict, open(os.path.join(ROOT_PATH, OUTPUT_DIR, "all{}_mappedDict.pi".format(element)), 'wb'))

    # QC plots 
    plot_mappingCoverage(element, allSeq, os.path.join(ROOT_PATH,INPUT_DIR,CONSENSUS_FASTA_FILE), allElemMapped_Dict,os.path.join(ROOT_PATH,OUTPUT_DIR,"QC_{}_hmmalign_coverage_{}.eps".format(element,time.strftime("%m_%d_%Y"))))

    stop = time.time() 
    print("---------------------------------------------")
    print("took {} minutes to finish\n".format(round((stop-start)/60,3)))
    print("check output in folder: {}/{}".format(ROOT_PATH,OUTPUT_DIR))
    print("---------------------------------------------")


# -----------
# main
# ----------- 

if __name__ == "__main__":
    # set up argparse 
    parser = argparse.ArgumentParser(description="will output a dictionary with mapping of TE to its consensus")
    parser.add_argument("element", default='None',
                    action='store', type=str,
                    help="Select a transposable element") 
    
    argsIN = parser.parse_args()
    element= argsIN.element
    # element='MER20'

    main(element)
    #create output folder and put contents there


