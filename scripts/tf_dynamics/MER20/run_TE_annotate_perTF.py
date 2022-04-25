#!/bin/python
# This script provides a wrapper to run TE_annotate_hmm.py. 
# Will run TE_annotate_hmm.py per TF for all TFs within a set of TEs 
# 
# 
# Abin Abraham
# created on: 2018-02-26 17:04:50
# 
import sys
sys.path.append("/dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align")



import TE_annotate_hmm as ta
import numpy as np
import datetime
import os
import pickle
import argparse
import glob 
import re



#-------
# functions
#-------

def saveDict(saveName, dictToSave):
     with open(saveName, 'wb') as fs: 
         pickle.dump(dictToSave, fs)

def checkFileExist(filepath, cmd):
    if os.path.isfile(filepath) == False:
        subprocess.run(cmd, shell=True)
        print("\ninput required not found: {}".format(""))
        print("ran cmd: {}".format(cmd)) 

# -----------
# section title
# ----------- 

# USER MUST MODIFY THESE VARIABLES 
# rootpath refers to where 1) mappingDict_PATH 2)intersectFile_PATH and 3) consensusFasta_PATH will be built on 
rootpath = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/"
outputDir = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/"
     
# USER MUST MODIFFY THESE VARIABLE!!!!
element = "MER20"
# annotation = "this_TF" #must be the same as the suffix of the intersect file (see intersectFile_PATH)
mappingDict_relativePATH = "all{}_mappedDict.pi".format(element)

# These input files will automatically update with teh appropriate names. 
mappingDict_PATH = os.path.join(rootpath, mappingDict_relativePATH)
consensusFasta_PATH = os.path.join(rootpath, "consensus_{}.fa".format(element))
# intersectFile_PATH = os.path.join(rootpath, "intersect/intersect-output_{}_{}".format(element,annotation))

intersect_file_list =glob.glob("/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/intersect/inter-TE-anno_by*.tsv")

# run annotateConsensus 
for this_file in intersect_file_list: 
    this_tf = re.search("[A-Z0-9]*.tsv",this_file).group()[:-4]
    annotDict = ta.annotateConsensus(element, this_tf, this_file, mappingDict_PATH, consensusFasta_PATH)
    saveConsArray_name = "{}_{}_mappedData_mult_{}".format(element,this_tf,str(datetime.date.today()))
    np.save(os.path.join(outputDir,saveConsArray_name), annotDict)

# ta.plot_MapToConsensus(annotDict, element, annotation)

