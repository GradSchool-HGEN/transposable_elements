#!/bin/python
# This script provides a wrapper to run TE_annotate_hmm.py
#
# Abin Abraham
# created on: 2017-12-28 09:40:49

import TE_annotate_hmm as ta
import numpy as np
import datetime
import os
import pickle
import argparse
`
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
rootpath = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/"
outputDir = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/TE_annotate_hmm_output/"
     
# USER MUST MODIFFY THESE VARIABLE!!!!
element = "THE1B"
annotation = "phyloP100" #must be the same as the suffix of the intersect file (see intersectFile_PATH)
mappingDict_relativePATH = "align-output_THE1B_01_09_2018/all{}_mappedDict.pi".format(element)

# These input files will automatically update with teh appropriate names. 
mappingDict_PATH = os.path.join(rootpath, mappingDict_relativePATH)
intersectFile_PATH = os.path.join(rootpath,"intersectBed-output", "intersect-output_{}_{}".format(element,annotation))
consensusFasta_PATH = os.path.join(rootpath, "align-input_{}/consensus_{}.fa".format(element,element))

# run annotateConsensus 
annotDict = ta.annotateConsensus(element, annotation, intersectFile_PATH, mappingDict_PATH, consensusFasta_PATH)

saveConsArray_name = "{}_{}_mappedData_mult_{}".format(element,annotation,str(datetime.date.today()))
np.save(os.path.joing(outputDir,saveConsArray_name), annotDict)

ta.plot_MapToConsensus(annotDict, element, annotation)

