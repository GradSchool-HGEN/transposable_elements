#!/bin/python
# This script will create a dictionary that maps a transposable element to its consensus based on alignment from nhmmer.
#
#
#
# Abin Abraham
# created on: 2017-12-25 15:31:22

import collections
import numpy as np

# import parse_nhmmerOutput as par
# parsedOutput = par.parseTermOut("nhmmer-output-terminal_allMER20")
# thisSeq = 'MER20::chr12:33159061-33159182(-)'

def createMapDict(thisSeq, parsedOutput):
    # load data for thisSeq
    thisElemSeq = parsedOutput[thisSeq]["elemSeq"]
    thisModelSeq = parsedOutput[thisSeq]['model_seq']
    thisElem_start = int(parsedOutput[thisSeq]["elemSeqStart"])
    thisElem_end = int(parsedOutput[thisSeq]["elemSeqEnd"])
    thisModel_start = int(parsedOutput[thisSeq]["model_start"])
    thisModel_end = int(parsedOutput[thisSeq]["model_end"])
    matchedBases = parsedOutput[thisSeq]["matchedBases"]

    elemMask = np.empty(len(thisElemSeq))
    modelMask = np.empty(len(thisModelSeq))

    ecounter = thisElem_start - 1 #change from 1 based to 0 based
    ecounterNeg = -1
    for i, x in enumerate(thisElemSeq):  # 
        if x.isupper() or x.islower():             
            elemMask[i] = ecounter
            ecounter += 1
        else:
            elemMask[i] = ecounterNeg      #negative is deletion
            ecounterNeg -= 1

    mcounter = thisModel_start - 1 #change from 1 based to 0 based
    mcounterNeg = -1
    for i, x in enumerate(thisModelSeq):
        if x.isupper() or x.islower():
            modelMask[i] = mcounter
            mcounter += 1
        else:
            modelMask[i] = mcounterNeg #negative value => insertion
            mcounterNeg -= 1

    elemTOcons = collections.OrderedDict((zip(elemMask, modelMask)))

    return elemTOcons

def getHomologBases_con(mapDict):
    saveKey = set()
    for k in mapDict: 
        if mapDict[k] > -1 and int(k) > -1: 
            saveKey.add(mapDict[k])

    return saveKey

def getInsertBases_con(mapDict):
    """ return the previous index before the first start of an insertion """
    storeval = list(mapDict.values())
    saveInsert = set() 
    runningCounter = 0; 
    for i,k in enumerate(mapDict):
        if mapDict[k] < 0:
            if i == 0: 
                saveInsert.add(0)
                runningCounter +=1 
            else:
                if runningCounter == 0: 
                    saveInsert.add(storeval[i-1]) #save the index before insertions
                    runningCounter +=1 
        else:
            runningCounter = 0; 

    return saveInsert

def getDeletedBases_con(mapDict): 
    saveDeleted = set() 
    for k in mapDict: 
        if mapDict[k] > -1 and int(k) < 0: 
            saveDeleted.add(mapDict[k])

    return saveDeleted


# for elemSeq
    # -2: deletion in elemSeq compared to consensus
    # -1: insertion in elemSeq compared to consensus
    # 0,1,2... element index for bases that match to consensus
            # update counter for each base in element

# for model Seq
    #  0,1,2... consensus index for ALL bases, whether there is match or not
    #  -1: deleted bases in consensus
    #  -2: insertion

# output will be a dictionary with element to consensus mapping
    # keys and value > 0 will represent the index of the element sequence and its mapping ot the index in the consensus seque
    # key >0, value ==-2 represent insertions in the element
    # key =-2, value > 0 represent deletions in the  element

