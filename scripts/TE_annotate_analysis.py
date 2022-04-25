#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2017-12-13 16:40:43


# analyze phastCon score on TE; compare to known TF binding motifs 
# first example, lets looks at MER20 

# map all instances of MER20 phastCon score onto consensus 
# in another plot, map all instances of MER20 onto consensus but this time look at TF motifs 


# for each TE, calculate the correlation between phastCon Score and TF motif at a sequence levels 
# TE_annotate.py "SVA_E" "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"

import TE_annotate_hmm as ta
import TE_annotate_consensusTF as tc 
import numpy as np
import datetime
import os
import pickle


#-------
# functions
#-------

def saveDict(saveName, dictToSave):
    with open(saveName, 'wb') as fs: 
        pickle.dump(dictToSave, fs)
    
# -----------
# main
# ----------- 
outputDir = "/dors/capra_lab/abraha1/projects/transposable_elements/results"
element = "MER20"


#-------
# create mapping to consensus for phastCon Score 
#-------
if False: 
    annotationfile_PATH = '/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test'
    multSeqAlig_PATH = '/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/multAlign/mafft-default_MER20'
    mapConDict, seqList = ta.annotateConsensus(element, annotationfile_PATH)

# -----------
# Analysis not using Mult Alignments
# ----------- 

if True: 
    annotationfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/phastCon_test"
    #annotationfile_PATH = "/dors/capra_lab/abraha1/data/hg38_phastCon100way/hg38.phyloP100way/hg38.phyloP100way.bed"
    [consArray_phastCon, nameList, consensusLength] = ta.annotateConsensus(element, annotationfile_PATH)
    saveName = element+"_phyloP100_"+str(datetime.date.today())+".p"
    saveName = os.path.join(outputDir,saveName)
    # np.save(saveName, consArray_phastCon)
    
    ta.plot_MapToConsensus(consensusLength, consArray_phastCon, element, 'phastCon', outputDir)
    
    #pickle save 
    # dict_toSave = dict(zip(nameList, consArray_phastCon))
    # saveDict(saveName, dict_toSave)
    # print(saveName)
    # ta.plot_MapToConsensus(consensusLength, consArray_phastCon, element, 'phyloP100_FULL', outputDir)

if False: 
    #-------
    # create mapping for overlap with enhancers
    #-------
    annotate_enhancer_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/all_facet-organ_enhancers_coordinates-formated.tsv"
    [consArray_enhancer, nameList_enhancer, consensusLength_enhancer] = ta.annotateConsensus(element, annotate_enhancer_PATH)
    saveName2 = element+"_allFacetEnhancer_"+str(datetime.date.today())+".p"
    saveName2 = os.path.join(outputDir,saveName2)
    # np.save(saveName2, consArray_enhancer)
if False: 

    #pickle save 
    dict_toSave2 = dict(zip(nameList_enhancer, consArray_enhancer))
    saveDict(saveName2, dict_toSave2)
    print(saveName2)
    # ta.plot_MapToConsensus(consensusLength_enhancer, consArray_enhancer, element, 'enhancerOverlap', outputDir)
if False: 
    #-------
    # create mapping for overlap with TFmotifs on facet_enhancers 
    #-------
    annotate_TF_on_enhancer_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/combined_fimo-out_all_facet_enhancers-FILTERED-withTFcoords_formatted.tsv"
    [consArray_TFenhancer, nameList_TFenhancer, consensusLength_TFenhancer] = ta.annotateConsensus(element, annotate_TF_on_enhancer_PATH)
    saveName3 = element+"_TF_on_Enhancer_"+str(datetime.date.today())+".p"
    saveName3 = os.path.join(outputDir,saveName3)
    # np.save(saveName3, consArray_TFenhancer)
    
    #pickle save 
    dict_toSave3 = dict(zip(nameList_TFenhancer, consArray_TFenhancer))
    saveDict(saveName3, dict_toSave3)
    print(saveName3)
    # ta.plot_MapToConsensus(consensusLength_TFenhancer, consArray_TFenhancer, element, 'TF_on_enhancer_Overlap', outputDir)
    
if False: 
    #-------
    # create mapping for TF motifs on consensus 
    #-------# 
    consArray_consensus = tc.mapTF_toConsensus(element)

    saveName4 = element+"_TF_on_consensus_"+str(datetime.date.today())
    saveName4 = os.path.join(outputDir,saveName4)
    # np.save(saveName4, consArray_consensus)
    tc.plot_TFtoConsensus(consArray_consensus, element, outputDir)
