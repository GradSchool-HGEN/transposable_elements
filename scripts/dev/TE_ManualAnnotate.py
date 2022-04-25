#!/bin/python
# This script will align all genomic instances of a transposable element to its consensus sequence. 
# Will output an array with start and end index of alignment occuring to a consensus sequence. 
#
#   
# Abin Abraham
# created on: Dec 4 2017
#
#-------
# This script depends on:   
#    TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"
#    consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
#    hg38Fasta_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 
#-------


import pandas as pd
import re
import numpy as np 
import time 

from multiprocessing import Pool 
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from functools import partial

# -----------
# FUNCTIONS
# ----------- 

def consensus_renameKey(id):
    newName = id.split("-")[0:-1]
    newName = '-'.join(newName)
    return newName

def TEfasta_renameKey(id):
    newName = id.split("(")[0]
    return newName

def getTEcoordinates(element, TEdataset, allTE_sequences):
    elemLoc_list =  TEdataset.loc[TEdataset.TE == element, ['chr','start','end']].values.tolist()
    strand_list = TEdataset.loc[TEdataset.TE == element, ['strand']].values.tolist()
    elemLoc_list = [[x[0]+ ":" + str(x[1]) + "-" + str(x[2])] for x in elemLoc_list]
    allElem_seq = []
    for i,x in enumerate(elemLoc_list): 
        if strand_list[i][0] =="+": 
            #allElem_seq[i] =  allTE_sequences[x[0]].seq.upper()
            allElem_seq.append(allTE_sequences[x[0]].seq.upper())
        else: 
            #allElem_seq[i] =  allTE_sequences[x[0]].seq.complement().upper()
            allElem_seq.append(allTE_sequences[x[0]].seq.complement().upper())
            
    # allElem_seq = [ allTE_sequences[x[0]].seq.upper() for x in elemLoc_list] --> ignoring strand
    
    allElem_seq_loc = {}
    for i in range(len(elemLoc_list)):
         allElem_seq_loc[str(allElem_seq[i])] = elemLoc_list[i][0]
   
    return [allElem_seq_loc, allElem_seq]

def indexMapping(elem_seq, cons_seq, aligned_elemSeq, aligned_consSeq):
    
    elem_keepIndex = [m.start() for m in re.finditer('[TGAC]', aligned_elemSeq)]
    cons_keepIndex = [m2.start() for m2 in re.finditer('[TGAC]', aligned_consSeq)]
    elem_mapDict = dict(zip( elem_keepIndex, list(range(len(elem_seq)))))
    cons_mapDict = dict(zip( cons_keepIndex, list(range(len(cons_seq)))))
    return [elem_mapDict, cons_mapDict]   


def getIndex(oneElem, consensusSeq, indexMappingFun, seq_loc_dict):
    align = pairwise2.align.localmd(consensusSeq, oneElem, 1, -1, -10, -5, -3, -1)
    [eDict, cDict] = indexMappingFun(oneElem, con_seq, align[0][1], align[0][0])
    [s , e] = [align[0][3] , align[0][4]]
    cons_START = cDict[s]
    cons_END = cDict[e-1] # e is the end index for slicing purposes, so have to subtract 1 
    elem_START = eDict[s]
    elem_END = eDict[e-1]
    return [cons_START, cons_END, elem_START, elem_END, seq_loc_dict[str(oneElem)]]   


# -----------
# LOAD DATA 
# ----------- 

# File Paths -
TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"
consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
hg38Fasta_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 

# Load Data 
TEdata = pd.read_table(TEbed_file, sep="\t", header=None)
TEdata.columns = ["chr", "start", "end", "TE", "TEfamily", "strand", "start_consensus", "end_consesus", "left_consensus"]

allCons_fasta = SeqIO.index(consensusFasta_file, "fasta", key_function = consensus_renameKey)
hg38_TE_fasta = SeqIO.index(hg38Fasta_file, "fasta", key_function = TEfasta_renameKey)

# -----------
# MAIN
# ----------- 
element='MER20'
# get all sequences for a specific TE  
#[allElem_seq, elemLoc_list] = getTEcoordinates(element, TEdata, hg38_TE_fasta)

[allElem_seq_loc, allElem_seq]  = getTEcoordinates(element, TEdata, hg38_TE_fasta)

#get consensus  sequences  
con_seq = str(allCons_fasta[element].seq).upper()

# start pooling 
p = Pool() 
partial_getIndex = partial(getIndex, consensusSeq=con_seq, indexMappingFun=indexMapping, seq_loc_dict = allElem_seq_loc)

t0 = time.time()
results = p.map(partial_getIndex, allElem_seq) 
#results =   [cons_START, cons_END, elem_START, elem_END, seq_loc_dict[oneElem]]   
total = time.time() - t0 

# write results to a tsv file 
# genomic coordinate (chr, start, end), element, consensus distance from start, TE distance from start
# 0 based start, last coordinate is inclusive 

chrm = [results[i][4].split(':')[0] for i in range(len(results))]
start = [results[i][4].split(':')[1].split("-")[0] for i in range(len(results))]
end  =  [results[i][4].split(':')[1].split("-")[1] for i in range(len(results))]
elem = [element]*len(results)
cons_start = [results[i][0] for i in range(len(results))]
cons_end = [results[i][1] for i in range(len(results))]
elem_start = [results[i][2] for i in range(len(results))]
elem_end = [results[i][3] for i in range(len(results))]

data = {'chr': chrm,\
        'start': start, \
        'end': end, \
        'elem': elem, \
        'cons_start': cons_start, \
        'cons_end': cons_end, \
        'elem_start': elem_start, \
        'elem_end': elem_end} 

colOrder = [ 'chr', 'start', 'end', 'elem', 'cons_start','cons_end', 'elem_start','elem_end']
frame = pd.DataFrame(data)
frame = frame.reindex(columns=colOrder)

frame.to_csv('MER20.dataset', sep='\t', header=True, index=False )
