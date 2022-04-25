#!/bin/python
# This script will annotate transposable element at a base pair level
#
#
#   
# Abin Abraham
# created on: Dec 4 2017


import pandas as pd
import numpy as np 
import re
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from io import StringIO
import subprocess, sys
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

def mapToAligned(elem_seq, aligned_seq):
    #check format of incoming variables 
    elem_seq.

    if not elem_seq.isupper()
        print('raise error, element sequences is upper case')
    if not aligned_seq.isupper()
        print('raise error, aligned sequences is upper case')
    if elem_seq.find('-') == -1:
        print("there are "-" in the consensus sequence")
        #change to better error handling 
    
    keepIndex = [m.start() for m in re.finditer('[TGAC]', aligned_seq)]
    mapDict = dict(zip( list(range(len(elem_seq))), keepIndex ))
    return mapDict    



# -----------
# LOAD DATA 
# ----------- 

# File Paths 
TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"
consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
hg38Fasta_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 

# Load Data 
TEdata = pd.read_table(TEbed_file, sep="\t", header=None)
TEdata.columns = ["chr", "start", "end", "TE", "TEfamily", "strand", "start_consensus", "end_consesus", "left_consensus"]

allCons_fasta = SeqIO.index(consensusFasta_file, "fasta", key_function = consensus_renameKey)
allTE_fasta = SeqIO.index(hg38Fasta_file, "fasta", key_function = TEfasta_renameKey)

# -----------
# MAIN
# ----------- 

# get TE coordinates & sequence
element = "Eulor11"
elemLoc = TEdata.loc[TEdata.TE == element, ['chr','start','end']]
elemLoc_list = elemLoc.values.tolist()
elemLoc_list = [[x[0]+ ":" + str(x[1]) + "-" + str(x[2])] for x in elemLoc_list]


[ allTE_fasta[x[0]].seq for x in elemLoc_list]

#get consensus and genomic sequences  
con_seq = allCons_fasta[element].seq
elem_seq = allTE_fasta[elemLoc_list[1][0]].seq # get one element fasta seq
storeSeq = [ allTE_fasta[x[0]] for x in elemLoc_list]

#pairwise alignment
align = pairwise2.align.localmd(con_seq.upper(), elem_seq.upper(), 1, -1, -10, -5, -3, -1)



#map TE coordinate into conensus aligned coordinates 
    #check for any gaps in consensus sequence --> if yes, raise error 
    #map genomic seq coordinates to align coordinates 
        # do this by creating a dict key:val assignment 






#align = pairwise2.align.localms(con_seq.upper(), elem_seq.upper(), 1, -1, -3, -1 )
#[[x[2], x[3], x[4] ] for x in align]
#align = pairwise2.align.localms(con_seq.upper(), elem_seq.upper(), 1, -1, -3, -1 )



# # set up MUSCLE for MSA 
# muscle_cline = MuscleCommandline()
# #stdout, stderr = muscle_cline()
# child = subprocess.Popen(str(muscle_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

# SeqIO.write(storeSeq, child.stdin, "fasta")
# child.stdin.close()

# align = AlignIO.read(child.stdout, "fasta")
# print(align)

# align_array = np.array( [r.seq for r in align])



# which TE are you interested in? 
# check to see if a consesus exists 
# biopython pairwise seq to align TE w/ consensus
    #obtain fasta sequence of TE in the genome 
        # use bedtools getfasta 
    #obtain fasta sequence of CONSENSUS TE sequence 
        # use biopython SEQI/O 
    #biopython pariwise align; obtain coordinate of fragment within consensus
    #do this for all instances of TE of interest in the genome

    #obtain hg38 genomic coordinates of TE - chr, start, end 


#-------
# TO DO 
# default align parameters: 
#   align = pairwise2.align.localxx(con_seq.upper(), elem_seq.upper())
#   consider optimizing this.

# References
# alignment output contains 1) both aligned sequences 2) score 3)start 4) end) 
# global = best align between all bases of two seq 
# loacl = best align between fragments of two seq

# storeSeq = [allTE_fasta[elemLoc_list[1][0]], allTE_fasta[elemLoc_list[1][0]], allTE_fasta[elemLoc_list[2][0]]]]
# SeqIO.write(storeSeq, "example.fa", "fasta")


#storeSeq = [ allTE_fasta[x[0]] for x in elemLoc_list]

