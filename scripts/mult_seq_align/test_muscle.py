#!/bin/python
# This script will quickly test of multiple aligments using MUSCLE via biopython. 
#
#
#
# Abin Abraham
# created on: 2017-12-19 15:28:52



# create MER20 fasta input file 
# run MUSCLE via biopython 
# store output into numpy arrays
# plot to assess alignment. 
# plot to overlay annotations 

import pandas as pd
import numpy as np 
import re
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

#-------
# function  
#-------

def mapToConsensus(seq):
    #seq should be a string    
    consIndices = np.array([m.start() for m in re.finditer('[TGAC]', seq)])
    ConsMap = dict(zip(np.arange(len(consIndices)), consIndices))

    return ConsMap

def addMapping(multSeq):
    for i,e in enumerate(multSeq):
        ConsMapDict = mapToConsensus(str(e.seq))


#-------
# main
#-------

## load, read data 
align_FILE = "/Volumes/capra_lab/abraha1/projects/transposable_elements/data/dfam/temp.muscle.out"
# align_FILE  ="/Users/abin-personal/Desktop/transfer/data/temp.muscle.out"
# align_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/temp.muscle.out"

multAlign = AlignIO.read(align_FILE, 'fasta')



oneAlign = multAlign[0]
StrSeq = str(oneAlign.seq)





# create a dict with alignedIndex:elementIndex


# # if False:

# def indexMapping(elem_seq, cons_seq, aligned_elemSeq, aligned_consSeq):
    
#     elem_keepIndex = [m.start() for m in re.finditer('[TGAC]', aligned_elemSeq)]
#     cons_keepIndex = [m2.start() for m2 in re.finditer('[TGAC]', aligned_consSeq)]
#     elem_mapDict = dict(zip( elem_keepIndex, list(range(len(elem_seq)))))
#     cons_mapDict = dict(zip( cons_keepIndex, list(range(len(cons_seq)))))
#     return [elem_mapDict, cons_mapDict]   






# TEbed_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"
# consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
# hg38Fasta_file = "hg38_dfam.nrph.hits_tidy.fasta" 

# # Load Data 
# TEdata = pd.read_table(TEbed_file, sep="\t", header=None)
# TEdata.columns = ["chr", "start", "end", "TE", "TEfamily", "strand", "start_consensus", "end_consesus", "left_consensus"]

# allCons_fasta = SeqIO.index(consensusFasta_file, "fasta", key_function = consensus_renameKey)
# allTE_fasta = SeqIO.index(hg38Fasta_file, "fasta", key_function = TEfasta_renameKey)

# # # set up MUSCLE for MSA 
# # muscle_cline = MuscleCommandline()d
# # #stdout, stderr = muscle_cline()
# # child = subprocess.Popen(str(muscle_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

# # SeqIO.write(storeSeq, child.stdin, "fasta")
# # child.stdin.close()

# # align = AlignIO.read(child.stdout, "fasta")
# # print(align)

# # align_array = np.array( [r.seq for r in align])
