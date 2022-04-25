#!/bin/python
# This script will calculate the Levenstein's distance (LVD) between 
# transposable elements (obtained from RepeatMasker) and their consesnsus sequence (from Dfam). 
# Writes a dictionary (element:LVDvalues) using json to current dir. 
#
# Dependencies:
#          1) TE BED FILE: "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"   
#          2) TE consensus sequence: "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
#          3) hg38 fasta: "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 
# Abin Abraham
# created on: 2017-12-04 08:31:31
# modified on: 2017-12-05 10:04:25


import pandas as pd 
import editdistance as ed 
import json

from Bio import SeqIO
from Bio import pairwise2
from Bio import AlignIO


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

def uniqTE(TEds):
    uniqTE_list = TEds.TE.unique() 
    return uniqTE_list

def getTEcoordinates(element,TEds):
    elemLoc_list =  TEds.loc[TEds.TE == element, ['chr','start','end']].values.tolist()
    elemLoc_list = [[x[0]+ ":" + str(x[1]) + "-" + str(x[2])] for x in elemLoc_list]
    allElem_seq = [ allTE_fasta[x[0]].seq.upper() for x in elemLoc_list]
    return allElem_seq


# -----------
# LOAD DATA 
# ----------- 

# File Paths 
TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"
consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
hg38Fasta_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 

# TE genomic coordinates  
TEdata = pd.read_table(TEbed_file, sep="\t", header=None)
TEdata.columns = ["chr", "start", "end", "TE", "TEfamily", "strand", "start_consensus", "end_consesus", "left_consensus"]

# Consensus Sequences
allCons_fasta = SeqIO.index(consensusFasta_file, "fasta", key_function = consensus_renameKey)
# TE sequences 
allTE_fasta = SeqIO.index(hg38Fasta_file, "fasta", key_function = TEfasta_renameKey)

#-------
# MAIN 
#-------

elemList = uniqTE(TEdata)

store_LVD = {} 
for elem in elemList: 
     if elem in allCons_fasta.keys():
        allElem_seq = getTEcoordinates(elem,TEdata)
        con_seq = allCons_fasta[elem].seq.upper()
        edValues = [ ed.eval(str(x), str(con_seq)) for x in allElem_seq ]
        store_LVD[elem] = edValues
        print(elem)
       
json.dump(store_LVD, open("LVDvalues.txt",'w'))

# to load LVDvalues.txt...
# reload_store_LVD = json.load(open("LVDvalues.txt"))


# PLOT 

# violin plots  of LVD to consensus by TE family

# violin plots  of LVD to consensus by TE element per family 

# %overlap w/ CRE by LVD for each TE family 
