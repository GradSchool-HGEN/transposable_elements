#!/bin/python
#This script will take alignments of several elements and stack them onto its consensus sequence.  
#
#
#
# Abin Abraham
# created on: 2017-12-07 10:00:38

import pandas as pd 
import numpy as np 
from Bio import SeqIO
import matplotlib.pyplot as plt 


#-------
# FUNCTIONS
#-------
def consensus_renameKey(id):
    newName = id.split("-")[0:-1]
    newName = '-'.join(newName)
    return newName


#-------
# LOAD DATA
#-------    
alignFile = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/MER20.dataset"
consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"


alignDF = pd.read_table(alignFile, sep="\t", header=0)
allCons_fasta = SeqIO.index(consensusFasta_file, "fasta", key_function = consensus_renameKey)


#-------
# MAIN
#-------
element=alignDF.elem[0]

#get consensus  sequences  
con_seq = str(allCons_fasta[element].seq).upper()

#create array to store alignment data 
startEnd = alignDF[['cons_start','cons_end']].values
hits = np.zeros( len(con_seq), dtype=np.int)

for x in startEnd: 
    hits[x[0]:x[1]] +=1


#plot data 
bases = np.arange(len(con_seq))

# Data for plotting
# Note that using plt.subplots below is equivalent to using
# fig = plt.figure and then ax = fig.add_subplot(111)
fig, ax = plt.subplots()
ax.plot(bases, hits)

ax.set(xlabel='bases in consensus', ylabel='number of hits',
       title='MER20 alignment to consensus')
ax.grid()

#fig.savefig("MER20_alignment.png")
plt.show()
