#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2017-11-29


import numpy as np 
import pandas as pd 
from Bio import SeqIO
import subprocess
# -----------
# FUNCTIONS
# -----------

def renameKey_consensus(identifier):
    newName = identifier.split("-")[0:-1]
    newName = '-'.join(newName)
    return newName

def renameKey_TEfasta(identifier):
    newName = identifier.split("(")[0]
    return newName


# -----------
# Load Data 
# -----------
#Paths 
TEbed_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/filtered_formated_hg38-TE.tsv"
consensusFasta_file = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
hg38Fasta_file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/hg38-TE.fasta" 


#Load Data 
TEds=pd.read_table(TEbed_file, sep="\t", header=None)
TEds.columns = ["chr", "start", "end", "TE", "TEfamily", "strand", "start_consensus", "end_consesus", "left_consensus"]

consensus_dict = SeqIO.index(consensusFasta_file, "fasta" , key_function = renameKey_consensus)
TE_fasta = SeqIO.index(hg38_fasta_file1, "fasta", key_function= renameKey_TEfasta)



# -----------
# main
# -----------
# get TE coords 
element_of_interest = "Eulor11"
element_coords = TEds.loc[TEds.TE ==element_of_interest, ['chr','start','end']]
element_coords_list = element_coords.values.tolist()
element_coords_list = [[x[0]+ ":" + str(x[1]) + "-" + str(x[2])] for x in element_coords_list]
#get fasta 
elem_seq = TE_fasta[element_coords_list[1][0]].seq

#get consensus sequence fasta 
con_seq = consensus_dict[element_of_interest].seq



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
# create a subset of hg38.fasta, that has only teh fasta for TE elements, hopefully make it go faster...

#element_coords.to_csv("elem_coord.tsv" ,sep='\t', index=False, header=False)
#element_fasta_list = subprocess.check_output(['bedtools', 'getfasta', '-s', '-fi', hg38_fasta_file, '-bed', "elem_coord.tsv" ], universal_newlines=True)
#need to parse ^ ouput
