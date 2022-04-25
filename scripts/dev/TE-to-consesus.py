#!/bin/python
# This script will map each instance of a transposable element in the hg38 genome to its corresponding consesus sequence. 
#
#
#
# Abin Abraham
# created on: 2017-11-30 14:39:44


import pandas as pd


#-------
# load data
TEPath="/dors/capra_lab/abraha1/projects/transposable_elements/data/unique_TE-list.txt"
consensusPath="/dors/capra_lab/abraha1/projects/transposable_elements/data/TE_cons_seqLength.out"

uniqTElist = pd.read_csv(TEPath, sep="\t", header=None)
con_df =  pd.read_csv(consensusPath, sep="\t", header=None)

conTElist = con_df[0]
conTElist = conTElist.str[:-10].tolist()
uniqTElist = uniqTElist[0]
uniqTElist = uniqTElist.tolist()

uniqTElist = [i[0] for i in uniqTElist]

#-------
matchList = []

for TE in uniqTElist:
    for consensus in conTElist:
        if TE == consensus:
            print(TE, consensus)
            matchList.append(TE)
            



#-------