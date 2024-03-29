#!/bin/python
# This script will check the mapping that are generated by createMapDict.py.           
#   Will extract input transposable element (TE) sequence and the model sequence from nhmmer output file. 
#   Will extract aligned result for TE and alignd consensus. 
#   Both input sequences and aligned sequences are printed to output file.
# 
# Abin Abraham
# created on: 2017-12-25 09:04:49

import parse_nhmmerOutput as par
import createMapDict as cmd
import numpy as np
from Bio import SeqIO
import time
import sys

def checkMapping(element, outputPath, inputElemSeqfile, consensusSeqfile, nhmmerParsed, allElemMappedDict):
    elemSeq = SeqIO.parse(inputElemSeqfile, 'fasta') #input TE fasta sequences for mapping 
    consensusSeq = SeqIO.read(consensusSeqfile,'fasta') #consensus sequences 
    
    for oneSeq in elemSeq: #loops through each TE in the genome
        if allElemMappedDict.get(oneSeq.id) == None: 
            print("In check_createMapDict.py: the following sequence did not \nhave mapping information to its consensus: "+oneSeq.id)
            continue
        else:
            ElemMappedDict = allElemMappedDict.get(oneSeq.id) # load the mapping dictionary 

        elemseq = str()
        consseq = str()
        for k in ElemMappedDict.keys(): #go through each base, keys is element index
            if k > -1: 
                elemseq = elemseq + str(oneSeq.seq[int(k)])
            else:
                elemseq = elemseq + '-'

        for v in ElemMappedDict.values():#go through each base, values is consensus index
            if v > -1: 
                consseq = consseq + str(consensusSeq.seq[int(v)])
            else: 
                consseq = consseq + '.'

        with open(outputPath+'checkMapping_'+element+'_'+time.strftime("%m_%d_%Y")+'.txt', 'a') as infile2:
            # print(oneSeq.id)
            infile2.write(str(oneSeq.id))
            infile2.write("\n")
            infile2.write("inputElemSeq  : ")
            infile2.write(nhmmerParsed[oneSeq.id]['elemSeq'])
            infile2.write("\n")
            infile2.write("checkElemSeq  : ")
            infile2.write(''.join(elemseq))
            infile2.write("\n")
            infile2.write("inputModeltSeq: ")
            infile2.write(nhmmerParsed[oneSeq.id]['model_seq'])
            infile2.write("\n")
            infile2.write("checkModelSeq : ")
            infile2.write(''.join(consseq))
            infile2.write("\n")
            infile2.write("\n")

'''
# load MER20 in genome (n=100)
MER20fa = SeqIO.parse('MER20_sample.fa','fasta')
# load consensus MER20 fasta 
conMER20fa = SeqIO.read('hmmemit-ouput-MER20.fa','fasta')
output = par.parseTermOut("nhmmer-output-terminal_allMER20")

for i, oneSeq in enumerate(MER20fa):
    
 
    elemTOcons = cmd.createMapDict(oneSeq.id, output)

    elemseq = str()
    consseq = str()
    for k in  elemTOcons.keys():
        if k > -1: 
            elemseq = elemseq + str(oneSeq.seq[int(k)])
        else:
            elemseq = elemseq + '-'

    for v in  elemTOcons.values():
        if v > -1: 
            consseq = consseq + str(conMER20fa.seq[int(v)])
        else: 
            consseq = consseq + '.'


    with open('check.txt', 'a') as infile2:
        print(oneSeq.id)
        infile2.write(str(oneSeq.id))
        infile2.write("\n")
        infile2.write(output[oneSeq.id]['elemSeq'])
        infile2.write("\n")
        infile2.write(''.join(elemseq))
        infile2.write("\n")
        infile2.write(output[oneSeq.id]['model_seq'])
        infile2.write("\n")
        infile2.write(''.join(consseq))
        infile2.write("\n")
        infile2.write("\n")

'''
