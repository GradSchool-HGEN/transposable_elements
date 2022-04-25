#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-03-03 19:37:44

from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
from Bio import SeqIO


# -----------
# PATH
# ----------- 


MAF_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/MER20_speciesFiltered.maf"
OUT_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/MER20_speciesFiltered.tsv"

# -----------
# MAIN
# ----------- 

out_file = open(OUT_FILE, 'w')

record_iterator = AlignIO.parse(MAF_FILE, "maf")

count = 0 
for this_paragraph in record_iterator: 
    count += 1
    store_homolog_line = []
    
    for this_row in this_paragraph: 

        ### this is where the parsing happens....
        this_species = this_row.id.split('.')[0]
        this_chr = this_row.id.split('.')[1]
        this_start = this_row.annotations['start']
        this_end =  this_row.annotations['start'] + this_row.annotations['size'] + 1 
        if this_row.annotations['strand'] == 1:
            this_strand = "+" 
        elif this_row.annotations['strand'] == -1: 
            this_strand = "-"
        else: 
            raise ValueError('strand parsing did not work')     
    
        store_homolog_line.append([this_species, this_chr, this_start, this_end, this_strand])
    
    store_homolog_line = [value for v in store_homolog_line for value in v]
    out_file.write('\t'.join(map(str,store_homolog_line)) + '\n')

    # if count == 2: 
    #     break 

out_file.close()
