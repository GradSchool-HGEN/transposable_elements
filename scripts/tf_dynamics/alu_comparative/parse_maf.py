#!/bin/python
# This script will parse a maf file and output a tsv with homologs coordinates per row 
#
#
#
# Abin Abraham
# created on: 2018-03-03 19:37:44

from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
from Bio import SeqIO
import os 

# -----------
# PATH
# ----------- 

### input_file_name:output_file_name
FILE_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/concatenated_multiz46/"
# FILE_DICT = { "AluJo_multiz46way_speciesFiltered.maf":"AluJo_multiz46way_speciesFiltered.tsv",  
#             "AluSp_multiz46way_speciesFiltered.maf":"AluSp_multiz46way_speciesFiltered.tsv",  
#             "AluYa5_multiz46way_speciesFiltered.maf":"AluYa5_multiz46way_speciesFiltered.tsv",
#             "AluSc_multiz46way_speciesFiltered.maf":"AluSc_multiz46way_speciesFiltered.tsv",  
#             "AluSq_multiz46way_speciesFiltered.maf":"AluSq_multiz46way_speciesFiltered.tsv",  
#             "AluYb8_multiz46way_speciesFiltered.maf":"AluYb8_multiz46way_speciesFiltered.tsv",
#             "AluSg_multiz46way_speciesFiltered.maf":"AluSg_multiz46way_speciesFiltered.tsv",  
#             "AluSx_multiz46way_speciesFiltered.maf":"AluSx_multiz46way_speciesFiltered.tsv", 
#             "AluY_multiz46way_speciesFiltered.maf":"AluY_multiz46way_speciesFiltered.tsv"} 

FILE_DICT = { "AluJb_mutiz46way_speciesFiltered.maf":"AluJb_multiz46way_speciesFiltered.tsv"} 




# -----------
# MAIN
# ----------- 
for input_MAF_file, output_file in FILE_DICT.items(): 
    out_file = open(os.path.join(FILE_PATH, output_file), 'w')
    record_iterator = AlignIO.parse(os.path.join(FILE_PATH, input_MAF_file), "maf")
    
    print("parsing {}...".format(input_MAF_file))

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

print("DONE! Check {}".format(FILE_PATH))
