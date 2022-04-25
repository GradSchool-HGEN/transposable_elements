#!/bin/python
# Given a species and appropriate input files, the fasta sequence for human and homolog species will be retrived. 
#
#
#
# Abin Abraham
# created on: 2018-03-06 11:42:39

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

import subprocess 
import os 
import sys 

#-------
# USER MUST MODIFY
#-------

BED_FILE_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/concatenated_multiz46"
# BED_FILE_DICT = {"AluJo_multiz46way_speciesFiltered.tsv":"AluJo",
# "AluSc_multiz46way_speciesFiltered.tsv":"AluSc",
# "AluSg_multiz46way_speciesFiltered.tsv":"AluSg",
# "AluSp_multiz46way_speciesFiltered.tsv":"AluSp",
# "AluSq_multiz46way_speciesFiltered.tsv":"AluSq",
# "AluSx_multiz46way_speciesFiltered.tsv":"AluSx",
# "AluYa5_multiz46way_speciesFiltered.tsv":"AluYa5",
# "AluYb8_multiz46way_speciesFiltered.tsv":"AluYb8",
# "AluY_multiz46way_speciesFiltered.tsv":"AluY"}

BED_FILE_DICT = {"AluJb_multiz46way_speciesFiltered.tsv":"AluJb"}

## WHICH HOMOLOGOUS SPECIES FASTA SHOULD BE RETRIEVED? 
# homolog_species = 'mm9'
homolog_species = 'rheMac2'

## PATHS TO FASTA GENOME FILES 
hg19_FASTA_FILE="/dors/capra_lab/abraha1/data/dna/hg19.fa"
mm9_FASTA_FILE="/dors/capra_lab/data/dna/mm9/mm9.fa"
rheMac2_FASTA_FILE="/dors/capra_lab/data/dna/rheMac2/rheMac2.fa"

if homolog_species == 'rheMac2': 
    homolog_fasta_file = rheMac2_FASTA_FILE
elif homolog_species == 'mm9': 
    homolog_fasta_file = mm9_FASTA_FILE

## OUTPUT DIR where fasta files will be saved 
WORKING_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/comparative" #store temp files
OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/fasta"
    

#-------
# funcitons
#-------
def create_coord_string(coords):
    coord_string = str()

    for this_coord in coords: 
        temp_string = '\t'.join(this_coord)+'\n'
        coord_string = coord_string+temp_string

    # coord_string = coord_string[:-2] # remove last newline 

    return coord_string

# -----------
# main
# ----------- 


for input_file, output_header in BED_FILE_DICT.items(): 

    # store homologs sequences between hg19 and homolog seq as selected above 
    store_human_seq = []
    store_homolog_seq = []
    print("Getting homologs between hg19 and {} for file {}....".format(homolog_species, input_file))

    with open(os.path.join(BED_FILE_PATH, input_file), 'r') as f:
        for line in f:
            strip_line = line.rstrip().split('\t')

            try: 
                if strip_line[5] == homolog_species:
                    store_human_seq.append(strip_line[1:4]) 
                    store_homolog_seq.append(strip_line[6:9])

                elif strip_line[10] == homolog_species:
                    store_human_seq.append(strip_line[1:4]) 
                    store_homolog_seq.append(strip_line[11:14]) 
            
            except IndexError as err: # ignore IndexErrors because that means no homolog was found 
                    pass


    # filter out sequences with less than 50 bases 
    filter_store_human_seq = []
    filter_store_homolog_seq = []


    # filter out sequence with < 50 bps 
    for this_index, this_seq in enumerate(store_human_seq): 
        if int(this_seq[2]) - int(this_seq[1]) >= 50: 
            filter_store_human_seq.append(this_seq)
            filter_store_homolog_seq.append(store_homolog_seq[this_index])
            

    # write coord to bed file 
    with open(os.path.join(WORKING_DIR,'temp_human_coord.bed'), 'w') as temp_fh: 
        temp_fh.write(create_coord_string(filter_store_human_seq))

    with open(os.path.join(WORKING_DIR,'temp_homolog_coord.bed'), 'w') as temp_fh: 
        temp_fh.write(create_coord_string(filter_store_homolog_seq))



    # run bedtools getFasta 
    cmd_hg19 = "bedtools getfasta -fi {} -bed {} -fo {}".format(hg19_FASTA_FILE, os.path.join(WORKING_DIR,'temp_human_coord.bed'), os.path.join(OUTPUT_DIR,'hg19_{}.fa'.format(output_header))).split()
    cmd_homolog= "bedtools getfasta -fi {} -bed {} -fo {}".format(homolog_fasta_file, os.path.join(WORKING_DIR,'temp_homolog_coord.bed'), os.path.join(OUTPUT_DIR,'{}_{}.fa'.format(homolog_species, output_header))).split()

    check_hg19 = subprocess.check_output(cmd_hg19, shell=False, universal_newlines=True)
    check_homolog = subprocess.check_output(cmd_homolog, shell=False, universal_newlines=True)

    # clean up 
    os.remove(os.path.join(WORKING_DIR,'temp_human_coord.bed'))
    os.remove(os.path.join(WORKING_DIR,'temp_homolog_coord.bed'))

print("DONE!")
print("check {}".format(OUTPUT_DIR))