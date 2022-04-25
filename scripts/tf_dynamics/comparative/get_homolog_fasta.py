#!/bin/python
# This script will ...
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

WORKING_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/comparative"
BED_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/MER20_speciesFiltered.tsv"
OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files"
hg19_FASTA_FILE="/dors/capra_lab/abraha1/data/dna/hg19.fa"
mm9_FASTA_FILE="/dors/capra_lab/data/dna/mm9/mm9.fa"


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

# store homologs sequences between hg19 and mm9 
store_human = []
store_mouse = []

with open(BED_FILE, 'r') as f:
    for line in f:
        strip_line = line.rstrip().split('\t')

        try:
            if strip_line[5] == 'mm9': 
                store_human.append(strip_line[1:4])
                store_mouse.append(strip_line[6:9])
            elif strip_line[10] == 'mm9':
                store_human.append(strip_line[1:4])
                store_mouse.append(strip_line[11:14])
        except IndexError as err:
                pass

# filter out sequences with less than 50 bases 
filter_store_human = []
filter_store_mouse = []

for this_index, this_seq in enumerate(store_human): 
    if int(this_seq[2]) - int(this_seq[1]) >= 50: 
        filter_store_human.append(this_seq)
        filter_store_mouse.append(store_mouse[this_index])
        

# write coord to bed file 
with open(os.path.join(WORKING_DIR,'temp_human_coord.bed'), 'w') as temp_fh: 
    temp_fh.write(create_coord_string(filter_store_human))

with open(os.path.join(WORKING_DIR,'temp_mouse_coord.bed'), 'w') as temp_fh: 
    temp_fh.write(create_coord_string(filter_store_mouse))



# run bedtools getFasta 
cmd_hg19 = "bedtools getfasta -fi {} -bed {} -fo {}".format(hg19_FASTA_FILE, os.path.join(WORKING_DIR,'temp_human_coord.bed'), os.path.join(OUTPUT_DIR,'hg19_MER20.fa')).split()
cmd_mm9 = "bedtools getfasta -fi {} -bed {} -fo {}".format(mm9_FASTA_FILE, os.path.join(WORKING_DIR,'temp_mouse_coord.bed'), os.path.join(OUTPUT_DIR,'mm9_MER20.fa')).split()

check_hg19 = subprocess.check_output(cmd_hg19, shell=False, universal_newlines=True)
# check_mm9 = subprocess.check_output(cmd_mm9, shell=False, universal_newlines=True)

# clean up 
os.remove(os.path.join(WORKING_DIR,'temp_human_coord.bed'))
os.remove(os.path.join(WORKING_DIR,'temp_mouse_coord.bed'))