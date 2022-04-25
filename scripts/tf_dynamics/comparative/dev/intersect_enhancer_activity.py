#!/bin/python
# This script intersect Villar enhancer activity with MER20 coordinates across Hsap, Macaques, and Mouse. 
#      
#
#
# Abin Abraham
# created on: 2018-03-05 09:13:42

import sys 
sys.path.append("/dors/capra_lab/abraha1/bin/global_scripts/py_scripts")

import subprocess 
import pandas as pd 
import numpy as np 
import tempfile
import os 

from  gen_utils import calc_overlap

WORKING_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/" 
MER20_MAPPED_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/MER20_speciesFiltered.tsv"
ENHANCER_HITS_DIR= "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/enhancer_activity"

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

#-------
# main  
#-------
species = ['hg19', 'rheMac2', 'mm9']
species_name = ['Hsap','Mmul','Mmus']

### create list of coordinates by species 
hg19_coords = list()
rheMac2_coords = list()
mm9_coords = list()
with open(MER20_MAPPED_FILE, 'r') as fh: 
    for line in fh: 
        strip_line = line.rstrip()
        hg19_coords.append(strip_line.split('\t')[1:5])

        try:
            if strip_line.split('\t')[5] == 'rheMac2': 
                rheMac2_coords.append(strip_line.split('\t')[6:10])
            elif strip_line.split('\t')[5] == 'mm9':
                mm9_coords.append(strip_line.split('\t')[6:10])

        except IndexError as err:
            pass # if this col does not exists, just ignore it and the index error 

        try:
            if strip_line.split('\t')[10] == 'rheMac2': 
                rheMac2_coords.append(strip_line.split('\t')[11:15])
            elif strip_line.split('\t')[10] == 'mm9':
                mm9_coords.append(strip_line.split('\t')[11:15])
        except IndexError as err:
            pass # if this col does not exists, just ignore it and the index error 

species_to_coord = {'hg19':hg19_coords, 'rheMac2':rheMac2_coords, 'mm9':mm9_coords}
species_to_species_name = {'hg19':'HSap', 'rheMac2':'Mmul', 'mm9':'Mmus'}

### run bedtools merge then bedtools intersect by species 
for this_species in species: 
    
    # create temp file to hold coordinates for a given species 
    temp_species_file = os.path.join(WORKING_DIR,'temp_coord_{}.bed'.format(this_species))
    with open(temp_species_file, 'w') as temp_fh:
        temp_fh.write(create_coord_string(species_to_coord[this_species]))

    # sort first 
    cmd_sort = "sort -k1,1 -k2,2n {} -o {}".format(temp_species_file,temp_species_file).split()
    sort_output = subprocess.check_output(cmd_sort, shell=False, universal_newlines=True)
    
    # bedtools merge 
    cmd_merge = "bedtools merge -i {}".format(temp_species_file).split()
    merge_output = subprocess.check_output(cmd_merge, shell=False, universal_newlines=True)

    temp_merged_species_file = os.path.join(WORKING_DIR,'temp_merged_coord_{}.bed'.format(this_species))
    with open(temp_merged_species_file, 'w') as temp_merged: 
        temp_merged.write(merge_output)

    # intersect hg19 with enhancer activity 
    enhancer_file =  os.path.join(ENHANCER_HITS_DIR, "merged_{}_H3K4me3.H3K27Ac_overlap_H3K27Aconly".format(species_to_species_name[this_species]))
    cmd = 'bedtools intersect -loj -a {} -b {}'.format(temp_merged_species_file, enhancer_file).split()
    check = subprocess.check_output(cmd, shell=False, universal_newlines=True)

    with open(os.path.join(WORKING_DIR,'{}_MER20_intersect_enhancer.tsv'.format(this_species)), 'w') as fh: 
        fh.write(check)

    ### remove temp files  
    os.remove(temp_species_file)
    os.remove(temp_merged_species_file)


### calc overlap between TE and enhancer overlap 
for this_species in species: 
    df = pd.read_csv(os.path.join(WORKING_DIR,'{}_MER20_intersect_enhancer.tsv'.format(this_species)), sep='\t',header=None, usecols=np.arange(7))
    df.columns = ['MER20_chr', 'MER20_start', 'MER20_end', 'Villar_chr', 'Villar_start','Villar_end', 'enhancer_activity']
    df['overlap_MER20'] = 0 
    df['overlap_EnhancerActivity'] = 0


    for row_index in np.arange(df.shape[0]):
        
        if df.iloc[row_index, 3] != ".":
            # print(df.iloc[row_index,:]) 
            overlap = calc_overlap(df.iloc[row_index,1], df.iloc[row_index,2], df.iloc[row_index,4], df.iloc[row_index,5])
            
            df.iloc[row_index, 7] = overlap[0]
            df.iloc[row_index, 8] = overlap[1]
    output_file_name = os.path.join(WORKING_DIR,'intersect_enhancer','{}_MER20_intersect_enhancer_wOverlap.tsv'.format(this_species))
    df.to_csv(output_file_name, sep='\t', header=True, index=False)

    ### remove files 
    os.remove(os.path.join(WORKING_DIR,'{}_MER20_intersect_enhancer.tsv'.format(this_species)))



### output summary 
print("------------")
print("Output stored in {} as <species>_MER20_intersect_enhancers_wOverlap.tsv".format(os.path.join(WORKING_DIR,'intersect_enhancer')))
print("------------")


# bedtools intersect -loj -a "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/temp_coord.bed"   -b "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/enhancer_activity/HSap_H3K4me3.H3K27Ac_overlap_H3K27Aconly"
# intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', anno_FILE, '-b', 'stdin', '-sorted'], stdin = filterElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)  
